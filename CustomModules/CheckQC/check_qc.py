#! /usr/bin/env python
# import statements, alphabetic order of main package.
import argparse
from pathlib import Path
import re
import sys

# third party libraries alphabetic order of main package.
from pandas import DataFrame, merge, read_csv
import yaml

# custom libraries alphabetic order
from utils import non_empty_existing_path


def parse_arguments_and_check(args_in):
    parser = argparse.ArgumentParser(
        description="Check and summarize sample quality using qc metrics and their thresholds."
    )
    parser.add_argument(
        "settings",
        type=non_empty_existing_path,
        help="QC settings specifying at least the filename, which qc_col, the threshold and the operator."
    )
    parser.add_argument(
        "output_path",
        type=non_empty_existing_path,
        default=".",
        help="QC check output path for all output files."
    )
    parser.add_argument(
        "output_prefix",
        type=str,
        default="qc_check",
        help="QC check output prefix for all output file names."
    )
    parser.add_argument(
        "input_files", nargs="+", type=non_empty_existing_path,
        help="Input files containing a QC metric."
    )
    args = parser.parse_args(args_in)
    return args


def read_yaml(yaml_file):
    yaml_loaded = yaml.safe_load(open(yaml_file))
    if not yaml_loaded:
        raise ValueError("File is empty.")
    return(yaml_loaded)


def check_allowed_operators(qc_operator):
    operators = ["<", "<=", ">", ">=", "==", "!=", "match"]
    if qc_operator not in operators:
        raise ValueError(f"Unsupported operator provided: {qc_operator}. Please select from: {operators}")


def check_required_keys_metrics(qc_settings):
    for req_key in ["filename", "qc_col", "threshold", "operator", "report_cols"]:
        if any([req_key not in setting.keys() for setting in qc_settings["metrics"]]):
            raise KeyError(f"Required key {req_key} not in all metrics settings.")


def select_metrics(filename, input_files):
    metrics = list(filter(re.compile(f".*{filename}").match, input_files))
    if not metrics:
        raise ValueError(f"No input file provided with filename pattern {filename}")
    return metrics


def get_columns_to_report(qc_report_cols, qc_metric_cols, qc_col):
    not_existing_cols = list(set(qc_report_cols) - set(qc_metric_cols))
    if qc_report_cols == "@all":
        qc_report_cols = qc_metric_cols
    elif not_existing_cols:
        raise ValueError(f"Some column names provided as report_cols do not exists: {not_existing_cols}")
    if not isinstance(qc_report_cols, str) and not isinstance(qc_report_cols, list):
        raise TypeError(f"{qc_report_cols} not string, list or '@all'")
    qc_report_cols = list(map(lambda x: x.replace(qc_col, "qc_value"), qc_report_cols))
    qc_report_cols.insert(0, "qc_title")
    return qc_report_cols


def add_and_rename_columns(qc_metric, qc_title, qc_col, qc_operator, qc_threshold):
    qc_metric_assigned = qc_metric.assign(
        qc_title=qc_title.lower(),
        qc_status="PASS",
        qc_check=f"{qc_threshold} {qc_operator} {qc_col}",
        qc_msg="",
    )
    qc_metric_out = qc_metric_assigned.rename(columns={qc_col: "qc_value"})
    return qc_metric_out


def get_failed_rows(qc_metric, qc_col, qc_operator, qc_threshold):
    # select failed rows using qc_threshold regex pattern and qc_operator 'match'
    if qc_operator == "match" and isinstance(qc_threshold, str):
        return qc_metric["qc_value"].str.match(qc_threshold)
    elif isinstance(qc_threshold, str):  # add quotes areound qc_threshold if it is string.
        # Note: using `query` has the advantage to dynamically build the comparison condition.
        # Disadvantage: no boolean indexing available. Assumed it is ok to use 'index'
        return qc_metric.query(f"`{qc_col}` {qc_operator} '{qc_threshold}'").index
    elif isinstance(qc_threshold, int) or isinstance(qc_threshold, float):  # query failed_rows using integers/floats
        # Note: using `query` has the advantage to dynamically build the comparison condition.
        # Disadvantage: no boolean indexing available. Assumed it is ok to use 'index'
        return qc_metric.query(f"`{qc_col}` {qc_operator} {qc_threshold}").index
    else:
        raise TypeError(f"QC threshold {qc_threshold} type not supported.")


def add_failed_samples_metric(qc_metric, failed_rows, report_cols, sample_cols):
    qc_metric_out = DataFrame(columns=["sample", "qc_check", "qc_status", "qc_msg", "qc_value"])
    failed_samples = []
    if failed_rows.to_list():
        failed_samples = list(qc_metric.loc[failed_rows, sample_cols].values.ravel())
        qc_metric.loc[failed_rows, "qc_status"] = "FAIL"
        # Concatenate columns to report into a single column 'qc_msg'
        qc_metric["qc_msg"] = qc_metric.loc[failed_rows, report_cols].astype(str).apply(" ".join, axis=1)
        # Add failed samples to output
        # A single qc metric could have multiple sample columns
        # If a qc check fails for a 'multiple sample check', each individual sample is flagged as "failed"
        for sample_col in sample_cols:
            qc_metric_out = qc_metric_out.append(
                (
                    qc_metric
                    .rename(columns={sample_col: "sample"})
                    .loc[failed_rows, qc_metric_out.columns.to_list()]
                    .groupby(["sample", "qc_check", "qc_status"])
                    .agg(lambda val: ';'.join(val.astype(str)))  # or .agg(lambda val: val.to_list())
                    .reset_index()
                )
            )
        # Drop failed samples current metric
        for sample_col in sample_cols:
            drop_index = qc_metric[qc_metric[sample_col].isin(set(failed_samples))].index
            if drop_index.to_list():
                qc_metric.drop(drop_index, inplace=True)
    return qc_metric, qc_metric_out


def add_passed_samples_metric(qc_metric, qc_metric_out, sample_cols):
    # Add passed samples to output
    for sample_col in sample_cols:
        qc_metric_out = qc_metric_out.append(
            (
                qc_metric
                .rename(columns={sample_col: "sample"})
                .loc[:, qc_metric_out.columns]
            )
        )
    # In case 'multiple sample qc check', 
    # output could contain duplicate rows for individual samples used in multiple comparisons.
    return qc_metric_out.sort_values(by=["qc_check", "qc_status"]).drop_duplicates(keep="first")


def create_and_write_output(qc_output, output_path, output_prefix):
    # Add qc_summary
    qc_output.insert(1, "qc_summary", "PASS")
    qc_output.loc[qc_output.isin(["FAIL"]).any(axis=1), "qc_summary"] = "FAIL"
    # Write summary output
    qc_output.to_csv(output_path + output_prefix + "_summary.csv", index=False, header=True)


def check_qc(input_files, settings, output_path, output_prefix):
    # A single qc metric file can be used multiple times, by defining a metric section for each check in the qc settings.
    qc_settings = read_yaml(settings)
    check_required_keys_metrics(qc_settings)
    for qc in qc_settings["metrics"]:
        check_allowed_operators(qc["operator"])
        metrics = select_metrics(qc["filename"], input_files)
        for qc_file in metrics:
            qc_metric_raw = read_csv(qc_file, comment=qc.get("comment", None), delimiter="\t", quotechar='"')
            report_cols = get_columns_to_report(qc["report_cols"], qc_metric_raw.columns.to_list(), qc["qc_col"])
            qc_metric_edit = add_and_rename_columns(qc_metric_raw, qc["title"], qc["qc_col"], qc["operator"], qc["threshold"])
            failed_rows = get_failed_rows(qc_metric_edit, "qc_value", qc["operator"], qc["threshold"])
            qc_metric_subset, qc_metric_judged = add_failed_samples_metric(
                qc_metric_edit, failed_rows, report_cols, qc["sample_cols"]
            )
            qc_metric_judged = add_passed_samples_metric(qc_metric_subset, qc_metric_judged, qc["sample_cols"])
            # Rename columns
            suffix = f"_{qc['title'].lower()}"
            qc_judged_renamed = qc_metric_judged.add_suffix(suffix).rename(columns={f"sample{suffix}": "sample"})
            # Concatenate/merge metric output
            try:
                output = merge(output, qc_judged_renamed, on="sample", how="outer")
            except NameError:  # first time:
                output = merge(
                    DataFrame(qc_metric_judged['sample'], columns=["sample"]),
                    qc_judged_renamed,
                    on="sample",
                    how="outer"
                )
    create_and_write_output(output, output_path, output_prefix)


if __name__ == "__main__":
    args = parse_arguments_and_check(args_in=sys.argv[1:])
    check_qc(args.input_files, args.settings, args.output_path, args.output_prefix)
