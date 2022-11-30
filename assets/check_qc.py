# import statements, alphabetic order of main package.
import argparse
from pathlib import Path
import re
import sys

# third party libraries alphabetic order of main package.
from pandas import DataFrame, read_csv
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
        default=str(Path(__file__).parent) + "/qc_settings.yaml",
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


def check_required_keys_metrics(qc_settings):
    for req_key in ["filename", "qc_col", "threshold", "operator", "report_cols"]:
        if any([req_key not in setting.keys() for setting in qc_settings['metrics']]):
            raise KeyError(f"Required key {req_key} not in all metrics settings.")


def select_metrics(filename, input_files):
    metrics = list(filter(re.compile(f".*{filename}").match, input_files))
    if not metrics:
        raise ValueError(f"No input file provided with filename pattern {filename}")
    return metrics


def retrieve_failed_qc(qc_metric, qc_col, qc_operator, qc_threshold, report_cols):
    if report_cols == "@all":
        report_cols = qc_metric.columns.to_list()
    if not isinstance(report_cols, str) and not isinstance(report_cols, list):
        raise ValueError(f"{report_cols} not string, list or '@all'")
    if isinstance(qc_threshold, str):
        qc_subset = qc_metric[report_cols].query(f"`{qc_col}` {qc_operator} '{qc_threshold}'")
    else:
        qc_subset = qc_metric[report_cols].query(f"`{qc_col}` {qc_operator} {qc_threshold}")
    return qc_subset


def check_qc(input_files, settings, output_file):
    qc_settings = read_yaml(settings)
    check_required_keys_metrics(qc_settings)
    with open(output_file, 'w') as f:
        f.write(",".join(["qc_title", "qc_check", "qc_input_raw"]) + "\n")
    for qc in qc_settings['metrics']:
        metrics = select_metrics(qc['filename'], input_files)
        for qc_file in metrics:  # "/Users/ejong19/Downloads/multiqc_picard_HsMetrics.txt"
            qc_metric = read_csv(qc_file, comment=qc.get("comment", None), delimiter="\t", quotechar='"')
            qc_metric_subset = retrieve_failed_qc(qc_metric, qc['qc_col'], qc['operator'], qc['threshold'], qc['report_cols'])
            if not qc_metric_subset.empty:
                qc_cols = qc_metric_subset.columns.tolist()
                # qc_metric_subset["qc_input"] = qc_metric_subset.astype(str).apply(' '.join, axis=1)
                qc_metric_subset["qc_title"] = qc['title']
                qc_metric_subset["qc_check"] = f"{qc['qc_col']} {qc['operator']} {qc['threshold']}"
                qc_metric_subset[["qc_title", "qc_check"] + qc_cols].to_csv(
                    output_file, mode='a', index=False, header=False
                )  # ["qc_title", "qc_check", "qc_input"]


if __name__ == '__main__':
    args = parse_arguments_and_check(args_in=sys.argv[1:])
    check_qc(args.input_files, args.settings, args.output_file)

# settings = "/Users/ejong19/repos/DxNextflowWES/assets/qc_settings.yaml"
# input_files = ["/Users/ejong19/Downloads/multiqc_picard_HsMetrics.txt", "/Users/ejong19/Downloads/multiqc_verifybamid.txt",
#                "/Users/ejong19/Downloads/221111_A00295_0680_BHGCM3DMXY_1.kinship_check.out"]
# output_file = '/Users/ejong19/Downloads/failed.csv'
