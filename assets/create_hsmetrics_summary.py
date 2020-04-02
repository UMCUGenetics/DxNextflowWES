#! /usr/bin/env python
import argparse
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create HSMetrics summary file')
    parser.add_argument('hsmetrics_files', type=argparse.FileType('r'), nargs='*', help='PED file')
    arguments = parser.parse_args()

    interval_files_pattern = re.compile("BAIT_INTERVALS=\[(\S*)\].TARGET_INTERVALS=\[(\S*)\]")
    summary_header = []
    summary_data = {}
    for hsmetrics_file in arguments.hsmetrics_files:
        sample = hsmetrics_file.name.split('.')[0]
        header = []
        bait_intervals = ''
        target_intervals = ''

        for line in hsmetrics_file:
            line = line.strip()
            if line.startswith('## HISTOGRAM'):  # Skip all histogram lines
                break

            elif interval_files_pattern.search(line):
                bait_intervals = interval_files_pattern.search(line).group(1)
                target_intervals = interval_files_pattern.search(line).group(2)

            elif line.startswith('BAIT_SET'):
                header = ['baitIntervals', 'targetIntervals']
                header += line.split('\t')
                if not summary_header:
                    summary_header = header
                elif header != summary_header:
                    print "ERROR"

            elif header and line:
                data = [bait_intervals, target_intervals]
                data += line.split('\t')
                summary_data[sample] = data

    with open('HSMetrics_summary.txt', 'w') as hsmetrics_summary:
        hsmetrics_summary.write('\t{}'.format('\t'.join(summary_data.keys())))
        hsmetrics_summary.write('\n')
        for index, item in enumerate(summary_header):
            hsmetrics_summary.write(item)
            hsmetrics_summary.write('\t')
            for sample in summary_data:
                try:
                    hsmetrics_summary.write(summary_data[sample][index])
                except IndexError:
                    hsmetrics_summary.write('')
                hsmetrics_summary.write('\t')
            hsmetrics_summary.write('\n')
