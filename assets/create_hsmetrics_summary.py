#! /usr/bin/env python
import argparse
import re
import sys

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

    sys.stdout.write('\t{}'.format('\t'.join(summary_data.keys())))
    sys.stdout.write('\n')
    for index, item in enumerate(summary_header):
        sys.stdout.write(item)
        sys.stdout.write('\t')
        for sample in summary_data:
            try:
                sys.stdout.write(summary_data[sample][index])
            except IndexError:
                sys.stdout.write('')
            sys.stdout.write('\t')
        sys.stdout.write('\n')
