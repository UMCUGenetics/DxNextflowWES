#! /usr/bin/env python
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create ExomeDepth Metric summary file')
    parser.add_argument('exomedepth_files', type=argparse.FileType('r'), nargs='*', help='Exomedepth files')
    arguments = parser.parse_args()

    print("Sample\tModel\tCorrelation\tDelDupRatio\tNumberOfCalls") 
    for hsmetrics_file in arguments.exomedepth_files:
        for line in hsmetrics_file: 
            print(line.rstrip())
