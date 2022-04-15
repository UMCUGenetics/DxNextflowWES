#! /usr/bin/env python
import argparse
from check_kinship import parse_ped

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check kinship output based on ped file.')
    parser.add_argument('ped_file', type=argparse.FileType('r'), help='PED file')
    parser.add_argument('samples_analysis', nargs='+', help='samples within the analysis (space seperated)')
    arguments = parser.parse_args()

    samples = parse_ped(arguments.ped_file)
    trio_sample = []
    for sample in samples:
        if len(samples[sample]['parents']) == 2:  # Sample = child with parents
            for sample_analysis in list(set(arguments.samples_analysis)):
                if sample in sample_analysis:
                    trio_sample.append(sample)

    print(",".join(trio_sample))
