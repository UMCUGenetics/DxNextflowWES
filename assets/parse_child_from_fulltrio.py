#! /usr/bin/env python
import argparse
from check_kinship import parse_ped

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check kinship output based on ped file.')
    parser.add_argument('ped_file', type=argparse.FileType('r'), help='PED file')
    parser.add_argument('vcf_files', nargs='+', help='VCF files (space seperated)')
    arguments = parser.parse_args()

    samples = parse_ped(arguments.ped_file)
    trio_sample = []
    for sample in samples:
        if len(samples[sample]['parents']) == 2:
            for vcf in list(set(arguments.vcf_files)):
                if sample in vcf:
                    trio_sample.append(sample)

    print(",".join(trio_sample))
