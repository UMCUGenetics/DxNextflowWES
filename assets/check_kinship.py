#! /usr/bin/env python
import argparse


def parse_ped(ped_file):
    samples = {}  # 'sample_id': {'family': 'fam_id', 'parents': ['sample_id', 'sample_id']}

    for line in ped_file:
        ped_data = line.strip().split()
        family, sample, father, mother, sex, phenotype = ped_data

        # Create samples
        if sample not in samples:
            samples[sample] = {'family': family, 'parents': [], 'children': []}
        if father != '0' and father not in samples:
            samples[father] = {'family': family, 'parents': [], 'children': []}
        if mother != '0' and mother not in samples:
            samples[mother] = {'family': family, 'parents': [], 'children': []}

        # Save sample relations
        if father != '0':
            samples[sample]['parents'].append(father)
            samples[father]['children'].append(sample)
        if mother != '0':
            samples[sample]['parents'].append(mother)
            samples[mother]['children'].append(sample)
    return samples


def check_kinship(kinship_file, samples, kinship_setting):
    kinship_errors = False
    print_kinship('sample_1', 'sample_2', 'kinship', 'related', 'type', 'status')  # header
    for line in kinship_file:
        # Parse kinship data
        if line.startswith('FID1'):
            continue  # skip header line
        kinship_data = line.strip().split()
        sample_1, sample_2, kinship = kinship_data[1], kinship_data[3], float(kinship_data[7])

        # Check kinship data
        # Related
        if samples[sample_1]['family'] == samples[sample_2]['family']:
            # Parent - child
            if sample_2 in samples[sample_1]['parents'] or sample_1 in samples[sample_2]['parents']:
                if kinship > kinship_setting[0] and kinship < kinship_setting[1]:
                    print_kinship(sample_1, sample_2, kinship, True, 'parent_child', 'OK')
                else:
                    print_kinship(sample_1, sample_2, kinship, True, 'parent_child', 'FAIL')
                    kinship_errors = True
            # Parent - Parent -> both samples have the same children
            elif samples[sample_1]['children'] and samples[sample_1]['children'] == samples[sample_2]['children']:
                if kinship <= kinship_setting[0]:
                    print_kinship(sample_1, sample_2, kinship, True, 'parent_parent', 'OK')
                else:
                    print_kinship(sample_1, sample_2, kinship, True, 'parent_parent', 'FAIL')
                    kinship_errors = True
            # Aassume siblings
            else:
                if kinship > kinship_setting[0] and kinship < kinship_setting[1]:
                    print_kinship(sample_1, sample_2, kinship, True, 'sibling_sibling', 'OK')
                else:
                    print_kinship(sample_1, sample_2, kinship, True, 'sibling_sibling', 'FAIL')
                    kinship_errors = True
        # Unrelated
        else:
            if kinship <= kinship_setting[0]:
                print_kinship(sample_1, sample_2, kinship, False, 'NA', 'OK')
            else:
                print_kinship(sample_1, sample_2, kinship, False, 'NA', 'FAIL')
                kinship_errors = True

    return kinship_errors


def print_kinship(sample_1, sample_2, kinship, fam_status, relation_status, kinship_status):
    print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(sample_1, sample_2, kinship, fam_status, relation_status, kinship_status))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check kinship output based on ped file.')
    parser.add_argument('kinship_file', type=argparse.FileType('r'), help='Kinship file')
    parser.add_argument('ped_file', type=argparse.FileType('r'), help='PED file')
    arguments = parser.parse_args()

    # settings
    kinship_setting = [0.177, 0.354]

    # Parse ped file and check kinship
    samples = parse_ped(arguments.ped_file)
    kinship_errors = check_kinship(arguments.kinship_file, samples, kinship_setting)

    # Print summary
    if kinship_errors:
        print("\n# WARNING: Kinship errors found.")
    else:
        print("\n# No kinship errors found.")
    print("# Used kinship check settings: {0}".format(kinship_setting))
