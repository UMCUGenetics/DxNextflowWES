#! /usr/bin/env python
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse flagstat files and print summary')
    parser.add_argument('flagstat_files', type=argparse.FileType('r'), nargs='*', help='Flagstat file')
    arguments = parser.parse_args()

    counts = {
        'total': 0,
        'mapped': 0,
        'paired': 0,
        'properly_paired': 0,
        'dups': 0,
        'files': len(arguments.flagstat_files)
    }

    for flagstat_file in arguments.flagstat_files:
        print("wokring on {0}...".format(flagstat_file.name))

        for line in flagstat_file:
            line = line.strip()
            print("\t{0}".format(line))

            line_count = int(line.split()[0])

            if 'in total' in line:
                counts['total'] += line_count

            elif 'paired in sequencing' in line:
                counts['paired'] += line_count

            elif 'properly paired' in line:
                counts['properly_paired'] += line_count

            elif 'duplicates' in line:
                sample_dups = float(line_count)
                counts['dups'] += line_count

            elif 'mapped (' in line:
                sample_mapped = float(line_count)
                counts['mapped'] += line_count

        print("\n\t{0} %duplication\n".format(100*sample_dups/sample_mapped))

    print("Total raw reads: {total:,} reads (Total throughput, 75bp={total_75bp:,} bp, 100bp={total_100bp:,} bp, 150bp={total_150bp:,} bp)".format(
        total=counts['total'], total_75bp=counts['total']*75, total_100bp=counts['total']*100, total_150bp=counts['total']*150
    ))
    print("Total mapped reads: {total:,} reads (Total throughput, 75bp={total_75bp:,} bp, 100bp={total_100bp:,} bp, 150bp={total_150bp:,} bp)".format(
        total=counts['mapped'], total_75bp=counts['mapped']*75, total_100bp=counts['mapped']*100, total_150bp=counts['mapped']*150
    ))
    print("Average mapped per lib: {:,} reads".format(int(round(float(counts['mapped'])/float(counts['files'])))))
    print("Average dups per lib: {:,} reads".format(int(round(float(counts['dups'])/float(counts['files'])))))
    print("Average dups % per lib: {:.2f} %".format(100*float(counts['dups'])/float(counts['mapped'])))
