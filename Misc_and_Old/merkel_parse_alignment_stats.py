#!/usr/bin/env python

import sys
import csv
import argparse

from ddb import configuration


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-d', '--min_depth', help='Minimum depth threshold for variant reporting', default=250)

    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    with open("merkel_alignment_stats.txt", 'w') as output:
        output.write("Sample\tNum Reads\tNum Aligned\tPct Aligned\tPct On-Target Bases\n")
        for sample in samples:
            with open("{}.capture_alignment_metrics".format(sample), 'r') as coverage_file:
                sys.stdout.write("Reading alignment metrics for {}\n".format(sample))
                reader = csv.reader(coverage_file, dialect='excel-tab')
                # Skip first 9 rows for a single sample to get rid of header info and first/second read stats.
                # Only read data for pair
                rows = list()
                for row in reader:
                    rows.append(row)

                row = rows[7]

                num_reads = row[5]
                aligned_reads = row[10]
                pct_aligned_reads = row[11]
                pct_selected_bases = row[18]
                output.write("{}\t{}\t{}\t{}\t{}\n"
                             "\n".format(sample, num_reads, aligned_reads, pct_aligned_reads, pct_selected_bases))
