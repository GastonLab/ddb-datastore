#!/usr/bin/env python

import sys
import csv
import utils
import argparse
import getpass

from ddb import configuration
from collections import defaultdict


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

    with open("merkel_sequenced_bases_passing_depth.txt", 'w') as output:
        output.write("Sample\tMean Coverage\tKB passing\tTotal KB\tFraction")
        for sample in samples:
            num_bases = 0.0
            total_bases = 0.0
            total_coverage = 0.0
            with open("{}.recalibrated.sorted.bam.bed".format(sample), 'r') as coverage_file:
                sys.stdout.write("Reading base-by-base coverage data for {}\n".format(sample))
                reader = csv.reader(coverage_file, dialect='excel-tab')
                reader.next()
                for row in reader:
                    total_bases += 1
                    total_coverage += int(row[2])
                    if int(row[2]) >= int(args.min_depth):
                        num_bases += 1.0
                mean_coverage = total_coverage / total_bases
                output.write("{}\t{}\t{}\t{}\n".format(sample, mean_coverage, num_bases / 1000, total_bases / 1000,
                                                       num_bases / total_bases))
