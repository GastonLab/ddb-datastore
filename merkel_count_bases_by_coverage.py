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
    parser.add_argument('-d', '--min_depth', help='Minimum depth threshold for variant reporting', default=250.0)

    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    with open("merkel_sequenced_bases_passing_depth.txt") as output:
        for sample in samples:
            num_bases = 0
            with open("{}.recalibrated.sorted.bam.bed".format(sample), 'r') as coverage_file:
                sys.stdout.write("Reading base-by-base coverage data for {}\n".format(sample))
                reader = csv.reader(coverage_file, dialect='excel-tab')
                reader.next()
                for row in reader:
                    if int(row[2]) >= args.min_depth:
                        num_bases += 1
                output.write("{}\t{}\n".format(sample, num_bases))
