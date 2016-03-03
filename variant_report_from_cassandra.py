#!/usr/bin/env python

import sys
import argparse
import utils
from ddb import configuration
from variantstore import Variant
from cassandra.cqlengine import connection

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    connection.setup(['127.0.0.1'], "variantstore")

    thresholds = {'max_aaf': 0.01}
    callers = ['mutect', 'freebayes', 'scalpel', 'vardict']

    for sample in samples:
        variants = Variant.objects(Variant.sample == samples[sample]['sample_name'])
        utils.write_sample_variant_report(args.report, sample, variants, callers, thresholds)
