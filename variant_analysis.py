#!/usr/bin/env python

import sys
import csv
import argparse
import utils
import getpass
import plotly
import plotly.plotly as py
import plotly.graph_objs as go

from collections import defaultdict

from ddb import configuration
from variantstore import SampleVariant
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-v', '--variant_callers', help="Comma-delimited list of variant callers used")

    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    thresholds = {'min_saf': 0.02,
                  'max_maf': 0.005,
                  'depth': 500,
                  'regions': config['actionable_regions']}

    callers = ['mutect', 'freebayes', 'scalpel', 'vardict', 'platypus', 'pindel']

    summary_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Running Cassandra query for library {}\n".format(samples[sample]['library_name']))
        variants = SampleVariant.objects.timeout(None).filter(SampleVariant.sample == samples[sample]['sample_name'],
                                                              SampleVariant.run_id == samples[sample]['run_id'],
                                                              SampleVariant.library_name == samples[sample]['library_name'],
                                                              SampleVariant.max_som_aaf >= thresholds['min_saf'],
                                                              SampleVariant.max_maf_all <= thresholds['max_maf'],
                                                              SampleVariant.max_depth >= thresholds['depth']
                                                              ).allow_filtering()

        ordered_variants = variants.order_by('reference_genome', 'chr', 'pos').limit(variants.count() + 1000)

        sys.stdout.write("Retrieved {} total variants\n".format(variants.count()))
        sys.stdout.write("Running filters on sample variants\n")
        passing_variants = list()
        passed = 0
        iterated = 0
        for variant in ordered_variants:
            iterated += 1
            flag, info = utils.variant_filter(variant, callers, thresholds)
            passing_variants.append((variant, flag, info))
            if info['clinvar'] == "Not Benign":
                summary_data["{}_{}".format(samples[sample]['sample_name'], samples[sample]['target_pool'])][samples[sample]['num_libraries_in_run']]['clinvar'] += 1
            if variant.severity != "LOW":
                summary_data["{}_{}".format(samples[sample]['sample_name'], samples[sample]['target_pool'])][samples[sample]['num_libraries_in_run']]['med_high'] += 1

        summary_data["{}_{}".format(samples[sample]['sample_name'], samples[sample]['target_pool'])][samples[sample]['num_libraries_in_run']]['total'] = iterated

    # Output
    with open("variant_report.txt", 'w') as output:
        for library in summary_data:
            output.write("{}".format(library))
