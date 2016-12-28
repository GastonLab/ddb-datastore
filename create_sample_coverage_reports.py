#!/usr/bin/env python

import argparse
import getpass
import sys
import csv

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

import utils
from coveragestore import SampleCoverage

from collections import defaultdict


def get_target_amplicons(filename):
    amplicons_list = list()
    sys.stdout.write("Opening file {} to retrieve reporting amplicons\n".format(filename))
    with open(filename, "r") as bedfile:
        reader = csv.reader(bedfile, dialect='excel-tab')
        for row in reader:
            amplicons_list.append(row[3])

    return amplicons_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)

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

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Processing coverage for sample {}\n".format(sample))
        report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}".format(samples[sample]['panel'],
                                                                                       samples[sample]['report'])
        target_amplicons = get_target_amplicons(report_panel_path)
        reportable_amplicons = list()
        for amplicon in target_amplicons:
            coverage_data = SampleCoverage.objects.timeout(None).filter(
                SampleCoverage.sample == samples[sample]['sample_name'],
                SampleCoverage.amplicon == amplicon,
                SampleCoverage.run_id == samples[sample]['run_id'],
                SampleCoverage.library_name == samples[sample]['library_name'],
                SampleCoverage.program_name == "sambamba"
            )
            ordered_variants = coverage_data.order_by('amplicon', 'run_id').limit(coverage_data.count() + 1000)
            for variant in ordered_variants:
                reportable_amplicons.append(variant)

        with open("{}_{}.txt".format(sample, args.report), "w") as coverage_report:
            coverage_report.write("Sample\tLibrary\tAmplicon\tNum Reads\tCoverage\n")
            for amplicon in reportable_amplicons:
                coverage_report.write("{}\t{}\t{}\n".format(amplicon.sample,
                                                            amplicon.library_name,
                                                            amplicon.amplicon,
                                                            amplicon.num_reads,
                                                            amplicon.mean_coverage))
