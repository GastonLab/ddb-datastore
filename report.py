#!/usr/bin/env python

import sys
import csv
import utils
import argparse
import getpass

from ddb import configuration
from collections import defaultdict

from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default='any_evidence_report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)

    parser.add_argument('-d', '--min_depth', help='Minimum depth threshold for variant reporting', default=200)
    parser.add_argument('-t', '--min_somatic_var_freq', help='Minimum reportable somatic variant frequency',
                        default=0.01)
    parser.add_argument('-p', '--max_pop_freq', help='Maximum allowed population allele frequency', default=0.005)

    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    libraries = configuration.configure_samples(args.samples_file, config)

    samples = configuration.merge_library_configs_samples(libraries)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    thresholds = {'min_saf': args.min_somatic_var_freq,
                  'max_maf': args.max_pop_freq,
                  'depth': args.min_depth}

    callers = ("mutect", "platypus", "vardict", "scalpel", "freebayes", "pindel")

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Processing variants for sample {}\n".format(sample))
        target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

        with open("{}.{}.log".format(sample, args.report), 'w') as logfile:
            logfile.write("Reporting Log for sample {}\n".format(sample))
            logfile.write("---------------------------------------------\n")

        with open("{}_coverage_{}.txt".format(sample, args.report), "w") as coverage_report:
            coverage_report.write("Sample:\t{}\n".format(sample))
            coverage_report.write("---------------------------------------------\n")

        utils.setup_report_header("{}_tier1_pass_variants_{}.txt".format(sample, args.report), callers)
        utils.setup_report_header("{}_tier1_fail_variants_{}.txt".format(sample, args.report), callers)

        utils.setup_report_header("{}_vus_pass_variants_{}.txt".format(sample, args.report), callers)
        utils.setup_report_header("{}_vus_fail_variants_{}.txt".format(sample, args.report), callers)

        utils.setup_report_header("{}_tier4_pass_variants_{}.txt".format(sample, args.report), callers)
        utils.setup_report_header("{}_tier4_fail_variants_{}.txt".format(sample, args.report), callers)

        for library in samples[sample]:
            report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}" \
                                "".format(samples[sample][library]['panel'], samples[sample][library]['report'])
            target_amplicons = utils.get_target_amplicons(report_panel_path)

            with open("{}.{}.log".format(sample, args.report), 'a') as logfile:
                sys.stdout.write("Processing amplicons for library {} from file {}\n".format(library,
                                                                                             report_panel_path))
                logfile.write("Processing amplicons for library {} from file {}\n".format(library, report_panel_path))

            ordered_variants = utils.get_variants(config, samples, sample, library, thresholds, args.report)
            reportable_amplicons, target_amplicon_coverage = utils.get_coverage_data(target_amplicons, samples, sample,
                                                                                     library, target_amplicon_coverage)

            filtered_var_data = utils.classify_and_filter_variants(sample, library, args.report, target_amplicons,
                                                                   callers, ordered_variants, thresholds)

            utils.write_reports(args.report, samples, sample, library, filtered_var_data, target_amplicon_coverage,
                                reportable_amplicons, thresholds, callers)
