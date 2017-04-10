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
        with open("{}.{}.log".format(sample, args.report), 'w') as logfile:
            logfile.write("Reporting Log for sample {}\n".format(sample))
            logfile.write("---------------------------------------------\n")

        on_target_caller_variants = list()
        filtered_no_amplicon = list()
        filtered_non_target_amplicon = list()
        filtered_no_requested_caller = list()

        target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

        for library in samples[sample]:
            off_target_amplicons = defaultdict(int)

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

            for variant in ordered_variants:
                if variant.amplicon_data['amplicon']:
                    amplicons = variant.amplicon_data['amplicon'].split(',')
                    for amplicon in amplicons:
                        if amplicon in target_amplicons:
                            for caller in callers:
                                if caller in variant.callers:
                                    on_target_caller_variants.append(variant)
                                    break
                        else:
                            filtered_non_target_amplicon.append(variant)
                            off_target_amplicons[amplicon] += 1
                elif variant.amplicon_data['amplicon'] is 'None':
                    filtered_no_amplicon.append(variant)
                    off_target_amplicons[amplicon] += 1
                else:
                    filtered_no_amplicon.append(variant)
                    off_target_amplicons[amplicon] += 1

            with open("{}.{}.log".format(sample, args.report), 'a') as logfile:
                logfile.write("---------------------------------------------\n")
                logfile.write("{}\n".format(library))
                logfile.write(
                    "Sent {} variants to reporting (filtered {} variants for no amplicon data and {} for being"
                    " in a non-targeted amplicon)\n".format(len(on_target_caller_variants), len(filtered_no_amplicon),
                                                            len(filtered_non_target_amplicon)))
                logfile.write("---------------------------------------------\n")
                logfile.write("Off Target Amplicon\tCounts\n")
                for off_target in off_target_amplicons:
                    logfile.write("{}\t{}\n".format(off_target, off_target_amplicons[off_target]))

            sys.stdout.write("Sending {} variants to reporting (filtered {} variants for no amplicon data and {} for "
                             "being in a non-targeted amplicon)\n".format(len(on_target_caller_variants),
                                                                          len(filtered_no_amplicon),
                                                                          len(filtered_non_target_amplicon)))

        utils.write_sample_variant_report_no_caller_filter(args.report, sample, on_target_caller_variants,
                                                           target_amplicon_coverage, callers)
