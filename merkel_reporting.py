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

    parser.add_argument('-d', '--min_depth', help='Minimum depth threshold for variant reporting', default=100.0)
    parser.add_argument('-t', '--min_somatic_var_freq', help='Minimum reportable somatic variant frequency',
                        default=0.10)
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
    project_variant_data = defaultdict(lambda: defaultdict(int))
    variant_count_data = defaultdict(lambda: defaultdict(int))

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Processing variants for sample {}\n".format(sample))
        target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

        report_names = {'log': "{}.{}.log".format(sample, args.report),
                        'coverage': "{}_coverage_{}.txt".format(sample, args.report),
                        'tier1_pass': "{}_tier1_pass_variants_{}.txt".format(sample, args.report),
                        'tier1_fail': "{}_tier1_fail_variants_{}.txt".format(sample, args.report),
                        'vus_pass': "{}_vus_pass_variants_{}.txt".format(sample, args.report),
                        'vus_fail': "{}_vus_fail_variants_{}.txt".format(sample, args.report),
                        'tier4_pass': "{}_tier4_pass_variants_{}.txt".format(sample, args.report),
                        'tier4_fail': "{}_tier4_fail_variants_{}.txt".format(sample, args.report),
                        'all_ordered': "{}_all_ordered_variants_{}.txt".format(sample, args.report),
                        'categories': "{}_variants_by_category{}.txt".format(sample, args.report)
                        }

        with open(report_names['log'], 'w') as logfile:
            logfile.write("Reporting Log for sample {}\n".format(sample))
            logfile.write("---------------------------------------------\n")

        with open(report_names['coverage'], "w") as coverage_report:
            coverage_report.write("Sample:\t{}\n".format(sample))
            coverage_report.write("---------------------------------------------\n")

        utils.setup_report_header(report_names['tier1_pass'], callers)
        utils.setup_report_header(report_names['tier1_fail'], callers)

        utils.setup_report_header(report_names['vus_pass'], callers)
        utils.setup_report_header(report_names['vus_fail'], callers)

        utils.setup_report_header(report_names['tier4_pass'], callers)
        utils.setup_report_header(report_names['tier4_fail'], callers)

        utils.setup_report_header(report_names['all_ordered'], callers)

        for library in samples[sample]:
            report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}" \
                                "".format(samples[sample][library]['panel'], samples[sample][library]['report'])
            target_amplicons = utils.get_target_amplicons(report_panel_path)

            with open(report_names['log'], 'a') as logfile:
                sys.stdout.write("Processing amplicons for library {} from file {}\n".format(library,
                                                                                             report_panel_path))
                logfile.write("Processing amplicons for library {} from file {}\n".format(library, report_panel_path))

            ordered_variants, num_var = utils.get_variants(config, samples, sample, library, thresholds, report_names)

            sys.stdout.write("Processing amplicon coverage data\n")
            reportable_amplicons, target_amplicon_coverage = utils.get_coverage_data(target_amplicons, samples, sample,
                                                                                     library, target_amplicon_coverage)

            sys.stdout.write("Filtering and classifying variants\n")
            filtered_var_data = utils.classify_and_filter_variants_proj(samples, sample, library, report_names,
                                                                        target_amplicons, callers, ordered_variants,
                                                                        config, thresholds, project_variant_data,
                                                                        variant_count_data)
            project_variant_data = filtered_var_data[-2]
            variant_count_data = filtered_var_data[-1]

            sys.stdout.write("Writing variant reports\n")
            utils.write_reports(report_names, samples, sample, library, filtered_var_data, ordered_variants,
                                target_amplicon_coverage, reportable_amplicons, num_var, thresholds, callers)

    sys.stdout.write("Writing project/run level data\n")
    with open("Summary_Data.txt", 'w') as summary:
        summary.write("Variant\tNum Tier1 Pass\tNum Tier1 Fail\tNum VUS Pass\tNum VUS Fail\tNum Tier4 Pass\t"
                      "Num Tier4 Fail\n")
        for variant_id in project_variant_data:
            summary.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(variant_id,
                                                                project_variant_data[variant_id]['tier1_pass'],
                                                                project_variant_data[variant_id]['tier1_fail'],
                                                                project_variant_data[variant_id]['vus_pass'],
                                                                project_variant_data[variant_id]['vus_fail'],
                                                                project_variant_data[variant_id]['tier4_pass'],
                                                                project_variant_data[variant_id]['tier4_fail']))

    sys.stdout.write("Writing project/run level category data\n")
    with open("Category_Data.txt", 'w') as summary:
        summary.write("Variant\tNum Pos\tNum Neg\tDiff\n")
        for variant_id in project_variant_data:
            diff = abs(project_variant_data[variant_id]['positive'] - project_variant_data[variant_id]['negative'])
            summary.write("{}\t{}\t{}\t{}\n".format(variant_id, project_variant_data[variant_id]['positive'],
                                                    project_variant_data[variant_id]['negative'], diff))

    sys.stdout.write("Writing Sample-level variant count data\n")
    with open("Sample_Variant_Counts.txt", 'w') as summary:
        summary.write("Sample\tGroup\tViral Status\tNumber Passing Variants\tNumber Passing C>T Variants")
        for sample in variant_count_data:
            summary.write("{}\t{}\t{}\t{}\t{}\n".format(sample, samples[sample]['category'], samples[sample]['viral'],
                                                        variant_count_data[sample]['pass_count'],
                                                        variant_count_data[sample]['CT_count']))

