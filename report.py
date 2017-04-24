#!/usr/bin/env python

import sys
import utils
import getpass
import argparse

from toil.job import Job
from ddb import configuration
from ddb_ngsflow import pipeline
from collections import defaultdict

from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


def process_reporting_sample(job, sample, samples, report_root, callers, config, thresholds, target_amplicon_coverage,
                             addresses, authenticator):
    connection.setup(addresses, "variantstore", auth_provider=authenticator)

    report_names = {'log': "{}.{}.log".format(sample, report_root),
                    'coverage': "{}_coverage_{}.txt".format(sample, report_root),
                    'tier1_pass': "{}_tier1_pass_variants_{}.txt".format(sample, report_root),
                    'tier1_fail': "{}_tier1_fail_variants_{}.txt".format(sample, report_root),
                    'vus_pass': "{}_vus_pass_variants_{}.txt".format(sample, report_root),
                    'vus_fail': "{}_vus_fail_variants_{}.txt".format(sample, report_root),
                    'tier4_pass': "{}_tier4_pass_variants_{}.txt".format(sample, report_root),
                    'tier4_fail': "{}_tier4_fail_variants_{}.txt".format(sample, report_root),
                    'all_ordered': "{}_all_ordered_variants_{}.txt".format(sample, report_root)
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
            job.fileStore.logToMaster("Processing amplicons for library {} from file {}\n".format(library,
                                                                                                  report_panel_path))
            logfile.write("Processing amplicons for library {} from file {}\n".format(library, report_panel_path))

        ordered_variants, num_var = utils.get_variants(config, samples, sample, library, thresholds, report_names)

        job.fileStore.logToMaster("Processing amplicon coverage data\n")
        reportable_amplicons, target_amplicon_coverage = utils.get_coverage_data(target_amplicons, samples, sample,
                                                                           library, target_amplicon_coverage)

        job.fileStore.logToMaster("Filtering and classifying variants\n")
        filtered_var_data = utils.classify_and_filter_variants(samples, sample, library, report_names, target_amplicons,
                                                         callers, ordered_variants, config, thresholds)

        job.fileStore.logToMaster("Writing variant reports\n")

        utils.write_reports(job, report_names, samples, sample, library, filtered_var_data, ordered_variants,
                      target_amplicon_coverage, reportable_amplicons, num_var, thresholds, callers)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default='report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-d', '--min_depth', help='Minimum depth threshold for variant reporting', default=250.0)
    parser.add_argument('-t', '--min_somatic_var_freq', help='Minimum reportable somatic variant frequency',
                        default=0.01)
    parser.add_argument('-p', '--max_pop_freq', help='Maximum allowed population allele frequency', default=0.005)
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    libraries = configuration.configure_samples(args.samples_file, config)

    samples = configuration.merge_library_configs_samples(libraries)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
    else:
        auth_provider = None

    thresholds = {'min_saf': args.min_somatic_var_freq,
                  'max_maf': args.max_pop_freq,
                  'depth': args.min_depth}

    callers = ("mutect", "platypus", "vardict", "scalpel", "freebayes", "pindel")

    sys.stdout.write("Processing samples\n")
    root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)
    for sample in samples:
        sys.stdout.write("Processing variants for sample {}\n".format(sample))
        target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

        process_job = Job.wrapJobFn(process_reporting_sample, sample, samples, args.report, callers, config,
                                    thresholds, target_amplicon_coverage, [args.address], auth_provider, cores=1)

        root_job.addChild(process_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
