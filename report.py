#!/usr/bin/env python

import sys
import getpass
import argparse
import toil_reporting_utils

from toil.job import Job
from ddb import configuration
from ddb_ngsflow import pipeline
from cassandra.auth import PlainTextAuthProvider


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
        sample_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)
        coverage_job = Job.wrapJobFn(toil_reporting_utils.get_coverage_data, samples, sample, [args.address],
                                     auth_provider, cores=1)

        # var_job = Job.wrapJobFn(toil_reporting_utils.get_variants, config, samples, sample, thresholds,
        #                         coverage_job.rv(), connection, [args.address], auth_provider, cores=1)
        #
        # report_job = Job.wrapJobFn(toil_reporting_utils.create_report, var_job.rv(), sample, samples, callers,
        #                            thresholds)

        root_job.addChild(sample_job)
        sample_job.addChild(coverage_job)
        # coverage_job.addChild(var_job)
        # sample_job.addFollowOn(report_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
