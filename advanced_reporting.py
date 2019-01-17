#!/usr/bin/env python

import sys
import xlwt
import utils
import getpass
import argparse

import numpy as np

from scipy import stats
from toil.job import Job
from ddb import configuration
from ddb_ngsflow import pipeline
from variantstore import Variant
from collections import defaultdict
from variantstore import SampleVariant
from coveragestore import SampleCoverage
from coveragestore import AmpliconCoverage

from cassandra.cluster import Cluster
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


def get_sample_coverage_data(sample, thresholds, authenticator){
    cluster = Cluster(['142.239.155.181', '142.239.155.182', '142.239.155.183',
                       '142.239.155.184'], auth_provider=authenticator)
    session = cluster.connect('coveragestore')
    rows = session.execute('SELECT sample, amplicon, run_id, library_name,
                           program_name, panel, num_reads, mean_coverage
                           FROM SampleCoverage')
    for amplicon_row in rows:
        print amplicon_row.sample, amplicon_row.amplicon, amplicon_row.num_reads
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file',
                        help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration',
                        help="Configuration file for various settings")
    parser.add_argument('-r', '--report',
                        help="Root name for reports (per sample)",
                        default='report')
    parser.add_argument('-a', '--address',
                        help="IP Address for Cassandra connection",
                        default='127.0.0.1')
    parser.add_argument('-u', '--username',
                        help='Cassandra username for login',
                        default=None)
    parser.add_argument('-d', '--min_depth',
                        help='Minimum depth threshold for variant reporting',
                        default=200.0)
    parser.add_argument('-g', '--good_depth',
                        help='Floor for good depth of coverage',
                        default=500.0)
    parser.add_argument('-t', '--min_somatic_var_freq',
                        help='Minimum reportable somatic variant frequency',
                        default=0.01)
    parser.add_argument('-p', '--max_pop_freq',
                        help='Maximum allowed population allele frequency',
                        default=0.005)

    # Job.Runner.addToilOptions(parser)
    # args = parser.parse_args()
    # args.logLevel = "INFO"

    config = configuration.configure_runtime(args.configuration)
    libraries = configuration.configure_samples(args.samples_file, config)
    samples = configuration.merge_library_configs_samples(libraries)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username,
                                              password=password)
    else:
        auth_provider = None

    thresholds = {'min_saf': args.min_somatic_var_freq,
                  'max_maf': args.max_pop_freq,
                  'depth': args.min_depth}

    callers = ("mutect", "platypus", "vardict", "scalpel", "freebayes",
               "pindel")

    sys.stdout.write("Processing samples\n")
    # root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)
    # amplicons_list_job = Job.wrapJobFn(get_all_amplicons, samples)
    # spawn_samples_job = Job.wrapJobFn(pipeline.spawn_variant_jobs)

    # root_job.addChild(amplicons_list_job)

    for sample in samples:
        get_sample_coverage_data(sample, thresholds)

    # Start workflow execution
    # Job.Runner.startToil(root_job, args)
