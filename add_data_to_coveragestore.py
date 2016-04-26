#!/usr/bin/env python

import re
import sys
import getpass
import argparse
import utils
import cyvcf2

from cyvcf2 import VCF
from datetime import datetime
from collections import defaultdict
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider

from coveragestore import AmpliconCoverage
from coveragestore import SampleCoverage
from ddb import configuration
from ddb_ngsflow import pipeline

from toil.job import Job


def process_sample_coverage(job, address, keyspace, auth, thresholds, report_root, sample, samples, config):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default=None)
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
    else:
        auth_provider = None

    thresholds = {'threshold1': 500,
                  'threshold2': 1000,
                  'regions': config['actionable_regions']}

    root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)

    for sample in samples:
        sample_job = Job.wrapJobFn(process_sample_coverage, [args.address], "coveragestore", auth_provider,
                                   thresholds, args.report, sample, samples, config,
                                   cores=1)
        root_job.addChild(sample_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)