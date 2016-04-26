#!/usr/bin/env python

import sys
import csv
import getpass
import argparse

from collections import defaultdict
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider

from coveragestore import AmpliconCoverage
from coveragestore import SampleCoverage
from ddb import configuration
from ddb_ngsflow import pipeline

from toil.job import Job


def process_sample_coverage(job, addresses, keyspace, auth, report_root, sample, program, samples, config):
    connection.setup(addresses, keyspace, auth_provider=auth)

    with open("{}.sambamba_coverage.bed".format(sample), 'rb') as coverage:
        reader = csv.reader(coverage, delimiter='\t')
        header = reader.next()
        for row in reader:

            sample_data = SampleCoverage.create(sample=sample,
                                                library_name=samples[sample]['library_name'],
                                                run_id=samples[sample]['run_id'],
                                                program_name=program,
                                                extraction=samples[sample]['extraction'],
                                                amplicon_name=row[3],
                                                amplicon_depth=row[4])

            amplicon_data = AmpliconCoverage.create(amplicon_name=row[3],
                                                    sample=sample,
                                                    library_name=samples[sample]['library_name'],
                                                    run_id=samples[sample]['run_id'],
                                                    program_name=program,
                                                    extraction=samples[sample]['extraction'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default=None)
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-p', '--program', help='Coverage estimation program', default='sambamba')
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

    root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)

    for sample in samples:
        sample_job = Job.wrapJobFn(process_sample_coverage, [args.address], "coveragestore", auth_provider,
                                   args.report, sample, args.program, samples, config,
                                   cores=1)
        root_job.addChild(sample_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)