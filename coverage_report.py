#!/usr/bin/env python

import sys
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
from coveragestore import SampleCoverage
from coveragestore import AmpliconCoverage
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', help="File containing list of amplicon names to check")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default='report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    args.logLevel = "INFO"

    amplicons_list = list()
    target_amplicons = utils.get_target_amplicons(args.list)
    for amplicon in target_amplicons:
        if amplicon not in amplicons_list:
            amplicons_list.append(amplicon)

    for amplicon in target_amplicons:
        coverage_data = SampleCoverage.objects.timeout(None).filter(
            SampleCoverage.sample == samples[sample][library]['sample_name'],
            SampleCoverage.amplicon == amplicon,
            SampleCoverage.run_id == samples[sample][library]['run_id'],
            SampleCoverage.library_name == samples[sample][library]['library_name'],
            SampleCoverage.program_name == "sambamba"
        )
        ordered_amplicons = coverage_data.order_by('amplicon', 'run_id').limit(coverage_data.count() + 1000)
        for result in ordered_amplicons:
            reportable_amplicons.append(result)
            target_amplicon_coverage[amplicon] = result
            ordered_amplicon_coverage.append(result)
