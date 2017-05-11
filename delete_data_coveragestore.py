#!/usr/bin/env python

import sys
import argparse
import getpass
import utils

from ddb import configuration
from coveragestore import AmpliconCoverage
from coveragestore import SampleCoverage

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
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
        sys.stdout.write("Deleting variants for sample {}\n".format(sample))
        report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}" \
                            "".format(samples[sample]['panel'], samples[sample]['report'])
        target_amplicons = utils.get_target_amplicons(report_panel_path)

        for amplicon in target_amplicons:
            coverage_data = SampleCoverage.objects.timeout(None).filter(
                SampleCoverage.sample == samples[sample]['sample_name'],
                SampleCoverage.amplicon == amplicon,
                SampleCoverage.run_id == samples[sample]['run_id'],
                SampleCoverage.library_name == samples[sample]['library_name'],
                SampleCoverage.program_name == "sambamba"
            ).allow_filtering()

            for coverage in coverage_data:
                coverage.delete()
