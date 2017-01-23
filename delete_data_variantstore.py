#!/usr/bin/env python

import sys
import argparse
import getpass

from ddb import configuration
from variantstore import SampleVariant

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

        variants = SampleVariant.objects.timeout(None).filter(
            SampleVariant.reference_genome == config['genome_version'],
            SampleVariant.sample == samples[sample]['sample_name'],
            SampleVariant.run_id == samples[sample]['run_id'],
            SampleVariant.library_name == samples[sample]['library_name'],
        ).allow_filtering()

        for variant in variants:
            variant.delete()
