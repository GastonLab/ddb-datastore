#!/usr/bin/env python

import argparse
from cassandra.cqlengine.management import sync_table
from cassandra.cqlengine import connection
from variantstore import Variant
from variantstore import SampleVariant

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    args = parser.parse_args()

    connection.setup([args.address], "variantstore")

    sync_table(Variant)
    sync_table(SampleVariant)
