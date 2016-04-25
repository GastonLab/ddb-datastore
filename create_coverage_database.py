#!/usr/bin/env python

import argparse
import getpass
from cassandra.cqlengine.management import sync_table
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider
from coveragestore import SampleCoverage
from coveragestore import AmpliconCoverage

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    args = parser.parse_args()

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "coveragestore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "coveragestore")

    sync_table(SampleCoverage)
    sync_table(AmpliconCoverage)
