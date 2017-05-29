#!/usr/bin/env python

import argparse
import getpass

from cassandra import query
from cassandra.auth import PlainTextAuthProvider
from cassandra.cluster import Cluster
from cassandra.cqlengine import connection
from cassandra.cqlengine.management import create_keyspace_simple
from cassandra.cqlengine.management import sync_table
from variantstore import SampleVariant
from variantstore import TargetVariant
from variantstore import Variant

from ddb_data.coveragestore import AmpliconCoverage
from ddb_data.coveragestore import SampleCoverage

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-r', '--replication_factor', help="Cassandra replication factor", default=3)
    args = parser.parse_args()

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        cluster = Cluster([args.address], auth_provider=auth_provider)
        session = cluster.connect()
        session.row_factory = query.dict_factory
    else:
        cluster = Cluster([args.address])
        session = cluster.connect()
        session.row_factory = query.dict_factory

    connection.set_session(session)
    create_keyspace_simple("variantstore_dev", args.replication_factor)
    create_keyspace_simple("coveragestore_dev", args.replication_factor)

    sync_table(Variant)
    sync_table(SampleVariant)
    sync_table(TargetVariant)
    sync_table(SampleCoverage)
    sync_table(AmpliconCoverage)
