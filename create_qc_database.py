#!/usr/bin/env python

import argparse
import getpass
from cassandra import query
from cassandra.cqlengine.management import sync_table
from cassandra.cqlengine.management import create_keyspace_simple
from cassandra.cluster import Cluster
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider
from qcstore import QC

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
    create_keyspace_simple("qcstore", args.replication_factor)

    sync_table(QC)
