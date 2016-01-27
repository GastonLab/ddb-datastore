from cassandra.cqlengine.management import sync_table
from cassandra.cqlengine import connection
from variantstore import Variant

connection.setup(['127.0.0.1'], "variantstore")
sync_table(Variant)
