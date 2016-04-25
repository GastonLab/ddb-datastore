from cassandra.cqlengine import columns
from cassandra.cqlengine.models import Model


class SampleCoverage(Model):
    sample = columns.Text(primary_key=True, partition_key=True)

    library_name = columns.Text(primary_key=True)
    run_id = columns.Text(primary_key=True)


class AmpliconCoverage(Model):
    amplicon = columns.Text(primary_key=True, partition_key=True)


