from cassandra.cqlengine import columns
from cassandra.cqlengine.models import Model


class Sample(Model):
    __keyspace__ = 'samplestore'
    sample = columns.Text(primary_key=True, partition_key=True)
    library_name = columns.Text(primary_key=True)
    run_id = columns.Text(primary_key=True)

    extraction = columns.Text()
    panel_name = columns.Text()
    target_pool = columns.Text()
    sequencer = columns.Text()

    pipeline_id = columns.Text()
    config = columns.Map(columns.Text, columns.Text)
