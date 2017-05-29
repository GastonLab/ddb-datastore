from cassandra.cqlengine import columns
from cassandra.cqlengine.models import Model


class Pipeline(Model):
    __keyspace__ = 'pipelinestore'
    status = columns.Text(primary_key=True, partition_key=True)
    type = columns.Text(primary_key=True, partition_key=True)
    name = columns.DateTime(primary_key=True, partition_key=True)
    date = columns.DateTime(primary_key=True)

    regions = columns.Text()
    snv_regions = columns.Text()
    indel_regions = columns.Text()
    reference = columns.Text()
    dict = columns.Text()
    indel1 = columns.Text()
    indel2 = columns.Text()
    dbsnp = columns.Text()
    cosmic = columns.Text()

    config = columns.Map(columns.Text, columns.Text)
