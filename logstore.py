from cassandra.cqlengine import columns
from cassandra.cqlengine.models import Model


class Log(Model):
    __keyspace__ = 'logstore'
    sample = columns.Text(index=True, primary_key=True)
    library_name = columns.Text(index=True, primary_key=True)
    run_id = columns.Text(index=True, primary_key=True)
    analysis_id = columns.Text(index=True, primary_key=True)

    sequencer = columns.Text(index=True)
    target_pool = columns.Text(index=True)
    panel_name = columns.Text(index=True)
    initial_report_panel = columns.Text(index=True)
    extraction = columns.Text(index=True)
    date_annotated = columns.DateTime()

    log_data = columns.Map(columns.Text, columns.Blob)
