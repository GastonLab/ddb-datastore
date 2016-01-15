from cassandra.cqlengine import columns
from cassandra.cqlengine.models import Model


class Variant(Model):
    chr = columns.Text(primary_key=True, partition_key=True)
    start = columns.Integer(primary_key=True, partition_key=True)
    end = columns.Integer(primary_key=True)
    ref = columns.Text(primary_key=True)
    alt = columns.Text(primary_key=True)
    sample = columns.Text(primary_key=True)

    # Simple Annotation Data
    callers = columns.List(columns.Text)
    type = columns.Text()
    subtype = columns.Text()

    somatic = columns.Boolean()
    germline = columns.Boolean()

    rs_id = columns.Text()
    rs_ids = columns.List(columns.Text)
    cosmic_ids = columns.List(columns.Text)

    gene = columns.Text()
    transcript = columns.Text()

    genes = columns.List(columns.Text)
    transcripts = columns.List(columns.Text)

    in_cosmic = columns.Boolean()
    in_clinvar = columns.Boolean()
    is_pathogenic = columns.Boolean()
    is_coding = columns.Boolean()
    is_lof = columns.Boolean()
    is_splicing = columns.Boolean()
    is_truncating = columns.Boolean()

    # Derived Annotation Data
    max_aaf = columns.Float()
    min_alt_depth = columns.Float()
    max_alt_depth = columns.Float()

    # Complex Annotation Data
    population_freqs = columns.Map(columns.Text, columns.Float)
    clinvar_data = columns.Map(columns.Text, columns.Text)

    # Variant Caller Data
    freebayes = columns.Map(columns.Text, columns.Text)
    mutect = columns.Map(columns.Text, columns.Text)
    scalpel = columns.Map(columns.Text, columns.Text)
    vardict = columns.Map(columns.Text, columns.Text)
