from cassandra.cqlengine import columns
from cassandra.cqlengine.models import Model


class Variant(Model):
    chr = columns.Text(primary_key=True, partition_key=True)
    start = columns.Integer(primary_key=True, partition_key=True)
    end = columns.Integer(primary_key=True)
    ref = columns.Text(primary_key=True)
    alt = columns.Text(primary_key=True)
    sample = columns.Text(primary_key=True)
    library_name = columns.Text(primary_key=True)
    target_pool = columns.Text(primary_key=True)
    panel_name = columns.Text(primary_key=True)
    extraction = columns.Text(primary_key=True)
    reference_genome = columns.Text(primary_key=True)
    date_annotated = columns.DateTime(primary_key=True)

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
    exon = columns.Text()
    codon_change = columns.Text()
    aa_change = columns.Text()
    biotype = columns.Text()
    severity = columns.Text()
    impact = columns.Text()
    impact_so = columns.Text()

    genes = columns.List(columns.Text)
    transcripts_data = columns.Map(columns.Text, columns.Text)

    in_cosmic = columns.Boolean()
    in_clinvar = columns.Boolean()
    is_pathogenic = columns.Boolean()
    is_coding = columns.Boolean()
    is_lof = columns.Boolean()
    is_splicing = columns.Boolean()

    # Complex Annotation Data
    population_freqs = columns.Map(columns.Text, columns.Float)
    clinvar_data = columns.Map(columns.Text, columns.Text)
    cosmic_data = columns.Map(columns.Text, columns.Text)
    max_aaf_all = columns.Float()
    max_aaf_no_fin = columns.Float()

    # Variant Caller Data
    freebayes = columns.Map(columns.Text, columns.Text)
    mutect = columns.Map(columns.Text, columns.Text)
    scalpel = columns.Map(columns.Text, columns.Text)
    vardict = columns.Map(columns.Text, columns.Text)
    scanindel = columns.Map(columns.Text, columns.Text)
    platypus = columns.Map(columns.Text, columns.Text)
    mutect2 = columns.Map(columns.Text, columns.Text)
    haplotypecaller = columns.Map(columns.Text, columns.Text)
    unifiedgenotype = columns.Map(columns.Text, columns.Text)
    itdseek = columns.Map(columns.Text, columns.Text)
    manta = columns.Map(columns.Text, columns.Text)
