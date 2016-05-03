from cassandra.cqlengine import columns
from cassandra.cqlengine.models import Model


class Variant(Model):
    __keyspace__ = 'variantstore'
    reference_genome = columns.Text(primary_key=True, partition_key=True)
    chr = columns.Text(primary_key=True, partition_key=True)
    pos = columns.Integer(primary_key=True, partition_key=True)

    # Cluster Keys
    ref = columns.Text(primary_key=True)
    alt = columns.Text(primary_key=True)

    # Sample and Panel/Run Level data annotations
    sample = columns.Text(index=True, primary_key=True)
    library_name = columns.Text(index=True, primary_key=True)
    run_id = columns.Text(index=True, primary_key=True)
    date_annotated = columns.DateTime(index=True)
    sequencer = columns.Text(index=True)
    target_pool = columns.Text(index=True)
    panel_name = columns.Text(index=True)
    extraction = columns.Text(index=True)

    # Simple Annotation Data
    end = columns.Integer()
    callers = columns.List(columns.Text)
    type = columns.Text()
    subtype = columns.Text()

    somatic = columns.Boolean()
    germline = columns.Boolean()

    # Variant IDs
    rs_id = columns.Text()
    rs_ids = columns.List(columns.Text)
    cosmic_ids = columns.List(columns.Text)

    gene = columns.Text(index=True)
    transcript = columns.Text()
    exon = columns.Text()
    codon_change = columns.Text()
    aa_change = columns.Text()
    biotype = columns.Text(index=True)
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
    max_maf_all = columns.Float()
    max_maf_no_fin = columns.Float()
    min_depth = columns.Float()
    max_depth = columns.Float()
    min_som_aaf = columns.Float()
    max_som_aaf = columns.Float()
    variant_filters = columns.List(columns.Text)
    variant_reportable = columns.Text()
    panel_coverage = columns.List(columns.Text)
    amplicon_data = columns.Map(columns.Text, columns.Text)

    # Variant Caller Data
    freebayes = columns.Map(columns.Text, columns.Text)
    mutect = columns.Map(columns.Text, columns.Text)
    scalpel = columns.Map(columns.Text, columns.Text)
    vardict = columns.Map(columns.Text, columns.Text)
    scanindel = columns.Map(columns.Text, columns.Text)
    pindel = columns.Map(columns.Text, columns.Text)
    platypus = columns.Map(columns.Text, columns.Text)
    mutect2 = columns.Map(columns.Text, columns.Text)
    haplotypecaller = columns.Map(columns.Text, columns.Text)
    unifiedgenotype = columns.Map(columns.Text, columns.Text)
    itdseek = columns.Map(columns.Text, columns.Text)
    manta = columns.Map(columns.Text, columns.Text)


class SampleVariant(Model):
    __keyspace__ = 'variantstore'
    sample = columns.Text(primary_key=True, partition_key=True)
    run_id = columns.Text(primary_key=True, partition_key=True)
    library_name = columns.Text(primary_key=True)

    reference_genome = columns.Text(primary_key=True)
    chr = columns.Text(primary_key=True)
    pos = columns.Integer(primary_key=True)

    # Cluster Keys
    ref = columns.Text(primary_key=True)
    alt = columns.Text(primary_key=True)

    # Sample and Panel/Run Level data annotations
    date_annotated = columns.DateTime(primary_key=True, index=True)
    target_pool = columns.Text(index=True)
    panel_name = columns.Text(index=True)
    extraction = columns.Text(index=True)
    sequencer = columns.Text(index=True)

    # Simple Annotation Data
    end = columns.Integer()
    callers = columns.List(columns.Text)
    type = columns.Text()
    subtype = columns.Text()

    somatic = columns.Boolean()
    germline = columns.Boolean()

    # Variant IDs
    rs_id = columns.Text()
    rs_ids = columns.List(columns.Text)
    cosmic_ids = columns.List(columns.Text)

    gene = columns.Text(index=True)
    transcript = columns.Text()
    exon = columns.Text()
    codon_change = columns.Text()
    aa_change = columns.Text()
    biotype = columns.Text(index=True)
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
    amplicon_data = columns.Map(columns.Text, columns.Text)
    max_maf_all = columns.Float()
    max_maf_no_fin = columns.Float()
    min_depth = columns.Float()
    max_depth = columns.Float()
    min_som_aaf = columns.Float()
    max_som_aaf = columns.Float()
    variant_filters = columns.List(columns.Text)
    variant_categorization = columns.Text()

    # Variant Caller Data
    freebayes = columns.Map(columns.Text, columns.Text)
    mutect = columns.Map(columns.Text, columns.Text)
    scalpel = columns.Map(columns.Text, columns.Text)
    vardict = columns.Map(columns.Text, columns.Text)
    scanindel = columns.Map(columns.Text, columns.Text)
    pindel = columns.Map(columns.Text, columns.Text)
    platypus = columns.Map(columns.Text, columns.Text)
    mutect2 = columns.Map(columns.Text, columns.Text)
    haplotypecaller = columns.Map(columns.Text, columns.Text)
    unifiedgenotype = columns.Map(columns.Text, columns.Text)
    itdseek = columns.Map(columns.Text, columns.Text)
    manta = columns.Map(columns.Text, columns.Text)
