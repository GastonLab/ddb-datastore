from cassandra.cqlengine import columns
from cassandra.cqlengine.models import Model


class Pipeline(Model):
    __keyspace__ = 'pipelinestore'
    run_type = columns.Text(primary_key=True, partition_key=True)
    status = columns.Text(primary_key=True, partition_key=True)
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

    general_config = columns.Map(columns.Text, columns.Text)
    sambamba_config = columns.Map(columns.Text, columns.Text)
    vcfanno_config = columns.Map(columns.Text, columns.Text)
    bwa_config = columns.Map(columns.Text, columns.Text)
    ensemble_config = columns.Map(columns.Text, columns.Text)
    gatk_config = columns.Map(columns.Text, columns.Text)
    scanindel_config = columns.Map(columns.Text, columns.Text)
    pindel_config = columns.Map(columns.Text, columns.Text)
    picard_config = columns.Map(columns.Text, columns.Text)
    freebayes_config = columns.Map(columns.Text, columns.Text)
    mutect_config = columns.Map(columns.Text, columns.Text)
    vardict_config = columns.Map(columns.Text, columns.Text)
    scalpel_config = columns.Map(columns.Text, columns.Text)
    indelminer_config = columns.Map(columns.Text, columns.Text)
    platypus_config = columns.Map(columns.Text, columns.Text)
    gemini_config = columns.Map(columns.Text, columns.Text)
    snpeff_config = columns.Map(columns.Text, columns.Text)
    bcftools_config = columns.Map(columns.Text, columns.Text)
    vcftools_config = columns.Map(columns.Text, columns.Text)
    samtools_config = columns.Map(columns.Text, columns.Text)
    fastqc_config = columns.Map(columns.Text, columns.Text)
    vt_config = columns.Map(columns.Text, columns.Text)
