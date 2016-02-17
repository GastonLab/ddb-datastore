import geneimpacts


def get_effects(variant, annotation_keys):
    effects = []
    effects += [geneimpacts.SnpEff(e, annotation_keys) for e in variant.INFO.get("ANN").split(",")]

    return effects


def get_top_impact(effects):
    top_impact = geneimpacts.Effect.top_severity(effects)

    if isinstance(top_impact, list):
        top_impact = top_impact[0]

    return top_impact


def get_genes(effects):
    genes_list = []

    for effect in effects:
        if effect.gene not in genes_list:
            genes_list.append(effect.gene)

    return genes_list


def get_transcript_effects(effects):
    transcript_effects = dict()

    for effect in effects:
        transcript_effects[effect.transcript] = "{biotype}|{effect}".format(biotype=effect.biotype,
                                                                            effect=effect.impact_severity)

    return transcript_effects


def get_clinvar_info(variant):
    clinvar_data = dict()

    clinvar_data['significance'] = variant.INFO.get('clinvar_sig')
    clinvar_data['pathogenic'] = variant.INFO.get('clinvar_pathogenic')
    clinvar_data['hgvs'] = variant.INFO.get('clinvar_hgvs')
    clinvar_data['revstatus'] = variant.INFO.get('clinvar_revstatus')
    clinvar_data['origin'] = variant.INFO.get('clinvar_origin')
    clinvar_data['disease'] = variant.INFO.get('clinvar_diseasename')
    clinvar_data['accession'] = variant.INFO.get('clinvar_accession')

    return clinvar_data


def get_cosmic_info(variant):
    cosmic_data = dict()

    cosmic_data['ids'] = variant.INFO.get('cosmic_ids')
    cosmic_data['num_samples'] = variant.INFO.get('cosmic_numsamples')
    cosmic_data['cds'] = variant.INFO.get('cosmic_cds')
    cosmic_data['aa'] = variant.INFO.get('cosmic_aa')
    cosmic_data['gene'] = variant.INFO.get('cosmic_gene')

    return cosmic_data
