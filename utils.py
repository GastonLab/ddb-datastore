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

    clinvar_data['significance'] = variant.INFO.get('clinvar_sig') or 'None'
    clinvar_data['pathogenic'] = variant.INFO.get('clinvar_pathogenic') or 'None'
    clinvar_data['hgvs'] = variant.INFO.get('clinvar_hgvs') or 'None'
    clinvar_data['revstatus'] = variant.INFO.get('clinvar_revstatus') or 'None'
    clinvar_data['origin'] = variant.INFO.get('clinvar_origin') or 'None'
    clinvar_data['disease'] = variant.INFO.get('clinvar_diseasename') or 'None'
    clinvar_data['accession'] = variant.INFO.get('clinvar_accession') or 'None'

    return clinvar_data


def get_cosmic_info(variant):
    cosmic_data = dict()

    cosmic_data['ids'] = variant.INFO.get('cosmic_ids') or 'None'
    cosmic_data['num_samples'] = variant.INFO.get('cosmic_numsamples') or 'None'
    cosmic_data['cds'] = variant.INFO.get('cosmic_cds') or 'None'
    cosmic_data['aa'] = variant.INFO.get('cosmic_aa') or 'None'
    cosmic_data['gene'] = variant.INFO.get('cosmic_gene') or 'None'

    return cosmic_data
