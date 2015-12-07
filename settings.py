__author__ = 'dgaston'

import os

if os.environ.get('PORT'):
    # Alternative configuration for a production server
    pass
else:
    MONGO_HOST = 'localhost'
    MONGO_PORT = 27017
    # MONGO_USERNAME = 'user'
    # MONGO_PASSWORD = 'user'
    MONGO_DBNAME = 'variantservice'
    DEBUG = True
    VERSIONING = True
    SOFT_DELETE = True

# Enable reads (GET), inserts (POST) and DELETE for resources/collections
# (if you omit this line, the API will default to ['GET'] and provide
# read-only access to the endpoint).
RESOURCE_METHODS = ['GET', 'POST', 'DELETE']

# Enable reads (GET), edits (PATCH), replacements (PUT) and deletes of
# individual items  (defaults to read-only item access).
ITEM_METHODS = ['GET', 'PATCH', 'PUT', 'DELETE']

schema = {
    # Schema definition, based on Cerberus grammar. Check the Cerberus project
    # (https://github.com/nicolaiarocci/cerberus) for details.
    # 'variant_id': {
    #     'type': 'string',
    #     'minlength': 1,
    #     'maxlength': 100,
    #     'required': True,
    #     'unique': True,
    # },
    'chr': {
        'type': 'string',
        'minlength': 1,
        'maxlength': 100,
        'required': True,
    },
    'start': {
        'type': 'int',
        'required': True,
    },
    'end': {
        'type': 'int',
        'required': True,
    },
    'type': {
        'type': 'string',
        'allowed': ["snv", "indel", "sv"],
        'required': True,
    },
    'subtype': {
        'type': 'string',
        'allowed': ["tr", "tv", "ins", "del"],
        'required': True,
    },
    'ref': {
        'type': 'string',
        'minlength': 1,
        'maxlength': 1000,
        'required': True,
    },
    'alt': {
        'type': 'string',
        'minlength': 1,
        'maxlength': 1000,
        'required': True,
    },
    'max_aaf': {
        'type': 'float',
    },
    'het_samples': {
        'type': 'list',
    },
    'hom_alt_samples': {
        'type': 'list',
    },
    'hom_ref_samples': {
        'type': 'list',
    },
    'sample_quals': {
        'type': 'dict',
        'schema': {
            'sample_id': {'type': 'string'},
            'variant_caller': {'type': 'string'},
            'qual': {'type': 'float'}
        },
    },
    'sample_filters': {
        'type': 'dict',
        'schema': {
            'sample_id': {'type': 'string'},
            'variant_caller': {'type': 'string'},
            'filter': {'type': 'string'}
        },
    },
    'sample_genotypes': {
        'type': 'dict',
        'schema': {
            'sample_id': {'type': 'string'},
            'variant_caller': {'type': 'string'},
            'genotype': {'type': 'float'},
            'genotype_qual': {'type': 'int'},
            'alt_depth': {'type': 'int'},
            'ref_depth': {'type': 'int'},
            'depth': {'type': 'int'}
        },
    },
    'cosmic_ids': {
        'type': 'list',
    },
    'gene': {
        'type': 'string',
    },
    'transcript': {
        'type': 'string',
    },
    'pfam_domain': {
        'type': 'string',
    },
    'exon': {
        'type': 'string',
    },
    'biotype': {
        'type': 'string',
    },
    'codon_change': {
        'type': 'string',
    },
    'aa_change': {
        'type': 'string',
    },
    'aa_length': {
        'type': 'int',
    },
    'func_effects': {
        'type': 'dict',
        'schema': {
            'impact': {'type': 'string'},
            'impact_so': {'type': 'string'},
            'impact_severity': {'type': 'string'},
            'polyphen_pred': {'type': 'string'},
            'polyphen_score': {'type': 'float'},
            'sift_pred': {'type': 'string'},
            'sift_score': {'type': 'float'},
            'gerp_bp_score': {'type': 'float'},
            'cadd_raw': {'type': 'float'},
            'cadd_scaled': {'type': 'float'},
            'fitcons': {'type': 'float'}
        },
    },
    'frequencies': {
        'type': 'dict',
        'schema': {
            'aaf_esp_ea': {'type': 'float'},
            'aaf_esp_aa': {'type': 'float'},
            'aaf_esp_all': {'type': 'float'},
            'aaf_1kg_amr': {'type': 'float'},
            'aaf_1kg_eas': {'type': 'float'},
            'aaf_1kg_sas': {'type': 'float'},
            'aaf_1kg_afr': {'type': 'float'},
            'aaf_1kg_eur': {'type': 'float'},
            'aaf_1kg_all': {'type': 'float'},
            'aaf_adj_exac_all': {'type': 'float'},
            'aaf_adj_exac_afr': {'type': 'float'},
            'aaf_adj_exac_amr': {'type': 'float'},
            'aaf_adj_exac_eas': {'type': 'float'},
            'aaf_adj_exac_fin': {'type': 'float'},
            'aaf_adj_exac_nfe': {'type': 'float'},
            'aaf_adj_exac_oth': {'type': 'float'},
            'aaf_adj_exac_sas': {'type': 'float'},
            'exac_num_het': {'type': 'int'},
            'exac_num_hom_alt': {'type': 'int'},
            'exac_num_chroms': {'type': 'int'}
        },
    },
    'clin_info': {
        'type': 'dict',
        'schema': {
            'clinvar_disease_name': {'type': 'string'},
            'clinvar_dbsource': {'type': 'string'},
            'clinvar_dbsource_id': {'type': 'string'},
            'clinvar_origin': {'type': 'string'},
            'clinvar_dsdb': {'type': 'string'},
            'clinvar_dsdbid': {'type': 'string'},
            'clinvar_in_locus_spec_db': {'type': 'string'},
            'clinvar_on_diag_assay': {'type': 'string'}
        },
    },
    'sv_info': {
        'type': 'dict',
        'schema': {
            'sv_cipos_start_left': {'type': 'int'},
            'sv_cipos_end_left': {'type': 'int'},
            'sv_cipos_start_right': {'type': 'int'},
            'sv_cipos_end_right': {'type': 'int'},
            'sv_length': {'type': 'int'},
            'sv_is_precise': {'type': 'boolean'},
            'sv_tool': {'type': 'string'},
            'sv_evidence_type': {'type': 'string'},
            'sv_event_id': {'type': 'string'},
            'sv_mate_id': {'type': 'string'},
            'sv_strand': {'type': 'string'}
        },
    },
    'error_assessment': {
        'type': 'dict',
        'schema': {
            'grc': {'type': 'string'},
            'gms_illumina': {'type': 'float'},
            'gms_solid': {'type': 'float'},
            'gms_iontorrent': {'type': 'float'},
            'in_cse': {'type': 'boolean'}
        },
    },
    'cigar': {
        'type': 'string',
    },
    'rmsk': {
        'type': 'string',
    },
    'qual_depth': {
        'type': 'float',
    },
    'haplotype_score': {
        'type': 'float',
    },
    'strand_bias': {
        'type': 'float',
    },
    'rms_map_qual': {
        'type': 'float',
    },
    'rms_bq': {
        'type': 'float',
    },
    'recomb_rate': {
        'type': 'float',
    },
    'frac_reads_w_dels': {
        'type': 'float',
    },
    'num_mapq_zero': {
        'type': 'int',
    },
    'is_conserved': {
        'type': 'boolean',
    },
    'is_coding': {
        'type': 'boolean',
    },
    'is_exonic': {
        'type': 'boolean',
    },
    'is_lof': {
        'type': 'boolean',
    },
    'is_splicing': {
        'type': 'boolean',
    },
    'is_somatic': {
        'type': 'boolean',
    },
    'is_missense_nonsense': {
        'type': 'boolean',
    },
    'in_hom_run': {
        'type': 'boolean',
    },
    'in_dbsnp': {
        'type': 'boolean',
    },
    'in_esp': {
        'type': 'boolean',
    },
    'in_1kg': {
        'type': 'boolean',
    },
    'in_exac': {
        'type': 'boolean',
    },
    'in_omim': {
        'type': 'boolean',
    },
    'in_cpg_island': {
        'type': 'boolean',
    },
    'in_segdup': {
        'type': 'boolean',
    },
    'in_cse': {
        'type': 'boolean',
    },
    'info_string': {
        'type': 'dict',
    },
}

variants = {
    # 'title' tag used in item links. Defaults to the resource title minus
    # the final, plural 's' (works fine in most cases but not for 'people')
    'item_title': 'variant',

    # by default the standard item entry point is defined as
    # '/people/<ObjectId>'. We leave it untouched, and we also enable an
    # additional read-only entry point. This way consumers can also perform
    # GET requests at '/people/<lastname>'.
    # 'additional_lookup': {
    #     'url': 'regex("[\w]+")',
    #     'field': 'variant_id'
    # },

    # We choose to override global cache-control directives for this resource.
    'cache_control': 'max-age=10,must-revalidate',
    'cache_expires': 10,

    # most global settings can be overridden at resource level
    'resource_methods': ['GET', 'POST'],

    'schema': schema
}

DOMAIN = {
    'variants': variants,
}
