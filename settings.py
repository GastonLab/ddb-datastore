import os

if os.environ.get('PORT'):
    # Alternative configuration for a production server
    pass
else:
    MONGO_HOST = 'localhost'
    MONGO_PORT = 27017
    # MONGO_USERNAME = 'user'
    # MONGO_PASSWORD = 'user'
    MONGO_DBNAME = 'variantstore'
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
    'chr': {
        'type': 'string',
        'minlength': 1,
        'maxlength': 100,
        'required': True,
        'unique': True,
    },
    'start': {
        'type': 'int',
        'required': True,
    },
    'end': {
        'type': 'int',
        'required': True,
    },
    'ref': {
        'type': 'string',
        'minlength': 1,
        'maxlength': 100,
        'required': True,
        'unique': True,
    },
    'alt': {
        'type': 'string',
        'minlength': 1,
        'maxlength': 100,
        'required': True,
        'unique': True,
    },
    'desc': {
        'type': 'string',
        'minlength': 1,
        'maxlength': 100,
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
    'samples': {
        'type': 'dict',
        'schema': {
            'id': {'type': 'string'},
            'type': {'type': 'string'},
        },
    },
    'het_samples': {
        'type' 'list',
    },
    'hom_alt_samples': {
        'type' 'list',
    },
    'hom_ref_samples': {
        'type' 'list',
    },
    'max_af': {
        'type': 'float',
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
    'additional_lookup': {
        'url': 'regex("[\w]+")',
        'field': 'name'
    },

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
