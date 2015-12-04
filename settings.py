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
    MONGO_DBNAME = 'apitest'
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
    'name': {
        'type': 'string',
        'minlength': 1,
        'maxlength': 50,
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
        'allowed': ["dna", "rna", "epigenetic", "other"],
        'required': True,
    },
    'status': {
        'type': 'string',
        'allowed': ["dev-template", "prod-template", "run"],
        'required': True,
    },
    'samples': {
        'type': 'dict',
        'schema': {
            'id': {'type': 'string'},
            'type': {'type': 'string'},
        },
    },
    'engine': {
        'type': 'string',
        'allowed': ["toil", "bcbio", "other"],
        'required': True,
    },
    'parameters': {
        'type': 'dict',
    },
    'workflow_file': {
        'type': 'media',
    },
    'date_added': {
        'type': 'datetime',
    },
    'date_modified': {
        'type': 'datetime',
    },
}

pipelines = {
    # 'title' tag used in item links. Defaults to the resource title minus
    # the final, plural 's' (works fine in most cases but not for 'people')
    'item_title': 'pipeline',

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
    'pipelines': pipelines,
}
