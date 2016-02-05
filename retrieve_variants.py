#!/usr/bin/env python

import sys
import argparse
from cassandra.cqlengine import connection

from variantstore import Variant


if __name__ == "__main__":
    connection.setup(['127.0.0.1'], "variantstore")
    for variant in Variant.objects.all():
        for key in variant.keys():
            if variant[key] is None:
                sys.stdout.write("{}: {} ({})\n".format(key, variant[key], type(variant[key])))
