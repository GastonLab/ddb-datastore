#!/usr/bin/env python

import sys
import dill
import logging
import argparse
import utils
import cyvcf2
import multiprocess

from multiprocess import Pool
from cyvcf2 import VCF


def worker(variant):
    # print "***"
    return variant['chrom'], variant['ref'], variant['alt']


if __name__ == "__main__":
    arguments = list()
    sys.stdout.write("Creating VCF object\n")
    vcf = VCF("57295.vcfanno.snpEff.GRCh37.75.vcf")
    sys.stdout.write("Iterating over VCF records\n")
    i = 0
    for record in vcf:
        record_dict = {'chrom': record.CHROM, 'ref': record.REF, 'alt': record.ALT}
        arguments.append(record_dict)
        i += 1

    pool = Pool(processes=6)
    sys.stdout.write("Running map_async for {} results\n".format(i))
    result = pool.map_async(worker, arguments)
    sys.stdout.write("Getting results\n")
    report_variants = result.get()
    sys.stdout.write("Closing\n")
    pool.close()
    pool.join()

    sys.stdout.write("Iterating over results\n")
    i = 0
    for result in report_variants:
        i += 1
    sys.stdout.write("Processed {} results\n".format(i))
