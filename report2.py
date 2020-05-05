#!/usr/bin/env python
import getpass
import xlsxwriter

from toil.common import Toil
from toil.job import Job

def sampleReport(sample, memory="4G", cores=4, disk="1G"):
    return "Success"


if __name__=="__main__":
    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument('-s', '--samples_file',
                        help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration',
                        help="Configuration file for various settings")
    parser.add_argument('-r', '--report',
                        help="Root name for reports (per sample)",
                        default='report')
    parser.add_argument('-a', '--address',
                        help="IP Address for Cassandra connection",
                        default='127.0.0.1')
    parser.add_argument('-u', '--username',
                        help='Cassandra username for login',
                        default=None)
    parser.add_argument('-d', '--min_depth',
                        help='Minimum depth threshold for variant reporting',
                        default=200.0)
    parser.add_argument('-g', '--good_depth',
                        help='Floor for good depth of coverage',
                        default=500.0)
    parser.add_argument('-t', '--min_somatic_var_freq',
                        help='Minimum reportable somatic variant frequency',
                        default=0.01)
    parser.add_argument('-p', '--max_pop_freq',
                        help='Maximum allowed population allele frequency',
                        default=0.005)

    options = parser.parse_args()
    options.logLevel = "INFO"
    options.clean = "always"

    root_job = Job.wrapFn()
    for sample in samples:
        sample_job = Job.wrapFn(sampleReport, sample)

    with Toil(options) as workflow:
        if not toil.options.restart:
            print("Starting Toil Reporting Workflow")
            toil.start(root_job)
        else:
            toil.restart()
