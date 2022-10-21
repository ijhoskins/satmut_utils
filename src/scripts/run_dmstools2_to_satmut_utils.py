#!/usr/bin/env python3
"""Preprocesses dms_tools2 dual barcoded subamplicon reads to satmut_utils compatible reads."""

import argparse
import logging
import os
import sys

from core_utils.file_utils import replace_extension
from satmut_utils.definitions import LOG_FORMATTER
from prototype.design_conversion import DmsTools2ToSatmutUtils

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

LOGFILE = replace_extension(os.path.basename(__file__), "log")
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(LOG_FORMATTER)
logger.addHandler(console_handler)


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-1", "--fastq1", type=str, required=True, help='R1 FASTQ.')

    parser.add_argument("-2", "--fastq2", type=str, required=True, help='R2 FASTQ.')

    parser.add_argument("-o", "--output_dir", type=str, default=DmsTools2ToSatmutUtils.DEFAULT_OUTDIR,
                        help='Optional output directory for converted FASTQs.')

    parser.add_argument("-l", "--umi_length", type=int, default=DmsTools2ToSatmutUtils.DEFAULT_UMI_LEN,
                        help='UMI/barcode length in dms_tools2 R1 and R2.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(r1_fastq, r2_fastq, umi_len=DmsTools2ToSatmutUtils.DEFAULT_UMI_LEN,
             outdir=DmsTools2ToSatmutUtils.DEFAULT_OUTDIR):
    """Runs the dms_tools2 to satmut_utils read conversion workflow.

    :param str r1_fastq: path of the R1 FASTQ
    :param str r2_fastq: path of the R2 FASTQ
    :param int umi_len: length of randomer UMI. Default 8.
    :param bool no_umis: do not add UMIs, only make reads flush with the pos_range. Default False
    :param str outdir: output directory. Default current directory.
    :param int nthreads: Number of threads to use for BAM operations. Default 0 (autodetect).
    """

    zipped_r1_fastq, zipped_r2_fastq = DmsTools2ToSatmutUtils(
        r1_fastq=r1_fastq, r2_fastq=r2_fastq, umi_len=umi_len, outdir=outdir).workflow()

    return zipped_r1_fastq, zipped_r2_fastq


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    outdir = parsed_args["output_dir"]
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    log_handler = logging.FileHandler(os.path.join(outdir, LOGFILE))
    log_handler.setFormatter(LOG_FORMATTER)
    logger.addHandler(log_handler)

    logger.info("Started %s" % sys.argv[0])

    workflow(r1_fastq=parsed_args["fastq1"], r2_fastq=parsed_args["fastq2"], umi_len=parsed_args["umi_length"],
             outdir=parsed_args["output_dir"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
