#!/usr/bin/env python3
"""Postprocesses satmut_utils VCF summary files to DiMSum nucleotide string-count matrix."""

import argparse
import logging
import os
import sys

from core_utils.file_utils import replace_extension
from satmut_utils.definitions import LOG_FORMATTER
from prototype.design_conversion import SatmutUtilsToDimSum

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
    parser.add_argument("-v", "--vcf_summary", type=str, required=True, help='satmut_utils *vcf.summary.txt file.')

    parser.add_argument("-r", "--reference", type=str, required=True, help='Reference FASTA.')

    parser.add_argument("-b", "--cds_bed", type=str, required=True, help='R2 FASTQ.')

    parser.add_argument("-o", "--output_dir", type=str, default=SatmutUtilsToDimSum.DEFAULT_OUTDIR,
                        help='Optional output directory for DiMSum count matrix.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(vcf_summary, ref, cds_bed, outdir=SatmutUtilsToDimSum.DEFAULT_OUTDIR):
    """Runs the satmut_utils to DiMSum postprocessing workflow.

    :param str vcf_summary: satmut_utils VCF summary.txt file to process
    :param str ref: reference FASTA
    :param str cds_bed: Single-line BED file with CDS annotation relative to the reference FASTA
    :param str outdir: output directory. Default current directory.
    :return str: filepath of the output table
    """

    out_mat = SatmutUtilsToDimSum(
        vcf_summary=vcf_summary, reference=ref, cds_bed=cds_bed, outdir=outdir).workflow()

    return out_mat


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

    workflow(vcf_summary=parsed_args["vcf_summary"], ref=parsed_args["reference"], cds_bed=parsed_args["cds_bed"],
             outdir=outdir)

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
