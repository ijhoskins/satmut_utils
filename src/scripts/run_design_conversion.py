#!/usr/bin/env python3
"""Preprocesses amplicon/tile alignments to meet dms_tools2 and Enrich2 input requirements."""

import argparse
import logging
import os
import sys

from core_utils.file_utils import replace_extension
from satmut_utils.definitions import LOG_FORMATTER
from prototype.design_conversion import DesignConverter

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
    parser.add_argument("-b", "--in_bam", type=str, required=True,
                        help='BAM file of reads in a single tile to convert. Intersect any alignments with a single '
                             'target tile using bedtools -f (if alignments consist of multiple PCR tiles).')

    parser.add_argument("-r", "--reference", type=str, required=True, help='Corresponding reference FASTA.')

    parser.add_argument("-p", "--pos_range", type=str, required=True,
                        help='Comma-delimited 1-based start and end positions of the amplicon/tile. '
                             'Reads will be made flush with these coordinates.')

    parser.add_argument("-o", "--output_dir", type=str, default=DesignConverter.DEFAULT_OUTDIR,
                        help='Optional output directory for FASTQs and globally re-aligned BAM. '
                             'Default current working directory.')

    parser.add_argument("-l", "--umi_length", type=int, default=DesignConverter.DEFAULT_UMI_LEN,
                        help='UMI/barcode length, for addition to 5\' end of each read in dms_tools2 conversion.')

    parser.add_argument("-n", "--no_umis", action="store_true",
                        help='Do not append UMIs to the 5\' end of reads. Use for conversion to Enrich2 format.')

    parser.add_argument("-j", "--nthreads", type=int, default=DesignConverter.DEFAULT_NTHREADS,
                        help='Number of threads to use for BAM sorting operations. Default %i, autodetect.'
                             % DesignConverter.DEFAULT_NTHREADS)

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(in_bam, ref, pos_range, umi_len=DesignConverter.DEFAULT_UMI_LEN, no_umis=DesignConverter.DEFAULT_NO_UMI,
             outdir=DesignConverter.DEFAULT_OUTDIR, nthreads=DesignConverter.DEFAULT_NTHREADS):
    """Runs the dms_tools2/Enrich2 conversion workflow.

    :param str in_bam: input alignments for a single tile
    :param str ref: reference FASTA
    :param str pos_range: 1-based positions flush with codons spanning the target
    :param int umi_len: length of randomer UMI. Default 8.
    :param bool no_umis: do not add UMIs, only make reads flush with the pos_range. Default False
    :param str outdir: output directory. Default current directory.
    :param int nthreads: Number of threads to use for BAM operations. Default 0 (autodetect).
    """

    zipped_r1_fastq, zipped_r2_fastq = DesignConverter(
        in_bam=in_bam, ref=ref, pos_range=pos_range, umi_len=umi_len, no_umis=no_umis,
        outdir=outdir, nthreads=nthreads).workflow()

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

    workflow(in_bam=parsed_args["in_bam"], ref=parsed_args["reference"], pos_range=parsed_args["pos_range"],
             umi_len=parsed_args["umi_length"], no_umis=parsed_args["no_umis"],
             outdir=parsed_args["output_dir"], nthreads=parsed_args["nthreads"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
