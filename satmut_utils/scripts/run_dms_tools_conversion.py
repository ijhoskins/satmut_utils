#!/usr/bin/env/python
"""Preprocesses amplicon/tile alignments to meet Enrich2 and dms_tools2 input requirements."""

import argparse
import logging
import sys

from prototype.design_conversion import ConvertToDmstools

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

__logger = logging.getLogger(__name__)

# Consider putting the following in a logging config file
__logger.setLevel(logging.DEBUG)
__fhandler = logging.FileHandler("stderr.log")
__fhandler.setLevel(logging.DEBUG)
__formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
__fhandler.setFormatter(__formatter)
__logger.addHandler(__fhandler)


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

    parser.add_argument("-o", "--output_dir", type=str, default=ConvertToDmstools.DEFAULT_OUTDIR,
                        help='Output directory for FASTQs and BAM.')

    parser.add_argument("-l", "--umi_length", type=int, default=ConvertToDmstools.DEFAULT_UMI_LEN,
                        help='UMI length, for addition to 5\' end of each read.')

    parser.add_argument("-n", "--no_umis", action="store_true",
                        help='Do not append UMIs to the 5\' end of reads.')

    parser.add_argument("-j", "--nthreads", type=int, default=ConvertToDmstools.DEFAULT_NTHREADS,
                        help='Number of threads to use for BAM sorting operations. Default %i, autodetect.'
                             % ConvertToDmstools.DEFAULT_NTHREADS)

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(in_bam, ref, pos_range, umi_len=ConvertToDmstools.DEFAULT_UMI_LEN, no_umis=ConvertToDmstools.DEFAULT_NO_UMI,
             outdir=ConvertToDmstools.DEFAULT_OUTDIR, nthreads=ConvertToDmstools.DEFAULT_NTHREADS):
    """Runs the dms_tools conversion workflow.

    :param str in_bam: input alignments for a single tile
    :param str ref: reference FASTA
    :param str pos_range: 1-based positions flush with codons spanning the target
    :param int umi_len: length of randomer UMI. Default 8.
    :param bool no_umis: do not add UMIs, only make reads flush with the pos_range. Default False
    :param str outdir: output directory. Default current directory.
    :param int nthreads: Number of threads to use for BAM operations. Default 0 (autodetect).
    """

    zipped_r1_fastq, zipped_r2_fastq = ConvertToDmstools(
        in_bam=in_bam, ref=ref, pos_range=pos_range, umi_len=umi_len, no_umis=no_umis,
        outdir=outdir, nthreads=nthreads).workflow()

    return zipped_r1_fastq, zipped_r2_fastq


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(in_bam=parsed_args["in_bam"], ref=parsed_args["reference"], pos_range=parsed_args["pos_range"],
             umi_len=parsed_args["umi_length"], no_umis=parsed_args["no_umis"],
             outdir=parsed_args["output_dir"], nthreads=parsed_args["nthreads"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
