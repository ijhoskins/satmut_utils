#!/usr/bin/env/python
"""Runs read base quality masking of synthetic primer regions."""

import argparse
import logging
import os
import sys

import core_utils.file_utils as fu
from analysis.read_preprocessor import ReadMasker, DEDUP_FLAG, CDEDUP_FLAG
from core_utils.string_utils import none_or_str

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.3"
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
    parser.add_argument("-b", "--bam", type=str, help='BAM file containing alignments to mask.')

    parser.add_argument("-p", "--primers", type=str,
                        help='Primer feature file, e.g. BED, GFF, or GTF. Must have strand.')

    parser.add_argument("-d", "--deduplicate", action="store_true",
                        help='Flag indicating the reads have ArcherDX AMP library UMIs and should be deduplicated.')

    parser.add_argument("-cd", "--consensus_deduplicate", action="store_true",
                        help='Flag indicating consensus reads should be generated during deduplication. '
                             'Only has an effect when supplied with --deduplicate. For basic deduplication, omit.')

    parser.add_argument("-od", "--output_dir", type=str, default=".", help='Output directory for masked BAM')

    parser.add_argument("-o", "--output_prefix", type=none_or_str,
                        help='Optional output prefix for masked BAM. If not supplied uses prefix of input BAM.')

    parser.add_argument("-n", "--nthreads", type=int, default=ReadMasker.DEFAULT_NTHREADS,
                        help='Number threads for SAM/BAM file manipulation operations.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(bam, primers=None, dedup=DEDUP_FLAG, consensus_dedup=CDEDUP_FLAG, output_dir=".", output_prefix=None,
             nthreads=0):
    """Runs the read masker workflow.

    :param str bam: BAM file to mask
    :param str primers: feature file of primer locations for read masking and primer detection
    :param bool dedup: should the reads be deduplicated based on ArcherDX AMP library UMIs? Default False.
    :param bool consensus_dedup: were the reads consensus-deduplicated? Default False.
    :param str output_dir: optional output directory to store generated FASTQs and BAM
    :param str | None output_prefix: optional output prefix for the FASTQ(s) and BAM; if None, use same prefix as VCF
    :param int nthreads: number threads to use for SAM/BAM file manipulations. Default 0 (autodetect).
    :return str: output BAM filename
    """

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    out_prefix = output_prefix
    if output_prefix is None:
        out_prefix = os.path.basename(fu.remove_extension(bam))

    out_bam = fu.add_extension(os.path.join(output_dir, out_prefix), "masked.bam")

    # The read masker runs the workflow in __init__
    rm = ReadMasker(
        in_bam=bam, feature_file=primers, dedup=dedup, consensus_dedup=consensus_dedup,
        out_bam=out_bam, nthreads=nthreads)

    return rm.out_bam


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(bam=parsed_args["bam"], primers=parsed_args["primers"],
             dedup=parsed_args["deduplicate"], consensus_dedup=parsed_args["consensus_deduplicate"],
             output_dir=parsed_args["output_dir"], output_prefix=parsed_args["output_prefix"],
             nthreads=parsed_args["nthreads"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
