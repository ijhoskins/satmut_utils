#!/usr/bin/env python3
"""Runs bowtie2 alignment."""

import argparse
import logging
import os
import sys

from analysis.aligners import BowtieConfig, Bowtie2, DEFAULT_TEMPDIR
from core_utils.file_utils import replace_extension
from core_utils.string_utils import none_or_str
from satmut_utils.definitions import DEFAULT_QUALITY_OFFSET, LOG_FORMATTER

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
    parser.add_argument("-1", "--fast1", type=str, required=True, help='R1 FASTA or FASTQ.')

    parser.add_argument("-2", "--fast2", type=str, required=False, default=None, help='R2 FASTA or FASTQ.')

    parser.add_argument("-r", "--ref", type=str, required=True,
                        help='Reference FASTA. Should contain a bowtie2 FM index file set.')

    parser.add_argument("-d", "--output_dir", type=str, required=False, default=Bowtie2.DEFAULT_OUTDIR,
                        help='Optional output directory. Default current working directory.')

    parser.add_argument("-o", "--outbam", type=none_or_str, required=False, default=Bowtie2.DEFAULT_OUTBAM,
                        help='Optional output BAM filename written to --output_dir. Default use basename of FASTQs.')

    parser.add_argument("-l", "--local", action="store_true",
                        help='Run a local alignment as opposed to a global alignment.')

    parser.add_argument("-j", "--nthreads", type=int, required=False, default=BowtieConfig.DEFAULT_NTHREADS,
                        help='Number of threads to use for alignment.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(f1, ref, f2=None, outdir=Bowtie2.DEFAULT_OUTDIR, outbam=Bowtie2.DEFAULT_OUTBAM,
             local=BowtieConfig.DEFAULT_LOCAL, nthreads=BowtieConfig.DEFAULT_NTHREADS):
    """Runs the Bowtie2 alignment workflow.

    :param str f1: path to FASTA or FASTQ 1
    :param str ref: path of indexed reference FASTA
    :param str | None f2: optional path to FASTA or FASTQ 2
    :param str outdir: optional output directory to write the output BAM to. Default current directory.
    :param str | None outbam: full file path of an output BAM to create; if None, use basename of f1, f2
    :param bool local: should a local alignment be done instead of global alignment (default True)
    :param int nthreads: number of threads to use for alignment
    :return analysis.aligners.Bowtie2: aligner object
    """

    # Need to change to tempdir in case temp BAMs are created, as bowtie2 does not have an option for setting
    # the temp directory; make sure to handle relative paths
    f1_full = os.path.abspath(f1)
    f2_full = os.path.abspath(f2)
    outdir_full = os.path.abspath(outdir)

    call_dir = os.getcwd()
    os.chdir(DEFAULT_TEMPDIR)

    bc = BowtieConfig(ref, local, nthreads, DEFAULT_QUALITY_OFFSET)
    bt = Bowtie2(config=bc, f1=f1_full, f2=f2_full, output_dir=outdir_full, output_bam=outbam)

    os.chdir(call_dir)

    return bt


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

    workflow(f1=parsed_args["fast1"], ref=parsed_args["ref"], f2=parsed_args["fast2"], outdir=parsed_args["output_dir"],
             outbam=parsed_args["outbam"], local=parsed_args["local"], nthreads=parsed_args["nthreads"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()

