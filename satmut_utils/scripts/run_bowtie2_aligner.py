#!/usr/bin/env/python
"""Runs bowtie2 alignment."""

import argparse
import logging
import os
import sys

from ..analysis.aligners import BowtieConfig, Bowtie2, DEFAULT_TEMPDIR
from ..core_utils.string_utils import none_or_str


__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

__logger = logging.getLogger(__name__)

# Consider putting the following in a logging config file
__logger.setLevel(logging.DEBUG)
__fhandler = logging.FileHandler("stderr.log")
__fhandler.setLevel(logging.DEBUG)
_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
__fhandler.setFormatter(_formatter)
__logger.addHandler(__fhandler)


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-f1", "--fast1", type=str, required=True, help='R1 FASTA or FASTQ.')

    parser.add_argument("-f2", "--fast2", type=str, required=False, default=None, help='R2 FASTA or FASTQ.')

    parser.add_argument("-r", "--ref", type=str, required=True,
                        help='Reference FASTA. Should contain a bowtie2 FM index file set.')

    parser.add_argument("-d", "--outdir", type=none_or_str, required=False, default=".",
                        help='Optional output directory.')

    parser.add_argument("-o", "--outbam", type=none_or_str, required=False, default=None,
                        help='Optional output BAM filename.')

    parser.add_argument("-l", "--local", action="store_true",
                        help='Run a local alignment as opposed to a global alignment.')

    parser.add_argument("-t", "--nthreads", type=int, required=False, default=BowtieConfig.DEFAULT_NTHREADS,
                        help='Number of threads to use for alignment.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(f1, ref, f2=None, outdir=".", outbam=None, local=True, nthreads=BowtieConfig.DEFAULT_NTHREADS):
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

    # Need to change dir to scratch in case temp BAMs are created, as bowtie2 does not have an option for setting
    # the temp directory
    call_dir = os.getcwd()
    os.chdir(DEFAULT_TEMPDIR)

    bc = BowtieConfig(ref, local, nthreads)
    bt = Bowtie2(config=bc, f1=f1, f2=f2, output_dir=outdir, output_bam=outbam)

    os.chdir(call_dir)

    return bt


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(f1=parsed_args["fast1"], ref=parsed_args["ref"], f2=parsed_args["fast2"], outdir=parsed_args["outdir"],
             outbam=parsed_args["outbam"], local=parsed_args["local"], nthreads=parsed_args["nthreads"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
