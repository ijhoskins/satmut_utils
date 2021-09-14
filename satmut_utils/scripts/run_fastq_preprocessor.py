#!/usr/bin/env/python
"""Runs FASTQ preprocessing (adapter, quality trimming)."""

import argparse
import logging
import sys

from analysis.read_preprocessor import FastqPreprocessor

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
    parser.add_argument("-f1", "--fastq1", type=str, required=True, help='R1 FASTQ')

    parser.add_argument("-f2", "--fastq2", type=str, required=True, help='R2 FASTQ')

    parser.add_argument("-a5", "--r1_fiveprime_adapters", required=True, type=str, help='Comma-delimited R1 5\' adapters.')

    parser.add_argument("-a3", "--r1_threeprime_adapters", required=True, type=str, help='Comma-delimited R1 3\' adapters.')

    parser.add_argument("-o", "--outdir", type=str, help='Output directory.')

    parser.add_argument("-c", "--cores", type=int, default=FastqPreprocessor.NCORES,
                        help='Number of cores to use for cutadapt. Default %i.' % FastqPreprocessor.NCORES)

    parser.add_argument("-n", "--ntrimmed", type=int, default=FastqPreprocessor.NTRIMMED,
                        help='Max number of adapters to trim from each read. Default %i.' % FastqPreprocessor.NTRIMMED)

    parser.add_argument("-q", "--min_bq", type=int, default=FastqPreprocessor.TRIM_QUALITY,
                        help='Min BQ for 3\' trimming. Default %i.' % FastqPreprocessor.TRIM_QUALITY)

    parser.add_argument("-nt", "--no_trim", action="store_true",
                        help='Flag to turn off adapter and base quality trimming.')

    parser.add_argument("-v", "--validate", action="store_true", help='QC FASTQs with FastQC?')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(f1, f2, r1_fiveprime_adapters, r1_threeprime_adapters,
             outdir=".", ncores=FastqPreprocessor.NCORES, min_quality=FastqPreprocessor.TRIM_QUALITY,
             ntrimmed=FastqPreprocessor.NTRIMMED, no_trim=FastqPreprocessor.TRIM_FLAG, validate=False):
    """Runs the FASTQ preprocessor workflow.

    :param str f1: path of the R1 FASTQ
    :param str f2: path of the R2 FASTQ
    :param str r1_fiveprime_adapters: comma-delimited 5' adapters to be trimmed from R1s
    :param str r1_threeprime_adapters: comma-delimited 3' adapters to be trimmed from R1s
    :param str outdir: Output directory to write preprocessed FASTQs to
    :param int ncores: Number of CPU cores to use for cutadapt. Default 0, autodetect.
    :param int min_quality: quality score for quality trimming at the 3' end. Default 15.
    :param int ntrimmed: Max number of adapters to trim from each read
    :param bool no_trim: flag to turn off adapter and 3' base quality trimming. Default False.
    :param bool validate: Validate FASTQs with FastQC?
    :return analysis.read_preprocessor.FastqPreprocessor: read preprocessor object
    """

    fqp = FastqPreprocessor(
        f1=f1, f2=f2, r1_fiveprime_adapters=r1_fiveprime_adapters, r1_threeprime_adapters=r1_threeprime_adapters,
        outdir=outdir, ncores=ncores, trim_bq=min_quality, ntrimmed=ntrimmed, no_trim=no_trim, validate=validate)

    return fqp


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(f1=parsed_args["fastq1"], f2=parsed_args["fastq2"],
             r1_fiveprime_adapters=parsed_args["r1_fiveprime_adapters"],
             r1_threeprime_adapters=parsed_args["r1_threeprime_adapters"],
             outdir=parsed_args["outdir"], ncores=parsed_args["cores"],
             min_quality=parsed_args["min_bq"], ntrimmed=parsed_args["ntrimmed"],
             no_trim=parsed_args["no_trim"], validate=parsed_args["validate"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
