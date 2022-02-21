#!/usr/bin/env/python
"""Runs FASTQ preprocessing (adapter, quality trimming)."""

import argparse
import logging
import os
import sys

from core_utils.file_utils import replace_extension
from core_utils.string_utils import none_or_str
from satmut_utils.definitions import LOG_FORMATTER
from analysis.read_preprocessor import FastqPreprocessor

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
    parser.add_argument("-1", "--fastq1", type=str, required=True, help='R1 FASTQ')

    parser.add_argument("-2", "--fastq2", type=str, required=True, help='R2 FASTQ')

    parser.add_argument("--r1_fiveprime_adapters", type=none_or_str, default="None",
                        help='Comma-delimited R1 5\' adapters, or None if no R1 5\' adapters exist.')

    parser.add_argument("--r1_threeprime_adapters", type=none_or_str, default="None",
                        help='Comma-delimited R1 3\' adapters, or None if no R1 3\' adapters exist.')

    parser.add_argument("--r2_fiveprime_adapters", type=none_or_str, default="None",
                        help='Comma-delimited R2 5\' adapters, or None if no R2 5\' adapters exist.')

    parser.add_argument("--r2_threeprime_adapters", type=none_or_str, default="None",
                        help='Comma-delimited R2 3\' adapters, or None if no R2 3\' adapters exist.')

    parser.add_argument("-o", "--output_dir", type=str, default=FastqPreprocessor.OUTDIR,
                        help='Optional output directory. Default current working directory.')

    parser.add_argument("-c", "--cores", type=int, default=FastqPreprocessor.NCORES,
                        help='Number of cores to use for cutadapt. Default %i.' % FastqPreprocessor.NCORES)

    parser.add_argument("-n", "--ntrimmed", type=int, default=FastqPreprocessor.NTRIMMED,
                        help='Max number of adapters to trim from each read. Useful for trimming terminal tiles '
                             'with vector-transgene alignment. Default %i.' % FastqPreprocessor.NTRIMMED)

    parser.add_argument("-l", "--overlap_length", type=int, default=FastqPreprocessor.OVERLAP_LEN,
                        help='Number of read bases overlapping the adapter sequence(s) to consider for cutadapt trimming. '
                             'Default %i.' % FastqPreprocessor.OVERLAP_LEN)

    parser.add_argument("-b", "--trim_bq", type=int, default=FastqPreprocessor.TRIM_QUALITY,
                        help='Base quality for cutadapt 3\' trimming. Default %i.' % FastqPreprocessor.TRIM_QUALITY)

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(f1, f2, r1_fiveprime_adapters=FastqPreprocessor.DEFAULT_ADAPTER,
             r1_threeprime_adapters=FastqPreprocessor.DEFAULT_ADAPTER,
             r2_fiveprime_adapters=FastqPreprocessor.DEFAULT_ADAPTER,
             r2_threeprime_adapters=FastqPreprocessor.DEFAULT_ADAPTER,
             outdir=FastqPreprocessor.OUTDIR,
             ncores=FastqPreprocessor.NCORES, ntrimmed=FastqPreprocessor.NTRIMMED,
             overlap_len=FastqPreprocessor.OVERLAP_LEN, trim_bq=FastqPreprocessor.TRIM_QUALITY):
    """Runs the FASTQ preprocessor workflow.

    :param str f1: path of the R1 FASTQ
    :param str f2: path of the R2 FASTQ
    :param str | None r1_fiveprime_adapters: comma-delimited 5' adapters to trim from R1. Default None.
    :param str | None r1_threeprime_adapters: comma-delimited 3' adapters to trim from R1. Default None.
    :param str | None r2_fiveprime_adapters: comma-delimited 5' adapters to trim from R2. Default None.
    :param str | None r2_threeprime_adapters: comma-delimited 3' adapters to trim from R2. Default None.
    :param str outdir: optional output dir for the results
    :param int ncores: Number of CPU cores to use for cutadapt. Default 0 (autodetect).
    :param int ntrimmed: Max number of adapters to trim from each read. Default 3.
    :param int overlap_len: number of bases to match in read to trim. Default 8.
    :param int trim_bq: quality score for cutadapt quality trimming at the 3' end. Default 15.
    :return tuple: (str, str) R1 and R2 trimmed FASTQ filepaths
    """

    fqp = FastqPreprocessor(
        f1=f1, f2=f2, r1_fiveprime_adapters=r1_fiveprime_adapters, r1_threeprime_adapters=r1_threeprime_adapters,
        r2_fiveprime_adapters=r2_fiveprime_adapters, r2_threeprime_adapters=r2_threeprime_adapters, outdir=outdir,
        ncores=ncores, trim_bq=trim_bq, ntrimmed=ntrimmed, overlap_len=overlap_len, no_trim=False)

    return fqp.trimmed_f1, fqp.trimmed_f2


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

    workflow(f1=parsed_args["fastq1"], f2=parsed_args["fastq2"],
             r1_fiveprime_adapters=parsed_args["r1_fiveprime_adapters"],
             r1_threeprime_adapters=parsed_args["r1_threeprime_adapters"],
             r2_fiveprime_adapters=parsed_args["r2_fiveprime_adapters"],
             r2_threeprime_adapters=parsed_args["r2_threeprime_adapters"],
             outdir=parsed_args["output_dir"], ncores=parsed_args["cores"], ntrimmed=parsed_args["ntrimmed"],
             overlap_len=parsed_args["overlap_length"], trim_bq=parsed_args["trim_bq"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
