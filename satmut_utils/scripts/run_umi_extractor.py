#!/usr/bin/env/python
"""Runs the workflow for UMI extraction using umi_tools extract."""

import argparse
import logging
import sys

from definitions import AMP_UMI_REGEX
from analysis.read_preprocessor import UMIExtractor
from core_utils.string_utils import none_or_str

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

    parser.add_argument("-f1", "--fastq1", required=True, type=str, help='R1 FASTQ.')

    parser.add_argument("-f2", "--fastq2", required=True, type=str, help='R2 FASTQ.')

    parser.add_argument("-ur", "--umi_regex", type=str, default=AMP_UMI_REGEX,
                        help='UMI regular expression to be passed to umi_tools extract command.')

    parser.add_argument("-pf", "--primer_fasta", type=none_or_str,
                        help='A primer FASTA file may optionally be provided to append originating primers to read names. '
                             'Useful for RACE-like (e.g. AMP) libraries to prohibit R2 merging. That is, without '
                             'this flag, tiled R2s sharing the same R1 UMI-position will be merged into '
                             'contigs during consensus deduplication.')

    parser.add_argument("-nm", "--primer_nm_allowance", type=int, default=UMIExtractor.PRIMER_NM_ALLOW,
                        help='If passing a primer FASTA file with -pf, match primers in R2 with up to the given '
                             'number of edit operations.')

    parser.add_argument("-o", "--outdir", type=str, default=".", help='Output directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(r1_fastq, r2_fastq, umi_regex=AMP_UMI_REGEX, primer_fasta=UMIExtractor.PRIMER_FASTA,
             primer_nm_allow=UMIExtractor.PRIMER_NM_ALLOW, outdir="."):
    r"""Runs the consensus deduplication workflow.

    :param str r1_fastq: R1 FASTQ
    :param str r2_fastq: R2 FASTQ
    :param str umi_regex: regex for matching the UMIs (see umi_tools for docs)
    :param str | None primer_fasta: primer FASTA
    :param int primer_nm_allow: Edit distance allowance for primer matching Default 3.
    :param str outdir: Output directory
    """

    UMIExtractor(r1_fastq=r1_fastq, r2_fastq=r2_fastq, umi_regex=umi_regex,
                 primer_fasta=primer_fasta, primer_nm_allow=primer_nm_allow, outdir=outdir)


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(r1_fastq=parsed_args["fastq1"], r2_fastq=parsed_args["fastq2"], umi_regex=parsed_args["umi_regex"],
             primer_fasta=parsed_args["primer_fasta"], primer_nm_allow=parsed_args["primer_nm_allowance"],
             outdir=parsed_args["outdir"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
