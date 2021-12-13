#!/usr/bin/env python3
"""Runs VCF subsampling."""

import argparse
import logging
import os
import sys

from core_utils.file_utils import replace_extension
from core_utils.string_utils import none_or_str, none_or_int
from core_utils.vcf_utils import VcfSubsampler
from definitions import LOG_FORMATTER

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

LOGFILE = replace_extension(os.path.basename(__file__), "log")
logger = logging.getLogger(__name__)


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-v", "--cf", required=True, type=str,
                        help='VCF/BCF to subsample. May be gzipped. '
                             'Any multiallelic records will be split into records for each ALT.')

    parser.add_argument("-n", "--nvars", required=True, type=int,
                        help='Number of variants to subsample. If this is greater than the number of variants in the '
                             'input VCF/BCF, a warning is raised and no subsampling is performed.')

    parser.add_argument("-b", "--nbases", type=none_or_int, default="None",
                        help='Subsample SNPs and MNP variants to achieve this number of mismatched bases. '
                             'Must be provided along with --nvars.')

    parser.add_argument("-r", "--random_seed", type=int, default=VcfSubsampler.DEFAULT_SEED,
                        help='Seed for random variant subsampling. Default %i.' % VcfSubsampler.DEFAULT_SEED)

    parser.add_argument("-s", "--snp_proportion", type=float,
                        help='SNP weight to determine number of SNPs versus MNPs to subsample. di- and tri-nt MNPs '
                             'are sampled uniformly from the complement of this proportion. Default %f.'
                             % VcfSubsampler.DEFAULT_SNP_PROP)

    parser.add_argument("-d", "--output_dir", type=str, default=VcfSubsampler.DEFAULT_OUTDIR,
                        help='Optional output directory. Default current working directory.')

    parser.add_argument("-o", "--output_vcf", type=none_or_str, default="None",
                        help='Optional output VCF name. Default use basename of input file.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(cf, nvars, nbases=VcfSubsampler.DEFAULT_NBASES, output_dir=VcfSubsampler.DEFAULT_OUTDIR,
             random_seed=VcfSubsampler.DEFAULT_SEED, snp_prop=VcfSubsampler.DEFAULT_SNP_PROP,
             output_vcf=VcfSubsampler.DEFAULT_OUTFILE):
    """Runs the VCF subsampling workflow.

    :param str cf: VCF or BCF to subsample. May be gzipped.
    :param int nvars: number of variants requested
    :param int | None nbases: number of mutant bases requested
    :param str output_dir: optional output directory for the results
    :param int random_seed: seed for variant sampling
    :param float snp_prop: SNP proportion; di-nt and tri-nt MNPs will be uniformly drawn from the complement
    :param str | None output_vcf: output filename
    """

    vs = VcfSubsampler(cf=cf, outdir=output_dir, random_seed=random_seed, snp_prop=snp_prop)

    if nbases is not None:
        vs.subsample_nbases(nvars=nvars, nbases=nbases, outfile=output_vcf)
    else:
        vs.subsample_variants(nvars=nvars, outfile=output_vcf)


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

    workflow(cf=parsed_args["cf"], nvars=parsed_args["nvars"], output_dir=parsed_args["output_dir"],
             random_seed=parsed_args["random_seed"], output_vcf=parsed_args["output_vcf"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
