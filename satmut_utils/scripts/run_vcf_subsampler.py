#!/usr/bin/env/python
"""Runs VCF subsampling."""

import argparse
import logging
import sys

from core_utils.string_utils import none_or_str
from core_utils.vcf_utils import VcfSubsampler

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
    parser.add_argument("-v", "--cf", required=True, type=str, help='BCF/VCF to subsample.')

    parser.add_argument("-n", "--nvars", type=int,
                        help='Number of variants to subsample. If this is greater than the number of variants in the '
                             'input, a warning is raised and no subsampling is performed.')

    parser.add_argument("-r", "--random_seed", type=int, default=VcfSubsampler.DEFAULT_SEED,
                        help='Seed for random qname sampling.')

    parser.add_argument("-d", "--output_dir", type=str, default=".", help='Output directory.')

    parser.add_argument("-o", "--output_vcf", type=none_or_str, default="None", help='Optional output VCF name.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(cf, nvars, output_dir=".", random_seed=VcfSubsampler.DEFAULT_SEED, output_vcf=None):
    """Runs the VCF subsampling workflow.

    :param str cf: BCF or VCF to subsample
    :param int nvars: number of variants requested
    :param str output_dir: optional output directory for the results
    :param int random_seed: seed for variant sampling
    :param str | None output_vcf: output filename
    """

    vs = VcfSubsampler(cf=cf, outdir=output_dir, random_seed=random_seed)
    vs.subsample_variants(nvars=nvars, outfile=output_vcf)


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(cf=parsed_args["cf"], nvars=parsed_args["nvars"], output_dir=parsed_args["output_dir"],
             random_seed=parsed_args["random_seed"], output_vcf=parsed_args["output_vcf"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
