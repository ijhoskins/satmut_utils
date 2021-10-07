#!/usr/bin/env/python
"""Runs in silico variant induction."""

import argparse
import datetime
import logging
import os
import sys

import core_utils.file_utils as fu
from core_utils.string_utils import none_or_str
import analysis.read_editor as ri

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.2"
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
    parser.add_argument("-a", "--alignments", type=str, help='SAM/BAM file containing alignments to edit.')

    parser.add_argument("-v", "--vcf", type=str,
                        help='VCF file specifying variants to edit. Should have an AF INFO field specifying a '
                             'fraction of fragments to edit into.')

    parser.add_argument("-r", "--reference", type=str, required=True, help='Corresponding reference FASTA.')

    parser.add_argument("-p", "--primers", type=none_or_str, default="None",
                        help='Primer feature file, e.g. BED, GFF, or GTF. Must have strand. '
                             'Without, no primer masking occurs.')

    parser.add_argument("-d", "--output_dir", type=str, default=".", help='Output directory for FASTQs and BAM.')

    parser.add_argument("-o", "--output_prefix", type=none_or_str,
                        default=".".join([datetime.datetime.now().strftime("%d%b%Y.%I.%M%p"), "read_inducer"]),
                        help='Output prefix for FASTQs and BAM.')

    parser.add_argument("-s", "--single_end", action="store_true", help='Flag indicating the input is single data.')

    parser.add_argument("-rs", "--random_seed", type=int, default=9, help='Seed for random qname sampling.')

    parser.add_argument("-n", "--no_realignment", action="store_true",
                        help='Flag indicating to not realign the induced FASTQs. Cuts runtime.')

    parser.add_argument("-f", "--filter_induced", action="store_true",
                        help='Flag indicating to filter and output induced reads/read pairs as BAM.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(am, variants, reference, primers=None, output_dir=".", output_prefix=None, single_end=False, realign=True,
             filter_induced=True, random_seed=9):
    """Runs the ReadInducer workflow.

    :param str am: SAM/BAM file to induce into
    :param str variants: VCF file specifying variants to induce
    :param str reference: reference FASTA.
    :param str | None primers: feature file of primer locations for read masking and primer detection
    :param str output_dir: Optional output directory to store generated FASTQs and BAM
    :param str | None output_prefix: Optional output prefix for the FASTQ(s) and BAM; if None, use same prefix as BAM.
    :param bool single_end: does the input SAM/BAM contain single end reads? Default False, paired end reads.
    :param bool realign: realign the induced FASTQs? Needed for visualization of the induced BAM.
    :param bool filter_induced: filter the induced BAM for those reads/read pairs that were induced? Default True.
    :param int random_seed: seed for random qname sampling
    :return tuple: (str | None, str, str | None) paths of the induced BAM, R1 FASTQ, R2 FASTQ
    """

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    out_prefix = output_prefix
    if output_prefix is None:
        out_prefix = fu.remove_extension(os.path.basename(am))

    # Run the induction workflow
    output_bam, zipped_r1_fastq, zipped_r2_fastq = ri.ReadEditor.workflow(
        am=am, vcf=variants, ref=reference, primers=primers, output_dir=output_dir, output_prefix=out_prefix,
        single_end=single_end, realign=realign, filter_edited=filter_induced, random_seed=random_seed)

    return output_bam, zipped_r1_fastq, zipped_r2_fastq


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(am=parsed_args["alignments"], variants=parsed_args["vcf"], reference=parsed_args["reference"],
             primers=parsed_args["primers"], output_dir=parsed_args["output_dir"], output_prefix=parsed_args["output_prefix"],
             single_end=parsed_args["single_end"], realign=not parsed_args["no_realignment"],
             filter_induced=parsed_args["filter_induced"], random_seed=parsed_args["random_seed"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
