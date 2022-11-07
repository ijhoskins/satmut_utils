#!/usr/bin/env python3
"""Runs variant generation using all possible codon permutations for a transcript."""

import argparse
import logging
import os
import sys

from satmut_utils.definitions import DEFAULT_MUT_SIG, LOG_FORMATTER
from prototype.variant_generator import VariantGenerator
from core_utils.file_utils import replace_extension
from core_utils.string_utils import none_or_str

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
    r"""Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-i", "--trx_id", type=str, required=True,
                        help='Full Ensembl transcript ID (with version) for which to generate variants.')

    parser.add_argument("-g", "--transcript_gff", required=True, type=str,
                        help='GFF file containing transcript metafeatures and exon, CDS, and stop_codon features, '
                             'from 5\' to 3\', regardless of strand.')

    parser.add_argument("-r", "--reference", type=str, required=True, help='Reference FASTA to generate variants in.')

    parser.add_argument("-k", "--gff_reference", type=str, required=True, help='Reference FASTA for the GFF.')

    parser.add_argument("-t", "--targets", type=none_or_str, default="None",
                        help='Optional target BED file. Only variants intersecting the targets will be generated. '
                             'Contig names in the target file should match the contig name in the reference FASTA.')

    parser.add_argument("-o", "--out_vcf", type=none_or_str, default="None",
                        help='Optional output VCF basename. Default use trx_id prefix.')

    parser.add_argument("-d", "--output_dir", type=str, default=VariantGenerator.DEFAULT_OUTDIR,
                        help='Optional output directory. Default current working directory.')

    parser.add_argument("-s", "--mutagenesis_signature", type=str, default=DEFAULT_MUT_SIG,
                        help='Mutagenesis signature. Only variants matching the mutagensis signature will be generated. '
                             'One of {NNN, NNK, NNS}. Default %s.' % DEFAULT_MUT_SIG)

    parser.add_argument("-v", "--var_type", type=str, default=VariantGenerator.DEFAULT_VAR_TYPE,
                        help='Variant types to generate. One of {snp, mnp, total}. Default %s.' %
                             VariantGenerator.DEFAULT_VAR_TYPE)

    parser.add_argument("-a", "--add_haplotypes", action="store_true",
                        help='Flag to add long range haplotypes up to read length.')

    parser.add_argument("-l", "--haplotype_length", type=int, default=VariantGenerator.DEFAULT_HAPLO_LEN,
                        help='Maximum length for which to create haplotypes. Default %i.'
                             % VariantGenerator.DEFAULT_HAPLO_LEN)

    parser.add_argument("-m", "--mnp_bases", type=int, default=VariantGenerator.DEFAULT_MNP_BASES,
                        help='For var_type==mnp or var_type==total, max number of bases in MNPs. Must be either 2 or 3. '
                             'Default %i.' % VariantGenerator.DEFAULT_MNP_BASES)

    parser.add_argument("-y", "--random_seed", type=int, default=9, help='Seed for random variant sampling.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(trx_id, transcript_gff, ref, gff_ref, targets=VariantGenerator.DEFAULT_TARGETS,
             out_vcf=VariantGenerator.DEFAULT_OUTFILE, outdir=VariantGenerator.DEFAULT_OUTDIR, mut_sig=DEFAULT_MUT_SIG,
             var_type=VariantGenerator.DEFAULT_VAR_TYPE, haplotypes=VariantGenerator.DEFAULT_HAPLO,
             haplotype_len=VariantGenerator.DEFAULT_HAPLO_LEN, random_seed=VariantGenerator.DEFAULT_HAPLO_SEED,
             mnp_bases=VariantGenerator.DEFAULT_MNP_BASES):
    r"""Runs the variant generation workflow.

    :param str trx_id: Transcript ID for which to generate variants
    :param str transcript_gff: GFF/GTF file containing transcript metafeatures and exon features, in 5' to 3' order, \
    regardless of strand
    :param str ref: reference FASTA used in alignment/variant calling
    :param str gff_ref: reference FASTA corresponding to the GFF features
    :param str | None targets: optional target feature file. Only variants intersecting the target will be generated.
    :param str | None out_vcf: output VCF name
    :param str outdir: Optional output directory. Default current directory.
    :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNN.
    :param str var_type: one of {"snp", "mnp", "total"}
    :param bool haplotypes: should haplotypes be created with uniform number to codon variants? Default True.
    :param int haplotype_len: max length to create haplotypes. No longer than read length.
    :param int random_seed: seed for variant sampling.
    :param int mnp_bases: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
    :return str: name of the output VCF
    """

    vg = VariantGenerator(
        gff=transcript_gff, ref=ref, gff_ref=gff_ref, mut_sig=mut_sig, haplotypes=haplotypes,
        haplotype_len=haplotype_len, outdir=outdir, random_seed=random_seed)

    out_vcf = vg.workflow(trx_id=trx_id, targets=targets, outfile=out_vcf, var_type=var_type, mnp_bases=mnp_bases)
    return out_vcf


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

    workflow(trx_id=parsed_args["trx_id"], transcript_gff=parsed_args["transcript_gff"], ref=parsed_args["reference"],
             gff_ref=parsed_args["gff_reference"], targets=parsed_args["targets"], out_vcf=parsed_args["out_vcf"],
             outdir=parsed_args["output_dir"],  mut_sig=parsed_args["mutagenesis_signature"],
             var_type=parsed_args["var_type"], haplotypes=parsed_args["add_haplotypes"],
             haplotype_len=parsed_args["haplotype_length"], random_seed=parsed_args["random_seed"],
             mnp_bases=parsed_args["mnp_bases"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
