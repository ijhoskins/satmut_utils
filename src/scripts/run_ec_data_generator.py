#!/usr/bin/env python3
"""Runs the error correction data generation workflow for training models."""

import argparse
import datetime
import logging
import os
import sys

from analysis.read_editor import ReadEditor
from core_utils.file_utils import replace_extension
from core_utils.string_utils import none_or_str
from definitions import DEFAULT_MUT_SIG, LOG_FORMATTER
from prototype.variant_generator import VariantGenerator
from prototype.error_correction import ErrorCorrectionDataGenerator

__author__ = "Ian Hoskins"
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
    parser.add_argument("-n", "--negative_summary", type=str, required=True,
                        help='.vcf.summary.txt file for the negative control alignments (--negative_bam).')

    parser.add_argument("-m", "--mutant_summary", type=str, required=True,
                        help='.vcf.summary.txt file for the mutant sample. Used to estimate frequency parameters '
                             'for each variant type (SNP, di-nt MNP, tri-nt MNP).')

    parser.add_argument("-b", "--negative_bam", type=str, required=True, help='BAM file for the negative control.')

    parser.add_argument("-i", "--trx_id", type=str, required=True, help='Transcript ID for which to generate variants.')

    parser.add_argument("-r", "--reference", type=str, required=True, help='Corresponding reference FASTA.')

    parser.add_argument("-g", "--transcript_gff", type=str, required=True,
                        help='GFF file containing transcript metafeatures and exon, CDS, start_codon, and stop_codon '
                             'features, from 5\' to 3\', regardless of strand.')

    parser.add_argument("-t", "--targets", type=none_or_str, default="None",
                        help='Optional target BED file. Only variants intersecting the targets will be generated. '
                             'Contig names in the target file should match the contig name in the reference FASTA.')

    parser.add_argument("-z", "--race_like", action="store_true",
                        help='Flag for data generated from a RACE-like chemistry with R1s starting at variable ends '
                             'and R2s starting at primers.')

    parser.add_argument("-p", "--primers", type=none_or_str, default="None",
                        help='Primer BED file. Must have strand field. Default no masking.')

    parser.add_argument("-d", "--output_dir", type=str, default=".",
                        help='Optional output directory. Default current working directory.')

    parser.add_argument("-s", "--mutagenesis_signature", type=str, default=DEFAULT_MUT_SIG,
                        help='Mutagenesis signature. Only variants matching the signature will be generated. '
                             'One of {NNN, NNK, NNS}. Default %s.' % DEFAULT_MUT_SIG)

    parser.add_argument("-o", "--output_prefix", type=none_or_str,
                        default=".".join([datetime.datetime.now().strftime("%d%b%Y.%I.%M%p"), "ec_data"]),
                        help='Output prefix for VCFs.')

    parser.add_argument("-v", "--var_type", type=str, default=VariantGenerator.DEFAULT_VAR_TYPE,
                        help='Variant types to generate. One of {snp, mnp, total}. Default %s.' %
                             VariantGenerator.DEFAULT_VAR_TYPE)

    parser.add_argument("-a", "--add_haplotypes", action="store_true",
                        help='Flag to add long range haplotypes up to read length.')

    parser.add_argument("-l", "--haplotype_length", type=int, default=VariantGenerator.DEFAULT_HAPLO_LEN,
                        help='Maximum length for which to create haplotypes.')

    parser.add_argument("-x", "--mnp_bases", type=int, default=VariantGenerator.DEFAULT_MNP_BASES,
                        help='For var_type==mnp or var_type==total, max bases for MNPs. Must be either 2 or 3. '
                             'Default %s.' % VariantGenerator.DEFAULT_MNP_BASES)

    parser.add_argument("-y", "--random_seed", type=int, default=VariantGenerator.DEFAULT_HAPLO_SEED,
                        help='Seed for random variant sampling.')

    parser.add_argument("-e", "--edit_buffer", type=int, default=ReadEditor.DEFAULT_BUFFER,
                        help='Buffer +/- the edit coordinate position(s) to check for pre-existing errors in reads. '
                             'Used to restrict phasing of true variants with errors, leading to false negatives. '
                             'Default %i' % ReadEditor.DEFAULT_BUFFER)

    parser.add_argument("-f", "--force_edit", action="store_true",
                        help='Flag to force editing of variants in the event of invalid variant configurations '
                             '(AF sum > 1). Editing of all variants is not ensured if this flag is provided, '
                             'especially for single-amplicon (PCR tile) alignments. If alignments contain many tiles, '
                             'this option enables more variants to be simulated, as they are diffuse in target space.')

    parser.add_argument("-j", "--nthreads", type=int, default=ErrorCorrectionDataGenerator.DEFAULT_NTHREADS,
                        help='Number of threads to use for BAM sorting operations. '
                             'Default %i, autodetect.' % ErrorCorrectionDataGenerator.DEFAULT_NTHREADS)

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(negative_summary, mutant_summary, negative_bam, trx_id, ref, gff, targets=VariantGenerator.DEFAULT_TARGETS,
             race_like=ErrorCorrectionDataGenerator.DEFAULT_RACE_LIKE, primers=ErrorCorrectionDataGenerator.DEFAULT_PRIMERS,
             output_dir=ErrorCorrectionDataGenerator.DEFAULT_OUTDIR, mut_sig=DEFAULT_MUT_SIG,
             output_prefix=ErrorCorrectionDataGenerator.DEFAULT_PREFIX, var_type=VariantGenerator.DEFAULT_VAR_TYPE,
             haplotypes=VariantGenerator.DEFAULT_HAPLO, haplotype_len=VariantGenerator.DEFAULT_HAPLO_LEN,
             mnp_bases=VariantGenerator.DEFAULT_MNP_BASES, random_seed=VariantGenerator.DEFAULT_HAPLO_SEED,
             buffer=ReadEditor.DEFAULT_BUFFER, force_edit=ReadEditor.DEFAULT_FORCE,
             nthreads=ErrorCorrectionDataGenerator.DEFAULT_NTHREADS):
    r"""Runs the error correction data generation workflow.

    :param str negative_summary: vcf.summary.txt file for the negative control library
    :param str mutant_summary: vcf.summary.txt file for a mutant library, used for estimating AF parameters
    :param str negative_bam: negative control BAM to edit into; should have alignments to trx_id only
    :param str trx_id: transcript ID to generate variants for
    :param str ref: reference FASTA for trx_id
    :param str gff: transcript GFF for trx_id
    :param str | None targets: optional target feature file. Only variants intersecting the target will be generated.
    :param bool race_like: is the data produced by RACE-like (e.g. AMP) data? Default False.
    :param str | None primers: feature file of primer locations for read masking and primer detection
    :param str output_dir: Optional output directory to write ReadEditor output files.
    :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNN.
    :param str | None output_prefix: Optional output prefix for the FASTQ(s) and BAM; if None, use same prefix as VCF
    :param str var_type: one of {"snp", "mnp", "total"}
    :param bool haplotypes: should haplotypes be created with uniform number to codon variants? Default True.
    :param int haplotype_len: max length to create haplotypes. No longer than read length.
    :param int mnp_bases: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
    :param int random_seed: seed for variant and qname random sampling
    :param int buffer: buffer about the edit span (position + REF len) to ensure lack of error before editing. Default 3.
    :param bool force_edit: flag to attempt editing of variants despite a NonconfiguredVariant exception.
    :param int nthreads: Number of threads to use for BAM operations. Default 0 (autodetect).
    :return tuple: (str, str | None, str, str | None) paths of the truth VCF, induced BAM, R1 FASTQ, and R2 FASTQ
    """

    ecdg = ErrorCorrectionDataGenerator(
        negative_summary=negative_summary, mutant_summary=mutant_summary, negative_bam=negative_bam,
        race_like=race_like, ref=ref, gff=gff, primers=primers, outdir=output_dir, nthreads=nthreads)

    truth_vcf, output_bam, zipped_r1_fastq, zipped_r2_fastq = ecdg.workflow(
        trx_id=trx_id, targets=targets, mut_sig=mut_sig, var_type=var_type, mnp_bases=mnp_bases,
        output_prefix=output_prefix, haplotypes=haplotypes, haplotype_len=haplotype_len, random_seed=random_seed,
        buffer=buffer, force_edit=force_edit)

    return truth_vcf, output_bam, zipped_r1_fastq, zipped_r2_fastq


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

    workflow(negative_summary=parsed_args["negative_summary"], mutant_summary=parsed_args["mutant_summary"],
             negative_bam=parsed_args["negative_bam"], trx_id=parsed_args["trx_id"],
             ref=parsed_args["reference"], gff=parsed_args["transcript_gff"],
             targets=parsed_args["targets"], race_like=parsed_args["race_like"],
             primers=parsed_args["primers"], output_dir=parsed_args["output_dir"],
             mut_sig=parsed_args["mutagenesis_signature"], output_prefix=parsed_args["output_prefix"],
             var_type=parsed_args["var_type"], haplotypes=parsed_args["add_haplotypes"],
             haplotype_len=parsed_args["haplotype_length"], mnp_bases=parsed_args["mnp_bases"],
             random_seed=parsed_args["random_seed"], buffer=parsed_args["edit_buffer"],
             force_edit=parsed_args["force_edit"], nthreads=parsed_args["nthreads"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
