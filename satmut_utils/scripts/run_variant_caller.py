#!/usr/bin/env/python
"""Runs the targeted variant caller for paired-end data."""

import argparse
import logging
import sys

from analysis.seq_utils import UNKNOWN_BASE
from analysis.variant_caller import VariantCaller

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
    parser.add_argument("-a", "--alignments", required=True, type=str, help='SAM/BAM file to call variants from.')

    parser.add_argument("-t", "--targets", required=True, type=str, help='Target feature file, e.g. BED, GFF, or GTF.')

    parser.add_argument("-r", "--reference", required=True, type=str, help='Reference FASTA for alignment.')

    parser.add_argument("-g", "--transcript_gff", required=True, type=str,
                        help='GFF file containing transcript metafeatures and exon features. The GFF must be '
                             'from 5\' to 3\', regardless of strand.')

    parser.add_argument("-gr", "--gff_reference", required=True, type=str,
                        help='Reference FASTA for the GFF. Needed for custom (e.g. vector) GFF files using the '
                             'same reference as alignment (--reference).')

    parser.add_argument("-p", "--primers", type=str, help='Primer feature file, e.g. BED, GFF, or GTF. Must have strand.')

    parser.add_argument("-o", "--outdir", type=str, help='Output directory.')

    parser.add_argument("-n", "--nthreads", type=int, default=VariantCaller.DEFAULT_NTHREADS,
                        help='Number of threads to use for SAM/BAM preprocessing manipulations. Default %i.' %
                             VariantCaller.DEFAULT_NTHREADS)

    parser.add_argument("-b", "--min_bq", type=int, default=VariantCaller.VARIANT_CALL_MIN_BQ,
                        help='Min BQ to consider a read for variant calling. Default %i.' %
                             VariantCaller.VARIANT_CALL_MIN_BQ)

    parser.add_argument("-e", "--max_nm", type=int, default=VariantCaller.VARIANT_CALL_MAX_NM,
                        help='Max edit distance to consider a read for variant calling. Default %i.' %
                             VariantCaller.VARIANT_CALL_MAX_NM)

    parser.add_argument("-s", "--min_supporting", type=int, default=VariantCaller.VARIANT_CALL_MIN_DP,
                        help='Min R1-R2 concordant counts for establishing candidate variants. Default %i.' %
                             VariantCaller.VARIANT_CALL_MIN_DP)

    parser.add_argument("-w", "--max_mnp_window", type=int, default=VariantCaller.VARIANT_CALL_MAX_MNP_WINDOW,
                        help='Max window to search for a MNP and merge phased SNPs. Must be >= 3. '
                             'Note any phased SNP within this window will be merged into a MNP, and its component SNPs '
                             'may not be called. Default %i.' % VariantCaller.VARIANT_CALL_MAX_MNP_WINDOW)

    parser.add_argument("-x", "--include_n", action="store_true",
                        help='Flag indicating to also call ALTs equal to, or containing, %s. '
                             'Potentially useful for training error models.' % UNKNOWN_BASE)

    parser.add_argument("-ms", "--mutagenesis_signature", type=str, default=VariantCaller.VARIANT_CALL_MUT_SIG,
                        help='Mutagenesis signature. One of NNN, NNK, or NNS.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(alignments, targets, ref, transcript_gff, gff_reference, outdir=None, primers=None,
             nthreads=VariantCaller.DEFAULT_NTHREADS, max_nm=VariantCaller.VARIANT_CALL_MAX_NM,
             min_bq=VariantCaller.VARIANT_CALL_MIN_BQ, min_supporting_qnames=VariantCaller.VARIANT_CALL_MIN_DP,
             max_mnp_window=VariantCaller.VARIANT_CALL_MAX_MNP_WINDOW,
             include_n=not VariantCaller.VARIANT_CALL_EXCLUDE_N,
             mut_sig=VariantCaller.VARIANT_CALL_MUT_SIG):
    r"""Runs the variant caller workflow.

    :param str alignments: SAM/BAM file to call variants in
    :param str targets: BED, GFF, GTF file containing targets to call variants in
    :param str ref: path to reference FASTA used in alignment. Must be samtools faidx indexed.
    :param str transcript_gff: GFF/GTF file containing transcript metafeatures and exon features, in 5' to 3' order, \
    regardless of strand. Ordering is essential.
    :param str gff_reference: reference FASTA corresponding to the GFF features
    :param str | None outdir: Optional output dir for the results
    :param str | None primers: BED or GFF file containing primers. For masking synthetic sequences for \
    accurate AF of variants overlapping primer regions. This feature file must contain the strand of the primer. \
    Set to None for no masking.
    :param int nthreads: number threads to use for SAM/BAM file manipulations. Default 0 (autodetect)
    :param int min_bq: min base quality; should be >= 1 such that masked primers do not contribute to depth
    :param int max_nm: max edit distance to consider a read for variant calling
    :param int min_supporting_qnames: min number of fragments with R1-R2 concordant coverage for which to keep a variant
    :param int max_mnp_window: max number of consecutive nucleotides to search for di-nt MNPs; must be >= 3
    :param bool include_n: include variant calls to an "N"
    :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNK.
    :return tuple: (VCF, BED) filepaths
    """

    vc = VariantCaller(
        am=alignments, targets=targets, ref=ref, trx_gff=transcript_gff, gff_ref=gff_reference,
        primers=primers, output_dir=outdir, nthreads=nthreads, exclude_n=not include_n)

    output_vcf, output_bed = vc.workflow(
        min_bq=min_bq, max_nm=max_nm, min_supporting_qnames=min_supporting_qnames,
        max_mnp_window=max_mnp_window, mut_sig=mut_sig)

    return output_vcf, output_bed


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(alignments=parsed_args["alignments"], targets=parsed_args["targets"], ref=parsed_args["reference"],
             transcript_gff=parsed_args["transcript_gff"], gff_reference=parsed_args["gff_reference"],
             outdir=parsed_args["outdir"], primers=parsed_args["primers"], nthreads=parsed_args["nthreads"],
             min_bq=parsed_args["min_bq"], max_nm=parsed_args["max_nm"],
             min_supporting_qnames=parsed_args["min_supporting"], max_mnp_window=parsed_args["max_mnp_window"],
             include_n=parsed_args["include_n"], mut_sig=parsed_args["mutagenesis_signature"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
