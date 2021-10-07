#!/usr/bin/env/python
"""Runs satmut_utils."""

import argparse
import datetime
import logging
import os
import sys
import tempfile

from analysis.read_preprocessor import FastqPreprocessor, UMIExtractor, ReadGrouper, \
    ConsensusDeduplicatorPreprocessor, ConsensusDeduplicator, ReadMasker

import analysis.read_editor as ri
from analysis.references import get_ensembl_references
from analysis.seq_utils import UNKNOWN_BASE
from analysis.variant_caller import VariantCaller

import core_utils.file_utils as fu
from core_utils.string_utils import none_or_str

from definitions import AMP_UMI_REGEX, GRCH38_FASTA

from scripts.run_bowtie2_aligner import workflow as baw


__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "1.0"
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

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    # Common satmut_utils options
    parser = argparse.ArgumentParser(description="%s arguments" % __file__)

    parser.add_argument("-i", "--ensembl_id", type=str,
                        help='Ensembl gene (ENSG) or transcript (ENST) ID to use for a reference.')

    parser.add_argument("-x", "--reference_dir", type=str, default="./references",
                        help='Directory containing curated reference files.')

    parser.add_argument("-r", "--reference", type=str, help='Reference FASTA for alignment.')

    parser.add_argument("-g", "--transcript_gff", type=str,
                        help='GFF file containing transcript metafeatures and exon features. The GFF must be from 5\' '
                             'to 3\', regardless of strand, and contain transcript, exon, CDS, and stop_codon features.')

    parser.add_argument("-p", "--primers", type=none_or_str,
                        help='Optional primer feature file, e.g. BED, GFF, or GTF. Must have strand field. '
                             'Recommended for removing synthetic error from alignments with base quality masking.')

    parser.add_argument("-o", "--outdir", type=str, default=".",
                        help='Optional output directory. Default current working directiory.')

    # Subcommands
    subparsers = parser.add_subparsers(help='sub-command help', title='subcommands')

    # sim subcommand
    parser_sim = subparsers.add_parser('sim', help='sim help')
    parser_sim.set_defaults(func=sim_workflow)

    parser_sim.add_argument("-a", "--alignments", required=True, type=str,
                            help='SAM/BAM file containing alignments to edit.')

    parser_sim.add_argument("-v", "--vcf",  required=True, type=str,
                        help='VCF file specifying variants to edit. Should have an AF INFO field specifying a '
                             'fraction/proportion of fragments to edit into.')

    parser_sim.add_argument("-s", "--random_seed", type=int, default=9, help='Seed for random read sampling.')

    parser_sim.add_argument("-n", "--no_realignment", action="store_true",
                        help='Flag indicating to not realign the edited FASTQs. Cuts runtime.')

    parser_sim.add_argument("-f", "--filter_edited", action="store_true",
                            help='Flag indicating to filter and output edited reads as BAM.')

    # call subcommand
    parser_call = subparsers.add_parser('call', help='call help')
    parser_call.set_defaults(func=call_workflow)

    parser_call.add_argument("-1", "--fastq1", required=True, type=str, help='R1 FASTQ.')

    parser_call.add_argument("-2", "--fastq2", required=True, type=str, help='R2 FASTQ.')

    parser_call.add_argument("-5", "--r1_fiveprime_adapters", required=True, type=str,
                             help='Comma-delimited R1 5\' adapters.')

    parser_call.add_argument("-3", "--r1_threeprime_adapters", required=True, type=str,
                             help='Comma-delimited R1 3\' adapters.')

    parser_call.add_argument("-t", "--targets", type=str,
                             help='Optional target BED file. Only variants intersecting the targets will be reported.')

    parser_call.add_argument("-k", "--gff_reference", type=str, help='Reference FASTA for the GFF.')

    parser_call.add_argument("-d", "--consensus_deduplicate", action="store_true",
                             help='Flag indicating consensus reads should be generated during deduplication.')

    parser_call.add_argument("-u", "--umi_regex", type=str, default=AMP_UMI_REGEX,
                             help='UMI regular expression to be passed to umi_tools extract command.')

    parser_call.add_argument("-c", "--contig_del_threshold", type=int, default=ConsensusDeduplicator.CONTIG_DEL_THRESH,
                             help='If -cd (consensus deduplicate), convert deletions to N that are larger than the threshold. '
                                  'For RACE-like (e.g. AMP) data, a small number of R2s may share the same R1 UMI-position '
                                  'but align to different coordinates, and will be merged, possibly with large gaps.')

    parser_call.add_argument("-f", "--primer_fa", type=none_or_str,
                             help='If -cd, this may optionally be set to append originating primers to read names.'
                                  'Useful for RACE-like (e.g. AMP) libraries to prohibit R2 merging. That is, without '
                                  'this flag, tiled R2s sharing the same R1 will be merged into '
                                  'contigs during consensus deduplication (-cd).')

    parser_call.add_argument("-a", "--primer_nm_allowance", type=int, default=UMIExtractor.PRIMER_NM_ALLOW,
                             help='If -pf, match primers in R2 with up to this many edit operations.')

    parser_call.add_argument("-q", "--min_bq", type=int, default=VariantCaller.VARIANT_CALL_MIN_BQ,
                             help='Min base quality to consider a read for variant calling. Default %i.' %
                                  VariantCaller.VARIANT_CALL_MIN_BQ)

    parser_call.add_argument("-e", "--max_nm", type=int, default=VariantCaller.VARIANT_CALL_MAX_NM,
                             help='Max edit distance to consider a read for variant calling. Default %i.' %
                                  VariantCaller.VARIANT_CALL_MAX_NM)

    parser_call.add_argument("-m", "--min_supporting", type=int, default=VariantCaller.VARIANT_CALL_MIN_DP,
                             help='Min R1-R2 concordant counts for establishing candidate variants. Default %i.' %
                                  VariantCaller.VARIANT_CALL_MIN_DP)

    parser_call.add_argument("-w", "--max_mnp_window", type=int, default=VariantCaller.VARIANT_CALL_MAX_MNP_WINDOW,
                             help='Max window to search for a MNP and merge phased SNPs. Must be <= 3. '
                                  'Any phased SNP within this window will be merged into a MNP, and its component SNPs '
                                  'will not be called. Default %i.' % VariantCaller.VARIANT_CALL_MAX_MNP_WINDOW)

    parser_call.add_argument("-z", "--include_n", action="store_true",
                             help='Flag indicating to also call ALTs equal to, or containing, the unknown base call %s. '
                                  'Potentially useful for training error models.' % UNKNOWN_BASE)

    parser_call.add_argument("-j", "--nthreads", type=int, default=FastqPreprocessor.NCORES,
                             help='Number of threads to use for bowtie2 alignment and BAM operations. '
                                  'Default %i, autodetect.' % FastqPreprocessor.NCORES)

    parser_call.add_argument("-n", "--ntrimmed", type=int, default=FastqPreprocessor.NTRIMMED,
                             help='Max number of adapters to trim from each read. Default %i.'
                                  % FastqPreprocessor.NTRIMMED)

    parser_call.add_argument("-b", "--trim_bq", type=int, default=FastqPreprocessor.TRIM_QUALITY,
                             help='Base quality for 3\' trimming. Default %i.' % FastqPreprocessor.TRIM_QUALITY)

    parser_call.add_argument("-v", "--omit_trim", action="store_true",
                             help='Flag to turn off adapter and 3\' base quality trimming. Useful for simulated data that '
                                  'has already had adapters and bases trimmed.')

    parser_call.add_argument("-s", "--mutagenesis_signature", type=str, default=VariantCaller.VARIANT_CALL_MUT_SIG,
                             help='Mutagenesis signature. One of NNN, NNK, or NNS.')

    parsed_args = parser.parse_args(args)
    return parsed_args


def get_sim_reference(reference_dir, ensembl_id, ref, outdir=ri.ReadEditor.DEFAULT_REFERENCE_DIR):
    """Get and build index files for the sim workflow reference FASTA.

    :param str reference_dir: directory containing curated APPRIS reference files
    :param str ensembl_id: Ensembl gene or transcript ID, with version number
    :param str ref: path to reference FASTA used in alignment. Must be bowtie2 FM-index and samtools faidx indexed.
    :param str outdir: optional output dir for the reference files
    :return str: reference FASTA filepath
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Determine if the provided Ensembl ID is found in the curated APPRIS references
    ref_fa = ref
    if ensembl_id is not None:
        if ref is not None:
            raise RuntimeError("Both an Ensembl ID and a reference FASTA were provided. Please choose one.")

        ref_fa, _ = get_ensembl_references(reference_dir=reference_dir, ensembl_id=ensembl_id, outdir=outdir)

    return ref_fa


def get_call_references(reference_dir, ensembl_id, ref, transcript_gff, gff_reference,
                        outdir=VariantCaller.VARIANT_CALL_REFERENCE_DIR):
    """Get and/or build index files for references.

    :param str reference_dir: directory containing curated APPRIS reference files
    :param str | None ensembl_id: Ensembl gene or transcript ID, with version number
    :param str | None ref: path to reference FASTA used in alignment. Must be bowtie2 FM-index and samtools faidx indexed.
    :param str transcript_gff: GFF/GTF file containing transcript metafeatures and exon features, in 5' to 3' order, \
    regardless of strand. Ordering is essential.
    :param str gff_reference: reference FASTA corresponding to the GFF features
    :param str outdir: optional output dir for the reference files
    :return tuple: (ref_fa, gff, gff_ref) filepaths
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    ref_fa = ref
    gff = transcript_gff
    gff_ref = gff_reference

    # Determine if the provided Ensembl ID is found in the curated APPRIS references
    if ensembl_id is not None:
        if ref is not None:
            raise RuntimeError("Both an Ensembl ID and a reference FASTA were provided. Please choose one.")

        ref_fa, gff = get_ensembl_references(reference_dir=reference_dir, ensembl_id=ensembl_id, outdir=outdir)
        gff_ref = os.path.join(reference_dir, GRCH38_FASTA)

    return ref_fa, gff, gff_ref


def sim_workflow(am, vcf, ensembl_id=ri.ReadEditor.DEFAULT_ENSEMBL_ID, reference_dir=ri.ReadEditor.DEFAULT_REFERENCE_DIR,
                 ref=ri.ReadEditor.DEFAULT_REF, primers=ri.ReadEditor.DEFAULT_PRIMERS, outdir=ri.ReadEditor.DEFAULT_OUTDIR,
                 output_prefix=ri.ReadEditor.DEFAULT_PREFIX, single_end=ri.ReadEditor.DEFAULT_SINGLE_END,
                 filter_edited=ri.ReadEditor.DEFAULT_FILTER, random_seed=ri.ReadEditor.DEFAULT_SEED):
    """Runs the satmut_utils sim workflow.

    :param str am: SAM/BAM file to edit into
    :param str vcf: VCF file specifying variants to edit
    :param str | None ensembl_id: Ensembl gene or transcript ID, with version number
    :param str reference_dir: directory containing curated APPRIS reference files. Default ./references.
    :param str ref: indexed reference FASTA
    :param str | None primers: feature file of primer locations for read masking and primer detection
    :param str outdir: Optional output directory to store generated FASTQs and BAM
    :param str | None output_prefix: Optional output prefix for the FASTQ(s) and BAM; if None, use same prefix as BAM.
    :param bool single_end: does the input SAM/BAM contain single end reads? Default False, paired end reads.
    :param bool filter_edited: filter the edited BAM for those reads/read pairs that were edited? Default True.
    :param int random_seed: seed for random qname sampling
    :return tuple: (str | None, str, str | None) paths of the edited BAM, R1 FASTQ, R2 FASTQ
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    ref_fa = get_sim_reference(
        reference_dir=reference_dir, ensembl_id=ensembl_id, ref=ref, outdir=os.path.join(outdir, "references"))

    out_prefix = output_prefix
    if output_prefix is None:
        out_prefix = fu.remove_extension(os.path.basename(am))

    # Run the editing workflow
    output_bam, zipped_r1_fastq, zipped_r2_fastq = ri.ReadEditor.workflow(
        am=am, vcf=vcf, ref=ref_fa, primers=primers, output_dir=outdir, output_prefix=out_prefix,
        single_end=single_end, filter_edited=filter_edited, random_seed=random_seed)

    return output_bam, zipped_r1_fastq, zipped_r2_fastq


def call_workflow(f1, f2, r1_fiveprime_adapters, r1_threeprime_adapters,
                  ensembl_id=VariantCaller.VARIANT_CALL_ENSEMBL_ID,
                  reference_dir=VariantCaller.VARIANT_CALL_REFERENCE_DIR,
                  ref=VariantCaller.VARIANT_CALL_REF, transcript_gff=VariantCaller.VARIANT_CALL_GFF,
                  gff_reference=VariantCaller.VARIANT_CALL_GFF_REF, targets=VariantCaller.VARIANT_CALL_TARGET,
                  outdir=VariantCaller.VARIANT_CALL_OUTDIR, primers=VariantCaller.VARIANT_CALL_PRIMERS,
                  primer_fa=UMIExtractor.PRIMER_FASTA, primer_nm_allowance=UMIExtractor.PRIMER_NM_ALLOW,
                  consensus_dedup=VariantCaller.VARIANT_CALL_CDEDUP, umi_regex=AMP_UMI_REGEX,
                  contig_del_thresh=ConsensusDeduplicator.CONTIG_DEL_THRESH,
                  min_bq=VariantCaller.VARIANT_CALL_MIN_BQ, max_nm=VariantCaller.VARIANT_CALL_MAX_NM,
                  min_supporting_qnames=VariantCaller.VARIANT_CALL_MIN_DP,
                  max_mnp_window=VariantCaller.VARIANT_CALL_MAX_MNP_WINDOW,
                  include_n=not VariantCaller.VARIANT_CALL_EXCLUDE_N,
                  nthreads=FastqPreprocessor.NCORES, ntrimmed=FastqPreprocessor.NTRIMMED,
                  trim_bq=FastqPreprocessor.TRIM_QUALITY, omit_trim=FastqPreprocessor.TRIM_FLAG,
                  mut_sig=VariantCaller.VARIANT_CALL_MUT_SIG):
    r"""Runs the satmut_utils call workflow.

    :param str f1: path of the R1 FASTQ
    :param str f2: path of the R2 FASTQ
    :param str r1_fiveprime_adapters: comma-delimited 5' adapters to trim from R1
    :param str r1_threeprime_adapters: comma-delimited 3' adapters to trim from R1
    :param str | None ensembl_id: Ensembl gene or transcript ID, with version number
    :param str reference_dir: directory containing curated APPRIS reference files. Default ./references.
    :param str | None ref: path to reference FASTA used in alignment. Must be bowtie2 FM-index and samtools faidx indexed.
    :param str | None transcript_gff: GFF/GTF file containing transcript metafeatures and exon features, in 5' to 3' \
    order, regardless of strand. Ordering is essential.
    :param str | None gff_reference: reference FASTA corresponding to the GFF features
    :param str | None targets: BED or GFF containing target regions to call variants in
    :param str outdir: optional output dir for the results
    :param str | None primers: BED or GFF file containing primers to mask. Must contain a strand field.
    :param str | None primer_fa: primer FASTA. Useful for annotating reads with originating tile. Default None.
    :param int primer_nm_allowance: Max edit distance a R2 can have to match a primer. Default 3.
    :param bool consensus_dedup: should consensus bases be generated during deduplication? Default False.
    :param str umi_regex: regex for matching the UMIs (see umi_tools for docs)
    :param int contig_del_thresh: max deletion length for which del/N gaps in the merged R2 contig are called
    :param int min_bq: min base qual; should be >= 1 such that masked primers do not contribute to depth
    :param int max_nm: max edit distance to consider a read for variant calling
    :param int min_supporting_qnames: min number of fragments with R1-R2 concordant coverage for which to keep a variant
    :param int max_mnp_window: max number of consecutive nucleotides to search for MNPs
    :param bool include_n: include variant calls to an "N"
    :param int nthreads: Number of threads to use for SAM/BAM operations. Default 0 (autodetect).
    :param int ntrimmed: Max number of adapters to trim from each read
    :param int trim_bq: quality score for cutadapt quality trimming at the 3' end. Default 15.
    :param bool omit_trim: flag to turn off adapter and 3' base quality trimming. Default False.
    :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNK.
    :return tuple: (VCF, BED) filepaths
    """

    # Get and index the references
    ref_fa, gff, gff_ref = get_call_references(
        reference_dir=reference_dir, ensembl_id=ensembl_id, ref=ref, transcript_gff=transcript_gff,
        gff_reference=gff_reference, outdir=os.path.join(outdir, "references"))

    if mut_sig not in VariantCaller.VARIANT_CALL_VALID_MUT_SIGS:
        raise RuntimeError("Mutation signature %s not one of NNN, NNK, or NNS" % mut_sig)

    fqp_r1 = f1
    fqp_r2 = f2
    if consensus_dedup:
        # Run the first part of umi_tools based deduplication: UMI extraction. Note this should take place
        # before we have trimmed adapters as for umi_tools the adapter is used for "anchoring" the UMI.
        ue = UMIExtractor(r1_fastq=f1, r2_fastq=f2, umi_regex=umi_regex,
                          primer_fasta=primer_fa, primer_nm_allow=primer_nm_allowance, outdir=outdir)
        fqp_r1 = ue.r1_out_fastq
        fqp_r2 = ue.r2_out_fastq

    # Run the FASTQ preprocessing workflow which includes adapter trimming and 3' BQ trimming
    fqp = FastqPreprocessor(
        f1=fqp_r1, f2=fqp_r2, r1_fiveprime_adapters=r1_fiveprime_adapters, r1_threeprime_adapters=r1_threeprime_adapters,
        outdir=outdir, ncores=nthreads, trim_bq=trim_bq, ntrimmed=ntrimmed, no_trim=omit_trim, validate=True)

    # Run local alignment; handle the ncores/nthreads option for cutadapt versus bowtie2 options
    bowtie2_nthreads = 1 if nthreads == 0 else nthreads
    bta = baw(f1=fqp.trimmed_f1, ref=ref, f2=fqp.trimmed_f2, outdir=outdir, outbam=None,
              local=True, nthreads=bowtie2_nthreads)

    # Run deduplication which may be done in at least two ways (standard and consensus)
    preproc_in_bam = bta.output_bam
    if consensus_dedup:
        # Run consensus deduplication (majority vote for each base call within a read's UMI network/group)
        rg = ReadGrouper(in_bam=bta.output_bam, outdir=outdir)
        cdp = ConsensusDeduplicatorPreprocessor(group_bam=rg.group_bam, outdir=outdir, nthreads=nthreads)
        cd = ConsensusDeduplicator(in_bam=cdp.preprocess_bam, ref=ref_fa, outdir=outdir, out_bam=None,
                                   nthreads=nthreads, contig_del_thresh=contig_del_thresh)
        preproc_in_bam = cd.out_bam

    # Optionally run primer masking
    vc_in_bam = preproc_in_bam
    if primers is not None:
        rm = ReadMasker(
            in_bam=preproc_in_bam, feature_file=primers, consensus_dedup=consensus_dedup, outdir=outdir)
        vc_in_bam = rm.out_bam

    # Initialize the VariantCaller and prepare the alignments
    vc = VariantCaller(
        am=vc_in_bam, targets=targets, ref=ref_fa, trx_gff=gff, gff_ref=gff_ref, primers=primers,
        output_dir=outdir,  nthreads=nthreads, exclude_n=not include_n, mut_sig=mut_sig)

    # Run variant calling
    output_vcf, output_bed = vc.workflow(min_bq, max_nm, min_supporting_qnames, max_mnp_window)

    return output_vcf, output_bed


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    # Call the assigned workflow
    parsed_args.func(parsed_args)


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
