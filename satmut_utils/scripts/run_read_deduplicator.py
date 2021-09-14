#!/usr/bin/env/python
"""Runs read deduplication."""

import argparse
import logging
import os
import sys

from analysis.read_preprocessor import FastqPreprocessor, ReadDeduplicator
from analysis.seq_utils import TRANSCRIPTOME_FASTA
import core_utils.file_utils as fu
from core_utils.string_utils import none_or_str
from scripts.run_fastq_preprocessor import workflow as fpw
from scripts.run_bowtie2_aligner import workflow as baw

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
    parser.add_argument("-f1", "--fast1", required=True, type=str, help='R1 FASTQ.')

    parser.add_argument("-f2", "--fast2", required=True, type=str, help='R2 FASTQ.')

    parser.add_argument("-r", "--reference", type=str, default=TRANSCRIPTOME_FASTA,
                        help='Corresponding reference FASTA. Default %s.' % TRANSCRIPTOME_FASTA)

    parser.add_argument("-pf", "--primer_fa", type=none_or_str,
                        help='If deduplicate=True, this may optionally be set to append originating primers to qnames.')

    parser.add_argument("-o", "--outdir", type=str, default=".", help='Output directory.')

    parser.add_argument("-c", "--ncores", type=int, default=FastqPreprocessor.NCORES,
                        help='Number of cores/threads to use for cutadapt, bowtie2 alignment, SAM/BAM operations. '
                             'Default %i.' % FastqPreprocessor.NCORES)

    parser.add_argument("-n", "--ntrimmed", type=int, default=FastqPreprocessor.NTRIMMED,
                        help='Max number of adapters to trim from each read. Default %i.' % FastqPreprocessor.NTRIMMED)

    parser.add_argument("-q", "--min_trim_bq", type=int, default=FastqPreprocessor.TRIM_QUALITY,
                        help='Min BQ for 3\' trimming. Default %i.' % FastqPreprocessor.TRIM_QUALITY)

    parser.add_argument("-a", "--amp", action="store_true",
                        help='Is the data generated with AMP technology (RACE-like library with 8 bp UMI on R1)?')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(f1, f2, ref=TRANSCRIPTOME_FASTA, primer_fa=None, outdir=None,
             ncores=FastqPreprocessor.NCORES, ntrimmed=FastqPreprocessor.NTRIMMED,
             min_trim_bq=FastqPreprocessor.TRIM_QUALITY, is_tileseq=FastqPreprocessor.IS_TILESEQ):
    r"""Runs basic umi_tools read deduplication.

    :param str f1: path of the R1 FASTA/FASTQ
    :param str f2: path of the R2 FASTA/FASTQ
    :param str ref: reference FASTA corresponding to the features
    :param str primer_fa: primer sequences FASTA; names must be unique.
    :param str outdir: Optional output dir for the results
    :param int ncores: Number of CPU cores to use for cutadapt. Default 1.
    :param int ntrimmed: Max number of adapters to trim from each read
    :param int min_trim_bq: quality score for quality trimming at the 3' end. Default 15.
    :param bool is_tileseq: flag indicating the data is Tileseq-generated; otherwise, assume ArcherDX UMI adapters. \
    Default True.
    """

    # Run the first part of umi_tools based deduplication: UMI extraction
    prepend = True if primer_fa is not None else False

    rd = ReadDeduplicator(r1_fastq=f1, r2_fastq=f2, ref=ref, prepend_primer=prepend, gsp2_fasta=primer_fa,
                          outdir=outdir, nthreads=ncores, single_strands=is_tileseq)

    # Run the FASTQ preprocessing workflow
    fqp = fpw(f1=rd.r1_out_fastq, f2=rd.r2_out_fastq, outdir=outdir,
              ncores=ncores, min_quality=min_trim_bq, ntrimmed=ntrimmed)

    # Run the local alignment workflow
    # Need to change dir to scratch in case temp BAMs are created, as bowtie2 does not have an option for setting
    # the temp directory
    call_dir = os.getcwd()
    scratch_dir = os.getenv("SCRATCH", "/tmp")
    os.chdir(scratch_dir)
    bta = baw(f1=fqp.trimmed_f1, ref=ref, f2=fqp.trimmed_f2, outdir=outdir,
              outbam=None, local=True, nthreads=ncores)
    os.chdir(call_dir)

    # Run the second half of deduplication
    dedup_bam = fu.replace_extension(bta.output_bam, rd.DEDUP_BAM_SUFFIX)
    rd.umitools_dedup(in_bam=bta.output_bam, out_bam=dedup_bam)

    return dedup_bam


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(f1=parsed_args["fast1"], f2=parsed_args["fast2"], ref=parsed_args["reference"],
             primer_fa=parsed_args["primer_fa"], outdir=parsed_args["outdir"],
             ncores=parsed_args["ncores"], ntrimmed=parsed_args["ntrimmed"],
             min_trim_bq=parsed_args["min_trim_bq"], is_tileseq=not parsed_args["amp"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
