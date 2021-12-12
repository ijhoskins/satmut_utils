#!/usr/bin/env python3
"""Runs basic in silico read generation."""

import argparse
import logging
import sys

from analysis.aligners import *
from analysis.seq_utils import BAM_SUFFIX
from definitions import LOG_FORMATTER
from prototype.read_generators import *
from scripts.run_bowtie2_aligner import workflow as ba_workflow

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
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
    parser.add_argument("-t", "--feature_file", type=str, help='Target feature file, e.g. BED, GFF, or GTF')

    parser.add_argument("-d", "--output_dir", type=str, default="read_gen_results",
                        help='Output directory for FASTQs and BAM')

    parser.add_argument("-o", "--output_prefix", type=str,
                        default=".".join([datetime.datetime.now().strftime("%d%b%Y.%I:%M%p"), "read_gen"]),
                        help='Output prefix for FASTQs and BAM')

    parser.add_argument("-n", "--nreads", type=int, default=DEFAULT_NREADS,
                        help='Number of fragments to generate. Default %i. Feature file scores override defaults.'
                             % DEFAULT_NREADS)

    parser.add_argument("-x", "--reference", type=str, help='Corresponding reference FASTA.')

    parser.add_argument("-p", "--paired", action="store_true", help='Generate PE reads instead of SE reads')

    parser.add_argument("-r", "--rna", action="store_true", help='Generate RNA reads instead of DNA reads.')

    parser.add_argument("-l", "--read_length", type=int, default=DEFAULT_READ_LEN,
                        help='Read length to generate. Default %i.' % DEFAULT_READ_LEN)

    parser.add_argument("-f", "--frag_length", type=int, default=DEFAULT_FRAG_LEN,
                        help='Mean fragment length to simulate. Default %i.' % DEFAULT_FRAG_LEN)

    parser.add_argument("-s", "--slop_length", type=int, default=DEFAULT_BP_SLOP,
                        help='Slop length around exon features for DNA read generation. Default %i.' % DEFAULT_BP_SLOP)

    parser.add_argument("-m", "--make_amplicons", action="store_true",
                        help='Flag to generate reads from the termini of the targets. Consider setting slop_length to 0.')

    parser.add_argument("-a", "--add_snps", action="store_true", help='Add SNP errors to reads.')

    parser.add_argument("-i", "--add_indels", action="store_true", help='Add InDel errors to reads.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(feature_file, output_dir, output_prefix, ref, nreads=DEFAULT_NREADS, paired=True, rna=False,
             read_length=DEFAULT_READ_LEN, frag_length=DEFAULT_FRAG_LEN, slop_length=DEFAULT_BP_SLOP,
             make_amplicons=False, add_snps=False, add_indels=False):
    r"""Runs the read generation workflow.

    :param str feature_file: Feature file containing features to make reads for. Target BED for DNA, transcript GTF for RNA
    :param str output_dir: Optional output directory to store generated FASTQs and BAM
    :param str output_prefix: Optional output prefix for the FASTQ(s) and BAM; if None, use same prefix as feature file
    :param str ref: reference FASTA corresponding to the features
    :param int nreads: number of reads to generate. This overrides scores present in the feature file.
    :param bool paired: should paired-end reads be created? Default True.
    :param bool rna: should RNA reads be created? Assumes features are exons and metafeatures are transcripts. \
    Default False, make DNA reads.
    :param int read_length: length of reads to generate
    :param int frag_length: average fragment length to simulate
    :param int slop_length: number of bases to slop around targets for DNA read generation
    :param bool make_amplicons: should reads be generated from the ends of the target regions? Default False. \
    Make uniform random ends.
    :param bool add_snps: add SNP errors to reads
    :param bool add_indels: add InDel errors to reads
    :return tuple: (r1, r2, bam) or (r1, None, bam) corresponding to FASTQ filepaths and indexed BAM
    """

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    read_gen = ReadGenerator(feature_file, ref, paired, rna, read_length, frag_length, slop_length, make_amplicons)
    fastqs = read_gen.fastqs_from_features(output_dir, output_prefix, nreads, add_snps, add_indels)

    # Now for bowtie2 alignment
    fq1, fq2 = fastqs
    if fq2 is None:
        outbam = replace_extension(fq1, BAM_SUFFIX)
    else:
        outbam = ".".join([os.path.basename(os.path.commonprefix([fq1, fq2])), BAM_SUFFIX])

    ba_workflow(ref=ref, local=True, f1=fq1, f2=fq2, outdir=output_dir, outbam=outbam)

    return fastqs + (outbam,)


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

    workflow(feature_file=parsed_args["feature_file"], output_dir=parsed_args["output_dir"],
             output_prefix=parsed_args["output_prefix"], nreads=parsed_args["nreads"], ref=parsed_args["reference"],
             paired=parsed_args["paired"], rna=parsed_args["rna"], read_length=parsed_args["read_length"],
             frag_length=parsed_args["frag_length"], slop_length=parsed_args["slop_length"],
             make_amplicons=parsed_args["make_amplicons"], add_snps=parsed_args["add_snps"],
             add_indels=parsed_args["add_indels"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
