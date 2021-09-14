#!/usr/bin/env/python
"""Runs the workflow for consensus deduplication."""

import argparse
import logging
import sys

from analysis.read_preprocessor import ConsensusDeduplicator, UMITOOLS_UG_TAG
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

    parser.add_argument("-b", "--in_bam", type=str, required=True,
                        help='Input BAM file. Should be grouped with updated tags and sorted by tag name.')

    parser.add_argument("-r", "--reference", type=str, required=True,
                        help='FM-indexed reference FASTA for realignment of consensus reads.')

    parser.add_argument("-gt", "--group_tag", type=str, default=UMITOOLS_UG_TAG,
                        help='BAM tag with the UMI group/network ID. Default %s.' % UMITOOLS_UG_TAG)

    parser.add_argument("-d", "--outdir", type=str, default=".", help='Output directory.')

    parser.add_argument("-o", "--outbam", type=none_or_str, default="None", help='Optional output BAM file path.')

    parser.add_argument("-t", "--nthreads", type=int, default=ConsensusDeduplicator.DEFAULT_NTHREADS,
                        help='Number of threads to use for bowtie2 alignment and BAM operations. Default %i.'
                             % ConsensusDeduplicator.DEFAULT_NTHREADS)

    parser.add_argument("-dt", "--contig_del_thresh", type=int, default=ConsensusDeduplicator.CONTIG_DEL_THRESH,
                        help='Max length of deletions to consider as true deletions. ')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def workflow(in_bam, ref, group_tag=UMITOOLS_UG_TAG, outdir=".", out_bam=None,
             nthreads=ConsensusDeduplicator.DEFAULT_NTHREADS,
             contig_del_thresh=ConsensusDeduplicator.CONTIG_DEL_THRESH):
    r"""Runs the consensus deduplication workflow.

    :param str in_bam: input alignments with UMI network/group ID in alignment tag
    :param str ref: reference FASTA to realign to
    :param str group_tag: BAM tag for the group ID. Default UG.
    :param str outdir: Output directory
    :param str | None out_bam: Optional filepath of the output BAM
    :param int nthreads: number of threads to use for alignment
    :param int contig_del_thresh: max del length for which gaps in a merged R2 contig are called as del instead of N. \
    This option is needed in cases where R2s originating from the same input fragment do not overlap.
    """

    ConsensusDeduplicator(in_bam=in_bam, ref=ref, group_tag=group_tag, outdir=outdir, out_bam=out_bam,
                          nthreads=nthreads, contig_del_thresh=contig_del_thresh)


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    workflow(in_bam=parsed_args["in_bam"], ref=parsed_args["reference"], group_tag=parsed_args["group_tag"],
             outdir=parsed_args["outdir"], out_bam=parsed_args["outbam"], nthreads=parsed_args["nthreads"],
             contig_del_thresh=parsed_args["contig_del_thresh"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
