#!/usr/bin/env python3
"""Postprocesses satmut_utils *summary.txt files to express MAVE-HGVS annotations relative to target regions."""

import argparse
import glob
import logging
import os
import pybedtools
import re
import string
import sys

from analysis.coordinate_mapper import AminoAcidMapper
from core_utils.file_utils import FILE_DELIM, FILE_NEWLINE, replace_extension
from satmut_utils.definitions import LOG_FORMATTER

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"


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
    parser.add_argument("-i", "--input_dir", type=str, required=True,
                        help='Input directory containing satmut_utils *summary.txt files.')

    parser.add_argument("-g", "--transcript_gff", required=True, type=str,
                        help='GFF file containing transcript metafeatures and exon, CDS, and stop_codon features, '
                             'from 5\' to 3\', regardless of strand.')

    parser.add_argument("-k", "--gff_reference", type=str, required=True, help='Reference FASTA for the GFF.')

    parser.add_argument("-t", "--targets", required=True,
                        help='Target BED file. Contig names in the target file should match the contig name '
                             'in the reference FASTA.')

    parser.add_argument("-o", "--output_dir", type=str, default=".",
                        help='Optional output directory. Default current working directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def correct_mave_hgvs_nt(mut_info_tuple, target_start):
    """Corrects positions of a MAVE-HGVS nt annotation to be relative to a target.

    :param collections.namedtuple mut_info_tuple: mutation info relative to whole reference
    :param int target_start: 0-based coordinate of the target start
    :return analysis.coordinate_mapper.MUT_INFO_TUPLE: corrected hgvs_nt annotation
    """

    mave_hgvs_nt_items = mut_info_tuple.mave_hgvs_nt.split(",")
    mave_hgvs_nt_items_corrected = []

    for annot in mave_hgvs_nt_items:
        mave_hgvs_nt_split = annot.split(".")
        mave_hgvs_nt_pos = re.sub("[->A-Za-z*]", "", mave_hgvs_nt_split[1])
        mave_hgvs_nt_char = re.sub("[_0-9]", "", mave_hgvs_nt_split[1])

        if re.search("_", mave_hgvs_nt_pos):
            mave_hgvs_nt_pos_split = mave_hgvs_nt_pos.split("_")
            mave_hgvs_nt_pos_adj = "_".join(list(map(str, [int(e) - target_start for e in mave_hgvs_nt_pos_split])))
        else:
            mave_hgvs_nt_pos_adj = str(int(mave_hgvs_nt_pos) - target_start)

        # For mave_hgvs_nt, integer and character tokens can be simply added back together
        mave_hgvs_nt_correct = mave_hgvs_nt_split[0] + "." + mave_hgvs_nt_pos_adj + mave_hgvs_nt_char

        mave_hgvs_nt_items_corrected.append(mave_hgvs_nt_correct)

    if len(mave_hgvs_nt_items_corrected) > 1:
        # Strip the "c." prefix
        mave_hgvs_nt_items_corrected_strip = [e.split(".")[1] for e in mave_hgvs_nt_items_corrected]
        res = "c.[" + ";".join(mave_hgvs_nt_items_corrected_strip) + "]"
    else:
        res = mave_hgvs_nt_items_corrected[0].split(":")[1]

    return res


def correct_mave_hgvs_pro(mut_info_tuple, target_aa_start):
    """Corrects positions of a MAVE-HGVS pro annotation to be relative to a target.

    :param collections.namedtuple mut_info_tuple: mutation info relative to whole reference
    :param int target_aa_start: 0-based amino acid of the target start
    :return analysis.coordinate_mapper.MUT_INFO_TUPLE: corrected hgvs_pro annotation
    """

    ascii_letters_punct = set(string.ascii_letters + string.punctuation)
    mave_hgvs_pro_items = mut_info_tuple.mave_hgvs_pro.split(",")
    mave_hgvs_pro_items_corrected = []

    for annot in mave_hgvs_pro_items:

        mave_hgvs_pro_correct = []
        last_is_digit = False
        curr_pos = []

        for e in annot:
            if e in ascii_letters_punct:
                if last_is_digit:
                    aa_pos = int("".join(curr_pos))
                    pos_correct = aa_pos - target_aa_start
                    mave_hgvs_pro_correct.extend([str(pos_correct), e])
                    last_is_digit = False
                else:
                    mave_hgvs_pro_correct.append(e)
            else:
                # Here we have a digit
                curr_pos.append(e)
                last_is_digit = True

        mave_hgvs_pro_items_corrected.append("".join(mave_hgvs_pro_correct))

    if len(mave_hgvs_pro_items_corrected) > 1:
        # Strip the "p." prefix
        mave_hgvs_pro_items_corrected_strip = [e[2:] for e in mave_hgvs_pro_items_corrected]
        res = "p.[" + ";".join(mave_hgvs_pro_items_corrected_strip) + "]"
    else:
        res = mave_hgvs_pro_items_corrected[0]

    return res


def add_annotation(summary_file, aam, inframe_interval_starts, target_positions, outdir="."):
    """Adds MAVE-HGVS annotations and writes to a new file.

    :param str summary_file: satmut_utils *summary.txt file
    :param analysis.coordinate_mapper.AminoAcidMapper aam: AminoAcidMapper object
    :param set inframe_interval_starts: 0-based interval starts in-frame with the CDS start
    :param set target_positions: all targeted positions
    :param str outdir: Optional output directory. Default current directory
    """

    outpath = os.path.join(outdir, replace_extension(os.path.basename(summary_file), "mave_hgvs.txt"))

    with open(summary_file, "r") as infile, open(outpath, "w") as outfile:

        for i, line in enumerate(infile):

            # Write the header
            if i == 0:
                linesplit = line.rstrip(FILE_NEWLINE).split(FILE_DELIM)
                new_header = linesplit + ["MAVE_HGVS_NT", "MAVE_HGVS_PRO"]
                outfile.write(FILE_DELIM.join(new_header) + FILE_NEWLINE)
                continue

            linesplit = line.rstrip(FILE_NEWLINE).split(FILE_DELIM)
            trx_id = linesplit[0].split("|")[0]
            pos = int(linesplit[1])

            # Filter variants that are not in a target position as all variants must have a valid MAVE-HGVS annotation
            if pos not in target_positions:
                continue

            _, cds_start_offset, _, _, _ = aam.cds_info[trx_id]

            # Get the target coordinates to adjust the mave_nt and mave_pro annotations
            target_start = min(inframe_interval_starts)
            target_cds_start = target_start - cds_start_offset

            mut_info_tuple = aam.get_codon_and_aa_changes(
                trx_id=trx_id, pos=int(linesplit[1]), ref=linesplit[3], alt=linesplit[4])

            # This is 0-based
            target_aa_start = round(target_cds_start / 3)

            mave_hgvs_nt_correct = correct_mave_hgvs_nt(mut_info_tuple, target_cds_start)
            mave_hgvs_pro_correct = correct_mave_hgvs_pro(mut_info_tuple, target_aa_start)

            annot_line = FILE_DELIM.join(linesplit + [mave_hgvs_nt_correct, mave_hgvs_pro_correct]) + FILE_NEWLINE

            outfile.write(annot_line)


def workflow(indir, transcript_gff, gff_ref, targets, outdir="."):
    r"""Runs the MAVE-HGVS annotation workflow.

    :param str indir: Input directory containing summary.txt files
    :param str transcript_gff: GFF/GTF file containing transcript metafeatures and exon features, in 5' to 3' order, \
    regardless of strand
    :param str gff_ref: reference FASTA corresponding to the GFF features
    :param str targets: target BED for adjusting MAVE-HGVS nt and protein annotations
    :param str outdir: Optional output directory. Default current directory
    """

    summary_files = glob.glob(os.path.join(indir, "*summary.txt"))
    aam = AminoAcidMapper(gff=transcript_gff, ref=gff_ref, outdir=outdir)
    bedtool = pybedtools.BedTool(targets)

    target_positions = set()
    inframe_interval_starts = set()
    for interval in bedtool:
        target_positions |= set(range(interval.start, interval.stop))

        # If the target does not start flush with a codon, expand start coordinate to include the full start codon
        _, cds_start_offset, _, _, _ = aam.cds_info[str(interval.chrom).split("|")[0]]

        interval_start = interval.start
        while ((interval_start - cds_start_offset) % 3) != 0:
            interval_start -= 1

        inframe_interval_starts.add(interval_start)

    for file in summary_files:
        add_annotation(file, aam, inframe_interval_starts, target_positions, outdir)


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    outdir = parsed_args["output_dir"]
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    logger.info("Started %s" % sys.argv[0])

    workflow(indir=parsed_args["input_dir"], transcript_gff=parsed_args["transcript_gff"],
             gff_ref=parsed_args["gff_reference"], targets=parsed_args["targets"], outdir=parsed_args["output_dir"])

    logger.info("Completed %s" % sys.argv[0])


if __name__ == "__main__":
    main()
