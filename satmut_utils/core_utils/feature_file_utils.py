#!/usr/bin/env/python
"""Collection of feature file (BED, GFF, etc.) manipulation utilities."""

import collections
import os
import pybedtools
import subprocess
import tempfile

from analysis.aligners import BowtieConfig
import analysis.seq_utils as su
from core_utils.file_utils import FILE_DELIM, FILE_NEWLINE, replace_extension

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")

BED_FILETYPE = "bed"
GFF_FILETYPE = "gff"
GTF_FILETYPE = "gtf"
GFF_FEATURE_TYPE_FIELD = 2
GFF_ATTR_FIELD_DELIM = "; "
GFF_ATTR_SUBFIELD_DELIM1 = "="
GFF_ATTR_SUBFIELD_DELIM2 = " "
GFF_ATTR_SUBFIELD_VAL_STR = '"'
GFF_ATTR_GENE_ID = "gene_id"
GFF_ATTR_TRANSCRIPT_ID = "transcript_id"
GFF_ATTR_EXON_ID = "exon_number"
GFF_ATTR_ID_ID = "ID"

PYBEDTOOLS_NULL_CHARS = {".", ""}
DEFAULT_BP_SLOP = 25

GFF_ATTR_FIELD = 8
GFF_NFIELDS = 9

# Indices from the start for SAM/BAM fields intersected in -wb mode
BED_INTERSECT_WB_READ_CONTIG_INDEX = 0
BED_INTERSECT_WB_READ_START_INDEX = 1
BED_INTERSECT_WB_READ_STOP_INDEX = 2
BED_INTERSECT_WB_READ_QNAME_INDEX = 3
BED_INTERSECT_WB_READ_STRAND_INDEX = 5
BED_INTERSECT_WB_READ_NFIELDS = 12

# Offsets from the B contig field; field length to be determined dynamically, as feature files may have various nfields
# Note the contig field index is just the -(nfields of B); the below offsets should be added to the nfields of B
# to get the respective B field, depending on the type of B feature file
BED_INTERSECT_WB_B_BED_START_OFFSET = 1
BED_INTERSECT_WB_B_BED_STOP_OFFSET = 2
BED_INTERSECT_WB_B_BED_NAME_OFFSET = 3
BED_INTERSECT_WB_B_BED_SCORE_OFFSET = 4
BED_INTERSECT_WB_B_BED_STRAND_OFFSET = 5

BED_INTERSECT_WB_B_GFF_TYPE_OFFSET = 2
BED_INTERSECT_WB_B_GFF_START_OFFSET = 3
BED_INTERSECT_WB_B_GFF_STOP_OFFSET = 4
BED_INTERSECT_WB_B_GFF_SCORE_OFFSET = 5
BED_INTERSECT_WB_B_GFF_STRAND_OFFSET = 6
BED_INTERSECT_WB_B_GFF_ATTR_OFFSET = 8

BED_GROUPBY_DELIM = ","

COORD_TUPLE = collections.namedtuple("coord_tuple", "contig, start, stop, name, strand, score, allowable_coords")


def intersect_features(ff1, ff2, outfile=None, as_bedtool=False, **kwargs):
    """bedtools intersect of two feature files. Feature file 2 is the smaller and held in memory.

    :param str ff1: first feature file, may be either VCF, BED, GFF, or BAM
    :param str ff2: second feature file, usually VCF, BED, GFF
    :param str | None outfile: optional output file
    :param bool as_bedtool: return results as a bedtool object?
    :return pybedtools.BedTool | str: if as_bedtool=True, a BedTool; otherwise the intersected result filename
    """

    out = outfile
    if outfile is None and not as_bedtool:
        out = tempfile.NamedTemporaryFile(mode="w", suffix=".intersect", delete=False).name

    ff1_bt = pybedtools.BedTool(ff1)
    ff2_bt = pybedtools.BedTool(ff2)

    # Always print the header as this can cause issues with opening the file downstream
    # (e.g. pysam.VariantFile on a header-less VCF)
    intersect_bt = ff1_bt.intersect(ff2_bt, header=True, **kwargs)

    if as_bedtool:
        return intersect_bt

    intersect_bt.moveto(out)
    return out


def sort_feature_file(feature_file, output=None, *args, **kwargs):
    """Lexicographically sorts a BED.

    :param str feature_file: unsorted input feature file
    :param str | None output: optional name of the output file
    :param args: additional args to pass to the bedtools sort command
    :param kwargs: key-value arguments to pass to the bedtools sort command
    :return str: name of the output file
    """

    bedtool = pybedtools.BedTool(feature_file)
    sorted_bedtool = bedtool.sort(*args, **kwargs)

    if output is not None:
        sorted_bedtool.moveto(output)
        return output
    else:
        return sorted_bedtool.fn


def get_genome_file(ref, output_file=None):
    """Returns a genome file for the reference (a tab-delimited file containing contig lengths).

    :param str ref: bowtie2 indexed reference FASTA
    :param str | None output_file: Optional name of output file. If None a file will be written to same dir as ref.
    :return str: name of the genome file
    :raises RuntimeError: if no FM index files are found for the reference FASTA
    """

    # Tests for index files by side effect
    _ = BowtieConfig(ref)

    inspect_cmd = ["bowtie2-inspect", "-s", ref]
    grep_cmd = ["grep", "Sequence"]
    cut_cmd = ["cut", "-f2-"]

    outfile = output_file
    if output_file is None:
        outfile = os.path.join(os.path.dirname(ref), replace_extension(os.path.basename(ref), "genome_file.txt"))

    if os.path.exists(outfile):
        return outfile

    with open(outfile, "w") as out_fh:

        inspect_p = subprocess.Popen(inspect_cmd, stdout=subprocess.PIPE)
        grep_p = subprocess.Popen(grep_cmd, stdin=inspect_p.stdout, stdout=subprocess.PIPE)
        cut_p = subprocess.Popen(cut_cmd, stdin=grep_p.stdout, stdout=out_fh)

        _ = cut_p.wait()
        inspect_p.stdout.close()
        grep_p.stdout.close()

    return outfile


class DuplicateFeatureException(Exception):
    """Exception for the existence of duplicate features (those with exact same coordinates)."""
    pass


def store_coords(feature_file, feature_slop=0, primer_allowable=False, use_name=True):
    r"""Stores primer coordinates from a feature file in a consistent dictionary format.

    :param str feature_file: BED, GFF, or GTF
    :param int feature_slop: Bases to slop the start or stop coordinates for allowable_starts attr of coord_tuple
    :param bool primer_allowable: Make the allowable_starts attribute contain up to the 3' end of the primer. \
    Default False. Keep at False for primer enumeration, and set to True to determine primer region bases to mask.
    :param bool use_name: Should the feature name field be used as the key? Default True. Otherwise, use the contig \
    and coordinate range as the key.
    :return collections.OrderedDict: {"contig:start-stop" | name: COORD_TUPLE}
    :raises RuntimeError: if the strand field is absent or contains an invalid character
    """

    feature_dict = collections.OrderedDict()
    observed_features = set()

    ff_bedtool = pybedtools.BedTool(feature_file)
    for feature in ff_bedtool:

        # Always generate a name even if it is redundant
        feature_coords = su.COORD_FORMAT_STRAND.format(str(feature.chrom), feature.start, feature.stop, feature.strand)

        if feature_coords in observed_features:
            raise DuplicateFeatureException("Please provide a feature file with duplicate features removed.")

        # Get some other basic attrs if they exist
        feature_name = str(feature.name) if feature.name not in PYBEDTOOLS_NULL_CHARS else None
        feature_strand = su.Strand(str(feature.strand)) if feature.strand not in PYBEDTOOLS_NULL_CHARS else None
        feature_score = float(feature.score) if feature.score not in PYBEDTOOLS_NULL_CHARS else 1.0
        feature_len = len(feature)

        # Compute allowable start coordinates for enumerating primers.
        downstream_offset = feature_len if primer_allowable else feature_slop

        # Allowable coordinates will contain 1-based coordinates;
        # Note +1 on the 1-based stop coords includes the terminal base to be in the set (given python range behavior)
        if feature_strand is None:
            # If we don't know the strand of the primer, do not proceed
            raise RuntimeError("Absent strand information for %s" % FILE_DELIM.join(feature.fields))
        elif feature_strand == su.Strand.PLUS:
            allowable_coords = frozenset(range(feature.start - feature_slop + 1, feature.start + downstream_offset + 1))
        elif feature_strand == su.Strand.MINUS:
            allowable_coords = frozenset(range(feature.stop - downstream_offset + 1, feature.stop + feature_slop + 1))
        else:
            raise RuntimeError("Unrecognized strand for feature %s" % FILE_DELIM.join(feature.fields))

        feature_key = feature_name if use_name else feature_coords

        feature_dict[feature_key] = COORD_TUPLE(str(feature.chrom), feature.start, feature.stop,
                                                feature_name, feature_strand, feature_score, allowable_coords)

        observed_features.add(feature_coords)

    return feature_dict


def filter_gff(gff, key_file, attr_field, outfile, keep_records=("exon",)):
    """Filters a GFF/GTF by matching keys from a key file.

    :param str gff: file containing data to filter
    :param str key_file: file with keys, one per line
    :param str attr_field: GFF attribute field name for the key in the filter file
    :param str outfile: output file to write to
    :param tuple keep_records: records to keep
    """

    bt = pybedtools.BedTool(gff)
    keep_records_set = set(keep_records)

    with open(key_file, "r") as key_fh, \
            open(outfile, "w") as out_fh:

        keys = {k.strip() for k in key_fh}

        for feature in bt:

            if str(feature.fields[2]) not in keep_records_set:
                continue

            field_id = feature.attrs[attr_field].split(".")[0]

            if field_id in keys:
                out_fh.write(FILE_DELIM.join(feature.fields) + FILE_NEWLINE)
