#!/usr/bin/env/python
"""Utilities for extracting reference files from Ensembl IDs."""

import logging
import pybedtools
import pysam
import re

from analysis.aligners import BowtieConfig
import analysis.seq_utils as su
import core_utils.file_utils as fu
import core_utils.feature_file_utils as ffu
from definitions import *

ENSEMBL_GENE_PREFIX = "ENSG"
ENSEMBL_TRX_PREFIX = "ENST"
GFF_FEATURE_GENE = "gene"
GFF_FEATURE_TRX = "transcript"
GFF_FEATURE_EXON = "exon"
GFF_FEATURE_FIELD = 2

_logger = logging.getLogger(__name__)


def extract_fasta_reference(ensembl_id, outdir="."):
    """Extracts the reference FASTA given an Ensembl gene or transcript ID.

    :param str ensembl_id: Ensembl gene or transcript ID, with version number
    :param str outdir: optional output directory for the reference files
    :return str: filepath of the output FASTA
    """

    # Remove the version number
    ensembl_id_base = ensembl_id.split(".")[0]

    output_fasta = os.path.join(outdir, "{}.{}".format(ensembl_id_base, su.FASTA_FILETYPE))

    with open(APPRIS_TRX_FASTA, "r") as in_fasta_fh, \
            open(output_fasta, "w") as out_fasta_fh:

        keep_search = False

        for line in in_fasta_fh:
            if not keep_search and line.startswith(su.FASTA_HEADER_CHAR) and re.search(ensembl_id, line):
                out_fasta_fh.write(line)
                keep_search = True
                continue
            if keep_search:
                if line.startswith(su.FASTA_HEADER_CHAR):
                    break
                out_fasta_fh.write(line)

    return output_fasta


def extract_gff_reference(ensembl_id, outdir="."):
    """Extracts the reference GFF records given an Ensembl gene or transcript ID.

    :param str ensembl_id: Ensembl gene or transcript ID, with version number
    :param str outdir: optional output directory for the reference files
    :return str: filepath of the output FASTA
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    ensembl_id_base = ensembl_id.split(".")[0]
    output_gff = os.path.join(outdir, fu.add_extension(ensembl_id_base, su.GFF_DEFAULT_EXT))

    save_features = {GFF_FEATURE_GENE, GFF_FEATURE_TRX}
    default_trx_gff_bedtool = pybedtools.BedTool(APPRIS_TRX_GFF)

    with open(output_gff, "w") as out_gff_fh:

        for interval in default_trx_gff_bedtool:

            if interval[GFF_FEATURE_FIELD] in save_features and \
                    ensembl_id in {interval.attrs[ffu.GFF_ATTR_GENE_ID], interval.attrs[ffu.GFF_ATTR_TRANSCRIPT_ID]}:

                out_gff_fh.write(str(interval))

    return output_gff


def get_ensembl_references(ensembl_id, outdir="."):
    """Extracts reference files given an Ensembl gene or transcript ID.

    :param str ensembl_id: Ensembl gene or transcript ID, with version number
    :param str outdir: optional output directory for the reference files
    :return tuple: (transcript_fasta, transcript_gff)
    :raises RuntimeError: if the ensembl ID is not valid or the ID was not found in the APPRIS set
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    fa = extract_fasta_reference(ensembl_id, outdir)
    gff = extract_gff_reference(ensembl_id, outdir)

    if not os.path.exists(fu.add_extension(fa, su.FASTA_INDEX_SUFFIX)):
        _logger.info("Indexing FASTA %s" % fa)
        pysam.faidx(fa)

    try:
        _ = BowtieConfig(fa).test_build()
    except RuntimeError:
        BowtieConfig(fa).build_fm_index()

    return fa, gff
