#!/usr/bin/env/python
"""Utilities for extracting reference files from Ensembl IDs."""

import logging
import pybedtools
import pysam
import re

from analysis.aligners import BowtieConfig
import analysis.seq_utils as su
from analysis.coordinate_mapper import MapperBase
import core_utils.file_utils as fu
import core_utils.feature_file_utils as ffu
from definitions import *


APPRIS_CONTIG_DELIM = "|"

_logger = logging.getLogger(__name__)


class EnsemblIdNotFound(Exception):
    """Exception for when an Ensembl transcript or gene ID was not found in curated APPRIS reference files."""
    pass


class InvalidEnsemblId(Exception):
    """Exception for when an Ensembl transcript or gene ID is not in valid format."""
    pass


def ensembl_id_exists(reference_dir, ensembl_id):
    """Determines if ensembl_id is present in curated APPRIS reference files.

    :param str reference_dir: directory containing curated APPRIS reference files
    :param str ensembl_id: Ensembl gene or transcript ID, with version number
    :return tuple: (bool, str | None) if ensembl_id is in the references, alternative ID if present
    :raises InvalidEnsemblId: if the Ensembl ID does not start with ENST or ENSG
    """

    if re.search(MapperBase.ENSEMBL_TRX_PREFIX, ensembl_id):
        ids = {line.rstrip(fu.FILE_NEWLINE) for line in open(os.path.join(reference_dir, APPRIS_TRX_IDS), "r")}
    elif re.search(MapperBase.ENSEMBL_GENE_PREFIX, ensembl_id):
        ids = {line.rstrip(fu.FILE_NEWLINE) for line in open(os.path.join(reference_dir, APPRIS_GENE_IDS), "r")}
    else:
        raise InvalidEnsemblId("Ensembl ID %s was passed. Please provide a valid Ensembl identifier." % ensembl_id)

    if ensembl_id in ids:
        return True, None
    else:
        ensembl_id_base = ensembl_id.split(".")[0]
        # See if we have the base transcript but a different version
        for appris_id in ids:
            base_id = appris_id.split(".")[0]
            if ensembl_id_base == base_id:
                _logger.info(
                    "Ensembl ID %s was not found; however, curated reference files contain %s" % (ensembl_id, appris_id))
                return False, appris_id
        else:
            return False, None


def extract_fasta_reference(reference_dir, ensembl_id, outdir="."):
    """Extracts the reference FASTA given an Ensembl gene or transcript ID.

    :param str reference_dir: directory containing curated APPRIS reference files
    :param str ensembl_id: Ensembl gene or transcript ID, with version number
    :param str outdir: optional output directory for the reference files
    :return str: filepath of the output FASTA
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    output_fasta = os.path.join(outdir, fu.add_extension(ensembl_id, su.FASTA_FILETYPE))

    contig_id = [line.rstrip(fu.FILE_NEWLINE) for line in open(os.path.join(reference_dir, APPRIS_CONTIG_IDS), "r")
                 if ensembl_id in set(line.split(APPRIS_CONTIG_DELIM)[:2])][0]

    fa = pysam.faidx(os.path.join(reference_dir, APPRIS_TRX_FASTA), contig_id)

    with open(output_fasta, "w") as out_fasta_fh:
        out_fasta_fh.write(fa)

    return output_fasta


def extract_gff_reference(reference_dir, ensembl_id, outdir="."):
    """Extracts the reference GFF records given an Ensembl gene or transcript ID.

    :param str reference_dir: directory containing curated APPRIS reference files
    :param str ensembl_id: Ensembl gene or transcript ID, with version number
    :param str outdir: optional output directory for the reference files
    :return str: filepath of the output FASTA
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    ensembl_id_base = ensembl_id.split(".")[0]
    output_gff = os.path.join(outdir, fu.add_extension(ensembl_id_base, su.GFF_DEFAULT_EXT))

    save_features = {MapperBase.GENE_ID, MapperBase.TRX_ID, MapperBase.CDS_ID, MapperBase.START_CODON_ID,
                     MapperBase.STOP_CODON_ID, MapperBase.UTR_ID}

    default_trx_gff_bedtool = pybedtools.BedTool(os.path.join(reference_dir, GENCODE_TRX_GFF))

    with open(output_gff, "w") as out_gff_fh:
        for interval in default_trx_gff_bedtool:

            if interval[MapperBase.GFF_FEATURE_FIELD] in save_features and \
                    ensembl_id in {interval.attrs[ffu.GFF_ATTR_GENE_ID], interval.attrs[ffu.GFF_ATTR_TRANSCRIPT_ID]}:

                out_gff_fh.write(str(interval))

    return output_gff


def index_reference(ref):
    """samtools and bowtie2-indexes the reference FASTA.

    :param str ref: reference FASTA
    """

    # Need to index the reference with samtools and create a FM-index with bowtie2 if it has not been done
    if not os.path.exists(fu.add_extension(ref, su.FASTA_INDEX_SUFFIX)):
        _logger.info("Indexing FASTA %s." % ref)
        pysam.faidx(ref)

    # Build the bowtie2 index if it doesn't exist
    try:
        _ = BowtieConfig(ref=ref).test_build()
    except RuntimeError:
        _logger.info("Building FM index files for %s." % ref)
        # This exception is passed if the reference is not FM-indexed
        _ = BowtieConfig(ref=ref).build_fm_index()


def get_ensembl_references(reference_dir, ensembl_id, outdir="."):
    """Extracts reference files given an Ensembl gene or transcript ID.

    :param str reference_dir: directory containing curated APPRIS reference files
    :param str ensembl_id: Ensembl gene or transcript ID, with version number
    :param str outdir: optional output directory for the reference files
    :return tuple: (transcript_fasta, transcript_gff)
    :raises EnsemblIdNotFound: if the ensembl ID is not valid or the ID was not found in the APPRIS set
    """

    id_exists, alternative_id = ensembl_id_exists(reference_dir, ensembl_id)

    if not id_exists:
        if alternative_id is not None:
            raise EnsemblIdNotFound(
                "Ensembl ID %s was not found in the references; however, an alternative ID is %s." %
                (ensembl_id, alternative_id))
        else:
            raise EnsemblIdNotFound(
                "Ensembl ID %s was not found in the references. Please provide custom reference files.")

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    fa = extract_fasta_reference(reference_dir, ensembl_id, outdir)
    gff = extract_gff_reference(reference_dir, ensembl_id, outdir)
    index_reference(fa)

    return fa, gff
