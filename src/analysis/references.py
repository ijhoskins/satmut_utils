#!/usr/bin/env python3
"""Utilities for extracting reference files from Ensembl IDs."""

import logging
import pybedtools
import pysam
import re
import subprocess

from analysis.aligners import BowtieConfig
import analysis.seq_utils as su
from analysis.coordinate_mapper import MapperBase
import core_utils.file_utils as fu
import core_utils.feature_file_utils as ffu
from satmut_utils.definitions import *

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

APPRIS_CONTIG_DELIM = "|"
APPRIS_TRX_INDEX = 0
APPRIS_GENE_INDEX = 1

logger = logging.getLogger(__name__)


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
        with open(os.path.join(reference_dir, APPRIS_TRX_IDS), "r") as id_file:
            ids = {line.rstrip(fu.FILE_NEWLINE) for line in id_file}
    elif re.search(MapperBase.ENSEMBL_GENE_PREFIX, ensembl_id):
        with open(os.path.join(reference_dir, APPRIS_GENE_IDS), "r") as id_file:
            ids = {line.rstrip(fu.FILE_NEWLINE) for line in id_file}
    else:
        raise InvalidEnsemblId("Ensembl ID %s was passed. Please provide a valid Ensembl identifier." % ensembl_id)

    if ensembl_id in ids:
        return True, None
    else:
        ensembl_id_base = ensembl_id.split(".")[0]
        # See if we have the base ID but a different version
        for appris_id in ids:
            base_id = appris_id.split(".")[0]
            if ensembl_id_base == base_id:
                logger.info(
                    "Ensembl ID %s was not found; however, curated reference files contain %s" % (ensembl_id, appris_id))
                return False, appris_id
        else:
            return False, None


def get_trx_id_from_gene_id(reference_dir, ensembl_gene_id):
    """Determines if ensembl_id is present in curated APPRIS reference files.

    :param str reference_dir: directory containing curated APPRIS reference files
    :param str ensembl_gene_id: Ensembl gene ID, with version number
    :return str | None: transcript ID corresponding to the ensembl_gene_id; None if the gene ID not found
    """

    with open(os.path.join(reference_dir, APPRIS_CONTIG_IDS), "r") as contig_file:
        for line in contig_file:
            gene_id = line.split(APPRIS_CONTIG_DELIM)[APPRIS_GENE_INDEX]

            # Use the first transcript found for the gene ID
            if ensembl_gene_id == gene_id:
                return line.split(APPRIS_CONTIG_DELIM)[APPRIS_TRX_INDEX]
        else:
            return None


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

    with open(os.path.join(reference_dir, APPRIS_CONTIG_IDS), "r") as contig_id_file:
        contig_id = [line.rstrip(fu.FILE_NEWLINE) for line in contig_id_file if ensembl_id in set(
            line.split(APPRIS_CONTIG_DELIM)[:2])][0]

    logger.info("Extracting transcript reference FASTA for %s." % ensembl_id)
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

    logger.info("Extracting GFF transcript annotations for %s." % ensembl_id)
    output_gff = os.path.join(outdir, fu.add_extension(ensembl_id, su.GFF_DEFAULT_EXT))

    extract_call = ("fgrep", ensembl_id, os.path.join(reference_dir, GENCODE_TRX_GFF))

    with open(output_gff, "w") as out_gff:
        subprocess.call(extract_call, stdout=out_gff)

    return output_gff


def faidx_ref(ref):
    """samtools faidx a reference file.

    :param str ref: reference FASTA
    """

    # Occasional issues with pysam.faidx(ref)
    subprocess.call(("samtools", "faidx", ref))


def index_reference(ref):
    """samtools and bowtie2-indexes the reference FASTA.

    :param str ref: reference FASTA
    """

    # Need to index the reference with samtools and create a FM-index with bowtie2 if it has not been done
    if not os.path.exists(fu.add_extension(ref, su.FASTA_INDEX_SUFFIX)):
        logger.info("Generating FASTA index file for %s." % ref)
        faidx_ref(ref)

    # Build the bowtie2 index if it doesn't exist
    try:
        _ = BowtieConfig(ref=ref).test_build()
    except RuntimeError:
        logger.info("Building FM index files for %s." % ref)
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

    # Always select a transcript ID if a gene ID was passed; this ensures we get only one transcript reference
    # when filtering the GFF
    ensembl_trx_id = ensembl_id
    if re.search(MapperBase.ENSEMBL_GENE_PREFIX, ensembl_id):
        ensembl_trx_id = get_trx_id_from_gene_id(reference_dir=reference_dir, ensembl_gene_id=ensembl_id)
        if ensembl_trx_id is None:
            raise EnsemblIdNotFound(
                "No Ensembl transcript associated with Ensembl ID %s" % ensembl_id)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    fa = extract_fasta_reference(reference_dir, ensembl_id, outdir)
    gff = extract_gff_reference(reference_dir, ensembl_trx_id, outdir)
    index_reference(fa)

    return fa, gff
