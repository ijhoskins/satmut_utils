#!/usr/bin/env python3
"""Tests for analysis.references"""

import shutil
import tempfile
import unittest

from analysis.references import *
from analysis.seq_utils import FASTA_INDEX_SUFFIX
from satmut_utils.definitions import *

tempfile.tempdir = DEFAULT_TEMPDIR


class TestReferences(unittest.TestCase):
    """Tests for reference extraction methods."""

    GENOME_REF = "GRCh38.chr21.fa.gz"
    CBS_REF = "CBS.fa"
    CBS_GFF = "CBS.gff"

    @classmethod
    def setUpClass(cls):
        """Set up for TestReferences."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.join(os.path.split(cls.test_dir)[0], "test_data")

        # This is a gzipped chr21
        cls.genome_ref_gz = os.path.join(cls.test_data_dir, cls.GENOME_REF)

        # This just has CBS and PKNOX1 annotations
        cls.gencode_gff = os.path.join(cls.test_data_dir, GENCODE_TRX_GFF)

        # For comparison to extracted GFF
        cls.cbs_gff = os.path.join(cls.test_data_dir, cls.CBS_GFF)

        cls.cbs_ref = os.path.join(cls.test_data_dir, cls.CBS_REF)
        cls.cbs_ref_copy = tempfile.NamedTemporaryFile(suffix=".copy.fa", delete=False, dir=cls.tempdir).name
        shutil.copyfile(cls.cbs_ref, cls.cbs_ref_copy)

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestReferences."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

    def test_ensembl_id_exists_invalid_id(self):
        """Tests that an InvalidEnsemblId exception is raised when an invalid ID is passed."""

        with self.assertRaises(InvalidEnsemblId):
            ensembl_id_exists(reference_dir=self.test_data_dir, ensembl_id="NM_000071.3")

    def test_ensembl_id_exists_true(self):
        """Tests that an Ensembl ID in the curated set is identified."""

        expected = True, None
        observed = ensembl_id_exists(reference_dir=self.test_data_dir, ensembl_id="ENST00000398165.7")
        self.assertEqual(expected, observed)

    def test_ensembl_id_exists_false(self):
        """Tests that an Ensembl ID not in the curated set is not identified."""

        expected = False, None
        observed = ensembl_id_exists(reference_dir=self.test_data_dir, ensembl_id="ENST1.1")
        self.assertEqual(expected, observed)

    def test_ensembl_id_exists_false_with_alternative(self):
        """Tests that an Ensembl ID not in the curated set is not identified, but an alternative returned."""

        expected = False, "ENST00000398165.7"
        observed = ensembl_id_exists(reference_dir=self.test_data_dir, ensembl_id="ENST00000398165.6")
        self.assertEqual(expected, observed)

    def test_extract_fasta_reference(self):
        """Tests that a FASTA with a single transcript is returned from an index FASTA."""

        extracted_fa = extract_fasta_reference(
            reference_dir=self.test_data_dir, ensembl_id="ENST00000398165.7", outdir=self.tempdir)

        with open(extracted_fa, "r") as observed_fa, \
                open(self.cbs_ref, "r") as expected_fa:

            observed = ""
            expected = ""

            # Compare the strings only as line width may not be static
            # -n default changed from 50 to 60 between samtools versions?
            for line in observed_fa:
                observed += line.strip(fu.FILE_NEWLINE)

            for line in expected_fa:
                expected += line.strip(fu.FILE_NEWLINE)

        self.assertEqual(expected, observed)

    def test_extract_gff_reference(self):
        """Tests that GFF annotations for a single transcript are properly filtered from the transcriptome GFF."""

        extracted_gff = extract_gff_reference(
            reference_dir=self.test_data_dir, ensembl_id="ENST00000398165.7", outdir=self.tempdir)

        # Just compare the non-attr fields
        with open(extracted_gff, "r") as observed_gff, \
                open(self.cbs_gff, "r") as expected_gff:

            # For now just check the non-attribute fields are equivalent
            # pybedtools seems to uniquify attributes that start with the same keyname (e.g. "tag")
            observed = [fu.FILE_DELIM.join(e.split(fu.FILE_DELIM)[0:8]) for e in observed_gff.readlines()]
            expected = [fu.FILE_DELIM.join(e.split(fu.FILE_DELIM)[0:8]) for e in expected_gff.readlines()]

        self.assertEqual(expected, observed)

    def test_index_reference_samtools(self):
        """Tests that a samtools FASTA index file is generated."""

        index_reference(self.cbs_ref_copy)
        self.assertTrue(os.path.exists(fu.add_extension(self.cbs_ref_copy, FASTA_INDEX_SUFFIX)))

    def test_index_reference_bowtie(self):
        """Tests that bowtie2 FM index files are generated."""

        index_reference(self.cbs_ref_copy)

        try:
            BowtieConfig(self.cbs_ref_copy).test_build()
        except RuntimeError:
            self.fail("bowtie2 FM index files were not generated for reference %s" % self.cbs_ref_copy)

    def test_get_ensembl_references_id_not_found_with_alternative(self):
        """Tests that an EnsemblIdNotFound exception is raised if passed an Ensembl ID not in the curated set."""

        with self.assertRaises(EnsemblIdNotFound):
            _, _ = get_ensembl_references(
                reference_dir=self.test_data_dir, ensembl_id="ENST00000398165.6", outdir=self.tempdir)

    def test_get_ensembl_references_id_not_found_no_alternative(self):
        """Tests that an EnsemblIdNotFound exception is raised if passed an Ensembl ID not in the curated set."""

        with self.assertRaises(EnsemblIdNotFound):
            _, _ = get_ensembl_references(reference_dir=self.test_data_dir, ensembl_id="ENST1.1", outdir=self.tempdir)

    # Don't make smoke test of get_ensembl_references() as we test extraction of both FASTA and GFF as well as indexing
