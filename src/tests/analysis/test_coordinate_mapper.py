#!/usr/bin/env python3
"""Tests for analysis.coordinate_mapper"""

from shutil import copyfile
import tempfile
import unittest

import analysis.coordinate_mapper as cm
from analysis.references import faidx_ref
import core_utils.file_utils as fu
from satmut_utils.definitions import *

tempfile.tempdir = DEFAULT_TEMPDIR


class TestAminoAcidMapper(unittest.TestCase):
    """Tests for AminoAcidMapper."""

    GFF = "CBS_pEZY3.gff"
    GFF_REF = "CBS_pEZY3.fa"
    PRIMERS = "CBS_pEZY3_primers.bed"
    TARGETS = "CBS_pEZY3_targets.bed"
    APPRIS_REF = "GRCh38.chr21.fa.gz"
    CBS_REF = "CBS.fa"
    PKNOX1_REF = "PKNOX1.fa"

    @classmethod
    def setUpClass(cls):
        """Set up for TestAminoAcidMapper."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.abspath(os.path.join(cls.test_dir, "..", "test_data"))

        # Validate coordinate mapping for both a custom vector reference and for default APPRIS annotations

        # Custom reference
        cls.cbs_pezy3_ref = os.path.join(cls.test_data_dir, cls.GFF_REF)
        cls.cbs_pezy3_gff = os.path.join(cls.test_data_dir, cls.GFF)

        cls.cbs_pezy3_aa_mapper = cm.AminoAcidMapper(
            gff=cls.cbs_pezy3_gff, ref=cls.cbs_pezy3_ref, outdir=cls.tempdir, use_pickle=False, make_pickle=False,
            overwrite_pickle=False, mut_sig=cm.AminoAcidMapper.MUT_SIG_ANY, filter_unexpected=False)

        cls.cbs_pezy3_trx_seq = ""
        with open(cls.cbs_pezy3_ref, "r") as ref_fa:
            for i, line in enumerate(ref_fa):
                if i == 0:
                    continue
                cls.cbs_pezy3_trx_seq += line.rstrip(fu.FILE_NEWLINE)

        cls.cbs_pezy3_cds_seq = cls.cbs_pezy3_trx_seq[973:2629]

        # APPRIS transcripts on both strands require a genome reference
        cls.appris_ref_gz = os.path.join(cls.test_data_dir, cls.APPRIS_REF)
        cls.appris_ref_copy_gz = os.path.join(cls.tempdir, cls.APPRIS_REF)
        copyfile(cls.appris_ref_gz, cls.appris_ref_copy_gz)
        fu.gunzip_file(cls.appris_ref_copy_gz)
        cls.temp_ref_fa = fu.remove_extension(cls.appris_ref_copy_gz)
        faidx_ref(cls.temp_ref_fa)

        # This just has CBS and PKNOX1 annotations
        cls.appris_gff = os.path.join(cls.test_data_dir, GENCODE_TRX_GFF)

        # For sequence comparison use the sequences extracted from the transcriptome
        cls.cbs_ref = os.path.join(cls.test_data_dir, cls.CBS_REF)
        cls.pknox1_ref = os.path.join(cls.test_data_dir, cls.PKNOX1_REF)

        cls.appris_aa_mapper = cm.AminoAcidMapper(
            gff=cls.appris_gff, ref=cls.temp_ref_fa, outdir=cls.tempdir, use_pickle=False, make_pickle=False,
            overwrite_pickle=False, mut_sig="NNK", filter_unexpected=False)

        cls.appris_cbs_trx_seq = ""
        with open(cls.cbs_ref, "r") as ref_fa:
            for i, line in enumerate(ref_fa):
                if i == 0:
                    continue
                cls.appris_cbs_trx_seq += line.rstrip(fu.FILE_NEWLINE)

        cls.appris_cbs_cds_seq = cls.appris_cbs_trx_seq[260:1916]

        cls.appris_pknox1_trx_seq = ""
        with open(cls.pknox1_ref, "r") as ref_fa:
            for i, line in enumerate(ref_fa):
                if i == 0:
                    continue
                cls.appris_pknox1_trx_seq += line.rstrip(fu.FILE_NEWLINE)

        cls.appris_pknox1_cds_seq = cls.appris_pknox1_trx_seq[211:1522]

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestAminoAcidMapper."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

    def test_add_cds_info(self):
        """Tests that transcript information is added to the CDS info dict."""

        # Just test for the dict value, as we must pass in the trx_id key
        expected_val = (7108, 973, 2629, self.cbs_pezy3_trx_seq, self.cbs_pezy3_cds_seq)

        observed = {}
        self.cbs_pezy3_aa_mapper._add_cds_info(
            trx_id="CBS_pEZY3", transcript_cds_info=observed,
            trx_seq=self.cbs_pezy3_trx_seq, cds_seq=self.cbs_pezy3_cds_seq)

        self.assertEqual(expected_val, list(observed.values())[0])

    def test_get_cds_info_custom_ref(self):
        """Tests that CDS information is properly extracted from the default APPRIS reference GFF."""

        # Just test for the dict value, as we must pass in the trx_id key
        expected_val = (7108, 973, 2629, self.cbs_pezy3_trx_seq, self.cbs_pezy3_cds_seq)

        observed = self.cbs_pezy3_aa_mapper._get_cds_info()

        self.assertEqual(expected_val, list(observed.values())[0])

    def test_get_cds_info_appris_positive(self):
        """Tests that CDS info for a transcript on positive strand is extracted from the genome."""

        # Include the stop codon in the CDS coordinates
        expected_val = (5003, 211, 1522, self.appris_pknox1_trx_seq, self.appris_pknox1_cds_seq)

        observed = self.appris_aa_mapper._get_cds_info()

        self.assertEqual(expected_val, observed["ENST00000291547.9"])

    def test_get_cds_info_appris_negative(self):
        """Tests that CDS info for a transcript on negative strand is extracted from the genome."""

        expected_val = (2605, 260, 1916, self.appris_cbs_trx_seq, self.appris_cbs_cds_seq)

        observed = self.appris_aa_mapper._get_cds_info()

        self.assertEqual(expected_val, observed["ENST00000398165.7"])

    def test_concat_multi_mut_info(self):
        """Tests that MUT_INFO tuples with multiple codon position changes are concatenated into strings."""

        input_mut_info = cm.MUT_INFO_TUPLE(
            location="CDS", wt_codons=["ATG", "AAA"], mut_codons=["AGG", "AGA"],
            wt_aas=["M", "K"], mut_aas=["R", "R"], aa_changes=["p.M1R", "p.K2R"],
            aa_positions=[], matches_mut_sig=[True, False])

        expected = cm.MUT_INFO_TUPLE(
            location="CDS", wt_codons="ATG,AAA", mut_codons="AGG,AGA",
            wt_aas="M,K", mut_aas="R,R", aa_changes="p.M1R,p.K2R",
            aa_positions="1,2", matches_mut_sig="True,False")

        observed = self.appris_aa_mapper._concat_multi_mut_info(input_mut_info)

        self.assertEqual(expected, observed)

    def test_get_mut_info(self):
        """Tests that mutation info is properly generate for a WT and mutant CDS."""

        # NNK means A and C are not expected in the wobble position, but only if it does not match the REF
        expected = cm.MUT_INFO_TUPLE(location="CDS", wt_codons="ATG,AAA", mut_codons="AGG,AGA",
                                     wt_aas="M,K", mut_aas="R,R", aa_changes="p.M1R,p.K2R",
                                     aa_positions="1,2", matches_mut_sig="True,True")

        wt_cds_seq = "ATGAAA"
        mut_cds_seq = "AGGAGA"

        observed = self.appris_aa_mapper._get_mut_info(wt_cds_seq, mut_cds_seq)

        self.assertEqual(expected, observed)

    def test_get_mut_info_mut_sig(self):
        """Tests that mutation info is properly generate for a WT and mutant CDS with mut_sig mismatch."""

        # NNK means A and C are not expected in the wobble position, but only if it does not match the REF
        expected = cm.MUT_INFO_TUPLE(location="CDS", wt_codons="ATG,AAA", mut_codons="AGG,AGC",
                                     wt_aas="M,K", mut_aas="R,S", aa_changes="p.M1R,p.K2S",
                                     aa_positions="1,2", matches_mut_sig="True,False")

        wt_cds_seq = "ATGAAA"
        mut_cds_seq = "AGGAGC"

        observed = self.appris_aa_mapper._get_mut_info(wt_cds_seq, mut_cds_seq)

        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_trx_not_found(self):
        """Tests that a TranscriptNotFound Exception is raised if requested transcript is not in the input annotations."""

        with self.assertRaises(cm.TranscriptNotFound):
            self.appris_aa_mapper.get_codon_and_aa_changes(trx_id="ENSTnonexistent", pos=2, ref="T", alt="G")

    def test_get_codon_and_aa_changes_integenic_left(self):
        """Tests that an empty MUT_INFO tuple is returned with intergenic location if pos < 1 passed."""

        expected = cm.MUT_INFO_TUPLE(location=cm.AminoAcidMapper.INTERGENIC, **cm.AminoAcidMapper.DEFAULT_KWARGS)
        observed = self.appris_aa_mapper.get_codon_and_aa_changes(trx_id="ENST00000398165.7", pos=0, ref="T", alt="G")
        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_integenic_right(self):
        """Tests that an empty MUT_INFO tuple is returned with intergenic location if pos > trx len is passed."""

        expected = cm.MUT_INFO_TUPLE(location=cm.AminoAcidMapper.INTERGENIC, **cm.AminoAcidMapper.DEFAULT_KWARGS)
        observed = self.appris_aa_mapper.get_codon_and_aa_changes(trx_id="ENST00000398165.7", pos=22756, ref="T", alt="G")
        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_ref_mismatches(self):
        """Tests that an RuntimeError is raised when the passed REF does not match the transcript reference."""

        with self.assertRaises(RuntimeError):
            self.appris_aa_mapper.get_codon_and_aa_changes(trx_id="ENST00000398165.7", pos=260, ref="T", alt="G")

    def test_get_codon_and_aa_changes_fiveprime_untranslated(self):
        """Tests no protein annotations are returned for a variant in the 5' UTR."""

        expected = cm.MUT_INFO_TUPLE(location=cm.AminoAcidMapper.FIVEPRIME_UTR, **cm.AminoAcidMapper.DEFAULT_KWARGS)
        observed = self.appris_aa_mapper.get_codon_and_aa_changes(trx_id="ENST00000398165.7", pos=260, ref="C", alt="G")
        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_threeprime_untranslated(self):
        """Tests no protein annotations are returned for a variant in the 3' UTR."""

        expected = cm.MUT_INFO_TUPLE(location=cm.AminoAcidMapper.THREEPRIME_UTR, **cm.AminoAcidMapper.DEFAULT_KWARGS)
        observed = self.appris_aa_mapper.get_codon_and_aa_changes(trx_id="ENST00000398165.7", pos=1917, ref="A", alt="G")
        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_custom_negative_single(self):
        """Tests correct protein annotations are returned for a single-codon variant in CBS."""

        expected = cm.MUT_INFO_TUPLE(
            location=cm.AminoAcidMapper.CDS_ID, wt_codons="ATG", mut_codons="GGG", wt_aas="M", mut_aas="G",
            aa_changes="p.M1G", aa_positions="1", matches_mut_sig="True")

        observed = self.cbs_pezy3_aa_mapper.get_codon_and_aa_changes(
            trx_id="CBS_pEZY3", pos=974, ref="ATG", alt="GGG")

        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_appris_negative_single(self):
        """Tests correct protein annotations are returned for a single-codon variant in CBS."""

        expected = cm.MUT_INFO_TUPLE(
            location=cm.AminoAcidMapper.CDS_ID, wt_codons="ATG", mut_codons="GGG", wt_aas="M", mut_aas="G",
            aa_changes="p.M1G", aa_positions="1", matches_mut_sig="True")

        observed = self.appris_aa_mapper.get_codon_and_aa_changes(
            trx_id="ENST00000398165.7", pos=261, ref="ATG", alt="GGG")

        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_appris_negative_multi(self):
        """Tests correct protein annotations are returned for a multi-codon variant in CBS."""

        expected = cm.MUT_INFO_TUPLE(
            location=cm.AminoAcidMapper.CDS_ID, wt_codons="CCT,TCT", mut_codons="CCA,GGT", wt_aas="P,S", mut_aas="P,G",
            aa_changes="p.P2P,p.S3G", aa_positions="2,3", matches_mut_sig="False,True")

        observed = self.appris_aa_mapper.get_codon_and_aa_changes(
            trx_id="ENST00000398165.7", pos=266, ref="TTC", alt="AGG")

        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_appris_positive_single(self):
        """Tests correct protein annotations are returned for a single-codon variant in PKNOX1."""

        expected = cm.MUT_INFO_TUPLE(
            location=cm.AminoAcidMapper.CDS_ID, wt_codons="GCT", mut_codons="ACA", wt_aas="A", mut_aas="T",
            aa_changes="p.A3T", aa_positions="3", matches_mut_sig="False")

        observed = self.appris_aa_mapper.get_codon_and_aa_changes(
            trx_id="ENST00000291547.9", pos=218, ref="GCT", alt="ACA")

        self.assertEqual(expected, observed)
