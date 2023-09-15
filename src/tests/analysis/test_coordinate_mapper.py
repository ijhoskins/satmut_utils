#!/usr/bin/env python3
"""Tests for analysis.coordinate_mapper"""

from shutil import copyfile
import tempfile
import unittest

import analysis.coordinate_mapper as cm
from analysis.references import faidx_ref
import core_utils.file_utils as fu
from core_utils.vcf_utils import VariantType
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

        cls.ref_codon_dict = cls.appris_aa_mapper._get_pos_codon_dict("ATGCCTTCTGAG")

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
            location="CDS", wt_codons=("ATG", "AAA",), mut_codons=("AGG", "AGA",),
            wt_aas=("M", "K",), mut_aas=("R", "R",), aa_changes=("p.M1R", "p.K2R",),
            aa_positions=(), matches_mut_sig=(True, False,),
            mave_hgvs_nt="ENST00000398165.7:c.2t>g,ENST00000398165.7:c.5a>g",
            mave_hgvs_tx="ENST00000398165.7:r.2u>g,ENST00000398165.7:r.5a>g",
            mave_hgvs_pro="p.Met1Arg,p.Lys2Arg")

        expected = cm.MUT_INFO_TUPLE(
            location="CDS", wt_codons="ATG,AAA", mut_codons="AGG,AGA",
            wt_aas="M,K", mut_aas="R,R", aa_changes="p.M1R,p.K2R",
            aa_positions="1,2", matches_mut_sig="True,False",
            mave_hgvs_nt="ENST00000398165.7:c.2t>g,ENST00000398165.7:c.5a>g",
            mave_hgvs_tx="ENST00000398165.7:r.2u>g,ENST00000398165.7:r.5a>g",
            mave_hgvs_pro="p.Met1Arg,p.Lys2Arg")

        observed = self.appris_aa_mapper._concat_multi_mut_info(input_mut_info)

        self.assertEqual(expected, observed)

    def test_get_pos_codon_dict(self):
        """Tests that a CDS sequence can be indexed."""

        observed = self.appris_aa_mapper._get_pos_codon_dict("ATGCCT")

        expected = dict(
            zip(range(0, 6),
                (cm.CODON_TUPLE(1, "ATG", 0), cm.CODON_TUPLE(1, "ATG", 1), cm.CODON_TUPLE(1, "ATG", 2),
                 cm.CODON_TUPLE(2, "CCT", 0), cm.CODON_TUPLE(2, "CCT", 1), cm.CODON_TUPLE(2, "CCT", 2))))

        self.assertDictEqual(expected, observed)

    def test_get_indel_downstream_remainder_ins(self):
        """Tests for correct downstream codons and amino acids following an insertion frameshift."""

        observed = self.appris_aa_mapper._get_indel_downstream_remainder(
            trx_seq=self.appris_cbs_trx_seq, pos=265, ref_len=1, alt="CT", curr_codon="CCT", base_index=1,
            var_type=VariantType.INS)

        # Remaining codons and remaining amino acids
        expected = (("CCT", "TTC", "TGA",), ("P", "F", "*"))

        self.assertEqual(expected, observed)

    def test_get_indel_downstream_remainder_del(self):
        """Tests for correct downstream codons and amino acids following an deletion frameshift."""

        observed = self.appris_aa_mapper._get_indel_downstream_remainder(
            trx_seq=self.appris_cbs_trx_seq, pos=263, ref_len=3, alt="G", curr_codon="ATG", base_index=2,
            var_type=VariantType.DEL)

        # Remaining codons and remaining amino acids
        expected = (("ATG", "TTC", "TGA",), ("M", "F", "*"))

        self.assertEqual(expected, observed)

    def test_annotate_snp(self):
        """Tests for correct custom annotation of a SNP."""

        observed = self.appris_aa_mapper._annotate_snp(
            ref_codon_dict=self.ref_codon_dict,
            alt_codon_dict=self.appris_aa_mapper._get_pos_codon_dict("ATGCCCTCTGAG"), start_index=5)

        expected = (("CCT",), ("CCC",), ("P",), ("P",), ("p.P2P",), (False,))

        self.assertEqual(expected, observed)

    def test_annotate_mnp(self):
        """Tests for correct custom annotation of a MNP."""

        observed = self.appris_aa_mapper._annotate_mnp(
            ref_codon_dict=self.ref_codon_dict,
            alt_codon_dict=self.appris_aa_mapper._get_pos_codon_dict("ATGCCCACTGAG"),
            start_index=5, end_index=6, cds_len=3315)

        expected = (("CCT", "TCT",), ("CCC", "ACT",), ("P", "S",), ("P", "T",), ("p.P2P", "p.S3T",), (False, True,))

        self.assertEqual(expected, observed)

    def test_annotate_ins(self):
        """Tests for correct custom annotation of an out-of-frame insertion."""

        observed = self.appris_aa_mapper._annotate_ins(
            trx_seq=self.appris_cbs_trx_seq, ref_codon_dict=self.ref_codon_dict, start_index=4,
            pos=265, ref="C", alt="CT")

        expected = (("CCT",), ("CCT", "TTC", "TGA",), ("P",), ("P", "F", "*",), ("p.P2P,F,*",), (False,))

        self.assertEqual(expected, observed)

    def test_annotate_ins_inframe(self):
        """Tests for correct custom annotation of an in-frame insertion."""

        observed = self.appris_aa_mapper._annotate_ins(
            trx_seq=self.appris_cbs_trx_seq, ref_codon_dict=self.ref_codon_dict, start_index=5,
            pos=266, ref="T", alt="TATG")

        expected = (("CCT",), ("CCT", "ATG",), ("P",), ("P", "M",), ("p.P2P,M",), (True,))

        self.assertEqual(expected, observed)

    def test_annotate_del(self):
        """Tests for correct custom annotation of a deletion."""

        observed = self.appris_aa_mapper._annotate_del(
            trx_seq=self.appris_cbs_trx_seq, ref_codon_dict=self.ref_codon_dict, start_index=2,
            pos=263, ref="GCC", alt="G")

        expected = (("ATG",), ("ATG", "TTC", "TGA",), ("M",), ("M", "F", "*",), ("p.M1M,F,*",), (False,))

        self.assertEqual(expected, observed)

    def test_annotate_stopvar(self):
        """Tests for correct custom annotation of a nonstop variant."""

        ref_codon_dict = self.appris_aa_mapper._get_pos_codon_dict(self.appris_cbs_cds_seq)

        # TGA>TGT SNP
        alt_codon_dict = self.appris_aa_mapper._get_pos_codon_dict(self.appris_cbs_cds_seq[:-1] + "T")

        observed = self.appris_aa_mapper._annotate_stopvar(
            trx_seq=self.appris_cbs_trx_seq, pos=1916, ref="A", alt="T", cds_start_offset=260, cds_len=3315,
            ref_codon_dict=ref_codon_dict, alt_codon_dict=alt_codon_dict)

        expected = (("TGA",), ("TGT",), ("*",), ("C",), ("p.*552C",), (True,))

        self.assertEqual(expected, observed)

    def test_annotate_snp_hgvs_cds(self):
        """Tests for correct MAVE-HGVS annotation of a SNP (substitution) in the CDS."""

        observed = self.appris_aa_mapper._annotate_snp_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, pos=268, ref="C", alt="T",
            cds_start_offset=260, cds_stop_offset=1916, start_index=7, ref_aa="S", aa_pos=3, alt_aa="F")

        expected = ("ENST00000398165.7:c.8c>t", "ENST00000398165.7:r.8c>u", "p.Ser3Phe",)

        self.assertEqual(expected, observed)

    def test_annotate_snp_hgvs_utr5(self):
        """Tests for correct MAVE-HGVS annotation of a SNP (substitution) in the 5' UTR."""

        observed = self.appris_aa_mapper._annotate_snp_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.FIVEPRIME_UTR, pos=259, ref="G", alt="T",
            cds_start_offset=260, cds_stop_offset=1916, start_index=None, ref_aa=None, aa_pos=None, alt_aa=None)

        expected = ("ENST00000398165.7:c.-2g>t", "ENST00000398165.7:r.-2g>u", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_snp_hgvs_utr3(self):
        """Tests for correct MAVE-HGVS annotation of a SNP (substitution) in the 3' UTR."""

        observed = self.appris_aa_mapper._annotate_snp_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.THREEPRIME_UTR, pos=1920, ref="C", alt="T",
            cds_start_offset=260, cds_stop_offset=1916, start_index=None, ref_aa=None, aa_pos=None, alt_aa=None)

        expected = ("ENST00000398165.7:c.*4c>t", "ENST00000398165.7:r.*4c>u", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_mnp_hgvs_cds_single_codon(self):
        """Tests for correct MAVE-HGVS annotation of a MNP (delins) affecting a single codon."""

        alt_codon_dict = self.appris_aa_mapper._get_pos_codon_dict("ATGCCTAAAGAG")

        observed = self.appris_aa_mapper._annotate_mnp_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, pos=267, alt="AAA", cds_start_offset=260,
            cds_stop_offset=1916, cds_len=3315, start_index=6, ref_codon_dict=self.ref_codon_dict,
            alt_codon_dict=alt_codon_dict)

        expected = ("ENST00000398165.7:c.7_9delinsaaa", "ENST00000398165.7:r.7_9delinsaaa", "p.Ser3delinsLys",)

        self.assertEqual(expected, observed)

    def test_annotate_mnp_hgvs_cds_multi_codon(self):
        """Tests for correct MAVE-HGVS annotation of a MNP (delins) affecting two codons."""

        alt_codon_dict = self.appris_aa_mapper._get_pos_codon_dict("ATGCCTTCACAG")

        observed = self.appris_aa_mapper._annotate_mnp_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, pos=269, alt="AC", cds_start_offset=260,
            cds_stop_offset=1916, cds_len=3315, start_index=8, ref_codon_dict=self.ref_codon_dict,
            alt_codon_dict=alt_codon_dict)

        expected = ("ENST00000398165.7:c.9_10delinsac", "ENST00000398165.7:r.9_10delinsac", "p.Ser3_Glu4delinsSerGln",)

        self.assertEqual(expected, observed)

    def test_annotate_mnp_hgvs_utr5_multi(self):
        """Tests for correct MAVE-HGVS annotation of a MNP (delins) in the 5' UTR."""

        observed = self.appris_aa_mapper._annotate_mnp_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.FIVEPRIME_UTR, pos=255, alt="GCG",
            cds_start_offset=260, cds_stop_offset=1916, cds_len=3315, start_index=None,
            ref_codon_dict=self.ref_codon_dict, alt_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.-6_-4delinsgcg", "ENST00000398165.7:r.-6_-4delinsgcg", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_mnp_hgvs_utr3_multi(self):
        """Tests for correct MAVE-HGVS annotation of a MNP (delins) in the 3' UTR."""

        observed = self.appris_aa_mapper._annotate_mnp_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.THREEPRIME_UTR, pos=1918, alt="TC",
            cds_start_offset=260, cds_stop_offset=1916, cds_len=3315, start_index=None,
            ref_codon_dict=self.ref_codon_dict, alt_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.*2_*3delinstc", "ENST00000398165.7:r.*2_*3delinsuc", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_group_mismatches(self):
        """Tests that mismatches in multivariants are appropriately grouped to determine substitution versus delins."""

        mm_positions = (500, 503, 504, 510, 511, 512, 520, 522)
        mm_indices = (0, 3, 4, 10, 11, 12, 20, 22)
        observed = self.appris_aa_mapper._group_mismatches(mm_positions, mm_indices)

        expected = ([(0, 500)], [(3, 503), (4, 504)], [(10, 510), (11, 511), (12, 512)], [(20, 520), (22, 522)],)

        self.assertEqual(expected, observed)

    def test_annotate_multivariant_hgvs_cds(self):
        """Tests that a multivariant is correctly annotated in the CDS."""

        alt_codon_dict = self.appris_aa_mapper._get_pos_codon_dict("ATGTTTTGTGAG")

        observed = self.appris_aa_mapper._annotate_multivariant_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, pos=264, ref="CCTTC",
            alt="TTTTG", cds_start_offset=260, cds_stop_offset=1916, cds_len=3315, ref_codon_dict=self.ref_codon_dict,
            alt_codon_dict=alt_codon_dict)

        expected = ("ENST00000398165.7:c.4_5delinstt,ENST00000398165.7:c.8c>g",
                    "ENST00000398165.7:r.4_5delinsuu,ENST00000398165.7:r.8c>g",
                    "p.Pro2delinsPhe,p.Ser3Cys",)

        self.assertEqual(expected, observed)

    def test_annotate_multivariant_hgvs_utr5_isolated(self):
        """Tests that a multivariant in the 5' UTR and not overlapping the CDS is correctly annotated."""

        # Need to pass codon dictionaries as the might be needed for 5' UTR variants, although here they are not req'd
        observed = self.appris_aa_mapper._annotate_multivariant_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.FIVEPRIME_UTR, pos=255, ref="CCCAGC",
            alt="TCCGAC", cds_start_offset=260, cds_stop_offset=1916, cds_len=3315, ref_codon_dict=self.ref_codon_dict,
            alt_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.-6c>t,ENST00000398165.7:c.-3_-2delinsga",
                    "ENST00000398165.7:r.-6c>u,ENST00000398165.7:r.-3_-2delinsga", "p.(=),p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_multivariant_hgvs_utr5_overlapping(self):
        """Tests that a multivariant in the 5' UTR and overlapping the CDS is correctly annotated."""

        alt_codon_dict = self.appris_aa_mapper._get_pos_codon_dict("TGGCCTTCTGAG")

        observed = self.appris_aa_mapper._annotate_multivariant_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.FIVEPRIME_UTR, pos=258, ref="AGCAT",
            alt="TGCTG", cds_start_offset=260, cds_stop_offset=1916, cds_len=3315, ref_codon_dict=self.ref_codon_dict,
            alt_codon_dict=alt_codon_dict)

        expected = ("ENST00000398165.7:c.-3a>t,ENST00000398165.7:c.1_2delinstg",
                    "ENST00000398165.7:r.-3a>u,ENST00000398165.7:r.1_2delinsug",
                    "p.(=),p.(=),p.Met1delinsTrp",)

        # For protein annotation, the first p.(=) is for the subst, and the second is associated with the MNP
        # that spans the UTR/CDS junction

        self.assertEqual(expected, observed)

    def test_annotate_multivariant_hgvs_utr3_isolated(self):
        """Tests that a multivariant in the 3' UTR and not overlapping the CDS is correctly annotated."""

        observed = self.appris_aa_mapper._annotate_multivariant_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.THREEPRIME_UTR, pos=1920, ref="CCGGAG",
            alt="TCGCCG", cds_start_offset=260, cds_stop_offset=1916, cds_len=3315, ref_codon_dict=self.ref_codon_dict,
            alt_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.*4c>t,ENST00000398165.7:c.*7_*8delinscc",
                    "ENST00000398165.7:r.*4c>u,ENST00000398165.7:r.*7_*8delinscc", "p.(=),p.(=)",)

        self.assertEqual(expected, observed)

    def test_is_duplication_true_single_nt(self):
        """Tests that a single-nt duplication is detected."""

        observed = self.appris_aa_mapper._is_duplication(trx_seq=self.appris_cbs_trx_seq, pos=265, alt="CC")
        self.assertTrue(observed)

    def test_is_duplication_false_single_nt(self):
        """Tests that a normal insertion is not determined to be a duplication."""

        observed = self.appris_aa_mapper._is_duplication(trx_seq="ENST00000398165.7", pos=265, alt="CT")
        self.assertFalse(observed)

    def test_is_duplication_true_three_nt(self):
        """Tests that a three-nt duplication is detected."""

        observed = self.appris_aa_mapper._is_duplication(trx_seq=self.appris_cbs_trx_seq, pos=265, alt="CGCC")
        self.assertTrue(observed)

    def test_is_duplication_true_four_nt(self):
        """Tests that a four-nt duplication is detected."""

        observed = self.appris_aa_mapper._is_duplication(trx_seq=self.appris_cbs_trx_seq, pos=265, alt="CTGCC")
        self.assertTrue(observed)

    def test_get_fs_first_alt_aa_ins(self):
        """Tests that the first amino acid not matching the reference amino acid for a frameshift ins is determined."""

        observed = self.appris_aa_mapper._get_fs_first_alt_aa(
            alt_aas=("P", "F", "*"), cds_stop_offset=1916, start_index=4, ref_codon_dict=self.ref_codon_dict)

        expected = ("S", 3)

        self.assertEqual(expected, observed)

    def test_get_fs_first_alt_aa_del(self):
        """Tests that the first amino acid not matching the reference amino acid for a frameshift del is determined."""

        observed = self.appris_aa_mapper._get_fs_first_alt_aa(
            alt_aas=("F", "*"), cds_stop_offset=1916, start_index=2, ref_codon_dict=self.ref_codon_dict)

        expected = ("P", 2)

        self.assertEqual(expected, observed)

    def test_annotate_dup_hgvs_cds_fs(self):
        """Tests for correct MAVE-HGVS annotation of a single-nt duplication in the CDS."""

        observed = self.appris_aa_mapper._annotate_dup_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, pos=265, alt_len=2, cds_start_offset=260,
            cds_stop_offset=1916, start_index=4, alt_aas=("P", "F", "*"), ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.5dup", "ENST00000398165.7:r.5dup", "p.Ser3fs",)

        self.assertEqual(expected, observed)

    def test_annotate_dup_hgvs_cds_inframe_single(self):
        """Tests for correct MAVE-HGVS annotation of an in-frame duplication of a single amino acid in the CDS."""

        observed = self.appris_aa_mapper._annotate_dup_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, pos=266, alt_len=4, cds_start_offset=260,
            cds_stop_offset=1916, start_index=5, alt_aas=("P", "P"), ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.4_6dup", "ENST00000398165.7:r.4_6dup", "p.Pro2dup",)

        self.assertEqual(expected, observed)

    def test_annotate_dup_hgvs_cds_inframe_multi(self):
        """Tests for correct MAVE-HGVS annotation of an in-frame duplication of a single amino acid in the CDS."""

        observed = self.appris_aa_mapper._annotate_dup_hgvs(
            trx_id="ENST00000398165.7", alt_len=7, location=self.appris_aa_mapper.CDS_ID, pos=266, cds_start_offset=260,
            cds_stop_offset=1916, start_index=5, alt_aas=("P", "M", "P"), ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.1_6dup", "ENST00000398165.7:r.1_6dup", "p.Met1_Pro2dup",)

        self.assertEqual(expected, observed)

    def test_annotate_dup_hgvs_utr5(self):
        """Tests for correct MAVE-HGVS annotation of a duplication in the 5' UTR."""

        observed = self.appris_aa_mapper._annotate_dup_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.FIVEPRIME_UTR, pos=255, alt_len=4,
            cds_start_offset=260, cds_stop_offset=1916, start_index=None, alt_aas=None, ref_codon_dict=None)

        # 255:C:CGTC
        expected = ("ENST00000398165.7:c.-8_-6dup", "ENST00000398165.7:r.-8_-6dup", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_dup_hgvs_utr3(self):
        """Tests for correct MAVE-HGVS annotation of a duplication in the 3' UTR."""

        observed = self.appris_aa_mapper._annotate_dup_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.THREEPRIME_UTR, pos=1920, alt_len=4,
            cds_start_offset=260, cds_stop_offset=1916, start_index=None, alt_aas=None, ref_codon_dict=None)

        # 1920:C:CGTC
        expected = ("ENST00000398165.7:c.*2_*4dup", "ENST00000398165.7:r.*2_*4dup", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_dup_hgvs_cds_utr3(self):
        """Tests for correct MAVE-HGVS annotation of a duplication spanning the CDS and 3' UTR."""

        observed = self.appris_aa_mapper._annotate_dup_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.THREEPRIME_UTR, pos=1917, alt_len=4,
            cds_start_offset=260, cds_stop_offset=1916, start_index=None, alt_aas=None, ref_codon_dict=None)

        # 1917:A:AGAA
        expected = ("ENST00000398165.7:c.1915_*1dup", "ENST00000398165.7:r.1915_*1dup", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_ins_hgvs_cds_inframe_single(self):
        """Tests for correct MAVE-HGVS annotation of an in-frame, single amino acid insertion in the CDS."""

        observed = self.appris_aa_mapper._annotate_ins_hgvs(
            trx_id="ENST00000398165.7", trx_seq=self.appris_cbs_trx_seq, location=self.appris_aa_mapper.CDS_ID,
            pos=266, alt="TATG", cds_start_offset=260, cds_stop_offset=1916, start_index=5, alt_aas=("P", "M",),
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.6_7insatg", "ENST00000398165.7:r.6_7insaug", "p.Pro2_Ser3insMet")

        self.assertEqual(expected, observed)

    def test_annotate_ins_hgvs_cds_inframe_multi(self):
        """Tests for correct MAVE-HGVS annotation of an in-frame, multi amino acid insertion in the CDS."""

        observed = self.appris_aa_mapper._annotate_ins_hgvs(
            trx_id="ENST00000398165.7", trx_seq=self.appris_cbs_trx_seq, location=self.appris_aa_mapper.CDS_ID,
            pos=266, alt="TATGAAA", cds_start_offset=260, cds_stop_offset=1916, start_index=5, alt_aas=("P", "M", "K"),
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.6_7insatgaaa", "ENST00000398165.7:r.6_7insaugaaa", "p.Pro2_Ser3insMetLys",)

        self.assertEqual(expected, observed)

    def test_annotate_ins_hgvs_utr5(self):
        """Tests for correct MAVE-HGVS annotation of an insertion in the 5' UTR."""

        observed = self.appris_aa_mapper._annotate_ins_hgvs(
            trx_id="ENST00000398165.7", trx_seq=self.appris_cbs_trx_seq, location=self.appris_aa_mapper.FIVEPRIME_UTR,
            pos=255, alt="CT", cds_start_offset=260, cds_stop_offset=1916, start_index=None, alt_aas=None,
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.-6_-5inst", "ENST00000398165.7:r.-6_-5insu", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_ins_hgvs_utr3(self):
        """Tests for correct MAVE-HGVS annotation of an insertion in the 3' UTR."""

        observed = self.appris_aa_mapper._annotate_ins_hgvs(
            trx_id="ENST00000398165.7", trx_seq=self.appris_cbs_trx_seq, location=self.appris_aa_mapper.THREEPRIME_UTR,
            pos=1920, alt="CT", cds_start_offset=260, cds_stop_offset=1916, start_index=None, alt_aas=None,
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.*4_*5inst", "ENST00000398165.7:r.*4_*5insu", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_del_hgvs_cds_inframe_single(self):
        """Tests for correct MAVE-HGVS annotation of an in-frame, single amino acid deletion in the CDS."""

        observed = self.appris_aa_mapper._annotate_del_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, pos=266, ref="TTCT",
            cds_start_offset=260, cds_stop_offset=1916, start_index=5, alt_aas=("P",),
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.7_9del", "ENST00000398165.7:r.7_9del", "p.Ser3del",)

        self.assertEqual(expected, observed)

    def test_annotate_del_hgvs_cds_inframe_multi(self):
        """Tests for correct MAVE-HGVS annotation of an in-frame, multi amino acid deletion in the CDS."""

        observed = self.appris_aa_mapper._annotate_del_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, pos=266, ref="TTCTGAG",
            cds_start_offset=260, cds_stop_offset=1916, start_index=5, alt_aas=("P",),
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.7_12del", "ENST00000398165.7:r.7_12del", "p.Ser3_Glu4del",)

        self.assertEqual(expected, observed)

    def test_annotate_del_hgvs_utr5(self):
        """Tests for correct MAVE-HGVS annotation of an isolated deletion in the 5' UTR."""

        observed = self.appris_aa_mapper._annotate_del_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.FIVEPRIME_UTR, pos=255, ref="CCC",
            cds_start_offset=260, cds_stop_offset=1916, start_index=None, alt_aas=None,
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.-5_-4del", "ENST00000398165.7:r.-5_-4del", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_del_hgvs_utr5_cds(self):
        """Tests for correct MAVE-HGVS annotation of a deletion spanning the 5' UTR - CDS junction."""

        observed = self.appris_aa_mapper._annotate_del_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.FIVEPRIME_UTR, pos=259, ref="GCATGCCT",
            cds_start_offset=260, cds_stop_offset=1916, start_index=None, alt_aas=None,
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.-1_6del", "ENST00000398165.7:r.-1_6del", "p.Met1_Pro2del",)

        self.assertEqual(expected, observed)

    def test_annotate_del_hgvs_utr3(self):
        """Tests for correct MAVE-HGVS annotation of an isolated deletion in the 3' UTR."""

        observed = self.appris_aa_mapper._annotate_del_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.THREEPRIME_UTR, pos=1920, ref="CCG",
            cds_start_offset=260, cds_stop_offset=1916, start_index=None, alt_aas=None,
            ref_codon_dict=self.ref_codon_dict)

        expected = ("ENST00000398165.7:c.*5_*6del", "ENST00000398165.7:r.*5_*6del", "p.(=)",)

        self.assertEqual(expected, observed)

    def test_annotate_stopvar_hgvs_synon(self):
        """Tests for correct MAVE-HGVS annotation of a synonymous SNP in the stop codon."""

        ref_codon_dict = self.appris_aa_mapper._get_pos_codon_dict(self.appris_cbs_cds_seq)
        alt_codon_dict = self.appris_aa_mapper._get_pos_codon_dict(self.appris_cbs_cds_seq[:-2] + "AA")

        # 1915:G:A synonymous SNP
        observed = self.appris_aa_mapper._annotate_stopvar_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, trx_seq=self.appris_cbs_trx_seq,
            pos=1915, ref="G", alt="A", ref_len=1, alt_len=1, cds_start_offset=260, cds_stop_offset=1916, cds_len=3315,
            start_index=1654, end_index=1654, ref_codon_dict=ref_codon_dict, alt_codon_dict=alt_codon_dict)

        expected = ("ENST00000398165.7:c.1655g>a", "ENST00000398165.7:r.1655g>a", "p.Ter552=",)

        self.assertEqual(expected, observed)

    def test_annotate_stopvar_hgvs_fs(self):
        """Tests for correct MAVE-HGVS annotation of a frameshift at the stop codon."""

        ref_codon_dict = self.appris_aa_mapper._get_pos_codon_dict(self.appris_cbs_cds_seq)
        alt_codon_dict = self.appris_aa_mapper._get_pos_codon_dict(self.appris_cbs_cds_seq[:-1] + "T")

        observed = self.appris_aa_mapper._annotate_stopvar_hgvs(
            trx_id="ENST00000398165.7", location=self.appris_aa_mapper.CDS_ID, trx_seq=self.appris_cbs_trx_seq,
            pos=1915, ref="G", alt="GT", ref_len=1, alt_len=2, cds_start_offset=260, cds_stop_offset=1916, cds_len=3315,
            start_index=1654, end_index=1654, ref_codon_dict=ref_codon_dict, alt_codon_dict=alt_codon_dict)

        expected = ("ENST00000398165.7:c.1655_1656inst", "ENST00000398165.7:r.1655_1656insu", "p.Ter552fs",)

        self.assertEqual(expected, observed)

    def test_get_codon_and_aa_changes_trx_not_found(self):
        """Tests that a TranscriptNotFound Exception is raised if transcript is not in the input annotations."""

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
