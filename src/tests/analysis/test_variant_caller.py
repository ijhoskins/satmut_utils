#!/usr/bin/env python3
"""Tests for analysis.variant_caller"""

import collections
import copy
import numpy as np
import pysam
import tempfile
import unittest

import analysis.variant_caller as vc
import analysis.seq_utils as su
import core_utils.file_utils as fu
from satmut_utils.definitions import *

tempfile.tempdir = DEFAULT_TEMPDIR

TEST_SAM = """@HD	VN:1.0	SO:unknown
@SQ	SN:CBS_pEZY3	LN:7108
@RG	ID:CBS1_35_comb_R
MG01HS02:1483:HG7MTBCX3:2:1213:4430:45262_ACTTTTGA	83	CBS_pEZY3	2290	44	105M	=	2290	-105	GGTGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAAGCCTGGGAATGGTGACGCTTGGGATCATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	111<1<11CCGCC<<@1CC</C</0CCG<F<1C<<C/C<<H@G<C<<111<1<1CEC<000<<1HHF<CG<1D1<</0/<<<1/11FFD<</D@D1HFHF?GHIG	AS:i:198	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:2A44T23A33	YS:i:167	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:ACTTTTGA	UG:Z:10000001_R1
MG01HS02:1483:HG7MTBCX3:2:1213:4430:45262_ACTTTTGA	163	CBS_pEZY3	2290	44	105M	=	2290	105	GGAGAAGGGCTTACAAGAGGCGCCCGAGGAGGATGCGGCGGGGGTAATCCTGGGAATGGTGAAGCTTGGGAACATGCTATCGTCCCTGCTTGCCGGGAAAGTGCA	0<00011<0111<1111<1D0/<C/<////0<0<<1/</<FEHH/0<CH11<1DH11DC111<DF1<1D?1D<11<1<<11<<10011111<1</C/C?11=1DD	AS:i:167	XN:i:0	XM:i:10	XO:i:0	XG:i:0	NM:i:10	MD:Z:12C0G1C0C9T2T5A26C15C20G5	YS:i:198	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10000001_R2
MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC	99	CBS_pEZY3	2397	44	64M1D66M	=	2408	163	CGTCAGACCAAGTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACGCGCTGGCAGGCTCTCGCACATCCTGGAGATGGACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTAC	HHHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHIIIHHHIIIIIIIHHIIIIHIIIGIHIIIIIIIIHHIIIIIIIIIIIIIIIIIIHIIIIIIIHIIHIHIIIEHHHIEHIIHIIC	AS:i:245	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:59A4^G66	YS:i:286	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:CGTTGATC	UG:Z:10015877_R1
MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC	147	CBS_pEZY3	2408	44	53M1D98M	=	2397	-163	GTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACGCGCTGGCAGGCTCTCGCACATCCTGGAGATGGACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTACCACAGCACCGGGAAGTCCAGTCAGCGGCAGAT	CIHHHEHGIIIHHFFHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHIHIIHIIIIIIIIIIHIIIHIIIIIIHCIIIIIIIIHIIIIHIIIIIIIHHHIHIGIIIIGIIIIIIIIIIIHHDHGGIIIIIIGGIHIIHIIIIIDCDDD	AS:i:286	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:48A4^G98	YS:i:245	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10015877_R2
MG01HS02:1483:HG7MTBCX3:2:1208:5073:30436_CGTTGATC	99	CBS_pEZY3	2397	44	64M1D66M	=	2409	163	CGTCAGACCAAGTTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACACGCTGGCAGGCTCTCGCACATCCTGGAGATGTACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTAC	HIIIIIIIIIIIIIIIIIIIIIIHHIIHIHIHHIIIIIIIIIIIIIIIIIHHIIIIIIIIIIIIIIIIHIIIIIIIIIIIGCGIHIIIGHIIIIIIIIIIIIIIIIIIIIIIHIIGHHIIHHHHHIHCH@	AS:i:244	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:64^G27G38	YS:i:285	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:CGTTGATC	UG:Z:10015877_R1
MG01HS02:1483:HG7MTBCX3:2:1208:5073:30436_CGTTGATC	147	CBS_pEZY3	2409	44	52M1D98M	=	2397	-163	TTGGCAAAGTCATCTACAAGCAGTTCAAACAGATCCGCCTCACGGACACGCTGGCAGGCTCTCGCACATCCTGGAGATGTACCACTTCGCCCTGGTGGTGCACGAGCAGATCCAGTACCACAGCACCGGGAAGTCCAGTCAGCGGCAGAT	HFCHHIHHIIHIHHIHHGIIIIIIIIIIIIHGIIGIHHHIHEHHFHHHHHIIHIIIIIIIIIIIIIIIHHEGIIHIIHIHIIIIIIIIIIHIIIIIIIHDFIIIIIIIHFHHEHIIIHHHHHIIIIIIIHIIIHHIIIIIGIIIIDDCDD	AS:i:285	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:52^G27G70	YS:i:244	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10015877_R2
"""

# Here we have a mismatch between the read names for a pair
# This is included because ReadMasker custom sort may be affected by format of the qname
TEST_INVALID_SAM = """@HD	VN:1.0	SO:unknown
@SQ	SN:CBS_pEZY3	LN:7108
@RG	ID:CBS1_35_comb_R
MG01HS02:1483:HG7MTBCX3:2:1213:4430:45262_ACTTTTGA	83	CBS_pEZY3	2290	44	105M	=	2290	-105	GGTGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAAGCCTGGGAATGGTGACGCTTGGGATCATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA	111<1<11CCGCC<<@1CC</C</0CCG<F<1C<<C/C<<H@G<C<<111<1<1CEC<000<<1HHF<CG<1D1<</0/<<<1/11FFD<</D@D1HFHF?GHIG	AS:i:198	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:2A44T23A33	YS:i:167	YT:Z:CP	RG:Z:CBS1_35_comb_R	BX:Z:ACTTTTGA	UG:Z:10000001_R1
MG01HS02:1483:HG7MTBCX3:2:1213:4430:45261_ACTTTTGA	163	CBS_pEZY3	2290	44	105M	=	2290	105	GGAGAAGGGCTTACAAGAGGCGCCCGAGGAGGATGCGGCGGGGGTAATCCTGGGAATGGTGAAGCTTGGGAACATGCTATCGTCCCTGCTTGCCGGGAAAGTGCA	0<00011<0111<1111<1D0/<C/<////0<0<<1/</<FEHH/0<CH11<1DH11DC111<DF1<1D?1D<11<1<<11<<10011111<1</C/C?11=1DD	AS:i:167	XN:i:0	XM:i:10	XO:i:0	XG:i:0	NM:i:10	MD:Z:12C0G1C0C9T2T5A26C15C20G5	YS:i:198	YT:Z:CP	RG:Z:CBS1_35_comb_R	UG:Z:10000001_R2
"""


class TestVariantCaller(unittest.TestCase):
    """Tests for VariantCaller."""

    GFF = "CBS_pEZY3.gff"
    GFF_REF = "CBS_pEZY3.fa"
    PRIMERS = "CBS_pEZY3_primers.bed"
    TARGETS = "CBS_pEZY3_targets.bed"

    @classmethod
    def setUpClass(cls):
        """Set up for TestVariantCaller."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.abspath(os.path.join(cls.test_dir, "..", "test_data"))

        cls.ref = os.path.join(cls.test_data_dir, cls.GFF_REF)
        cls.gff = os.path.join(cls.test_data_dir, cls.GFF)
        cls.gff_ref = os.path.join(cls.test_data_dir, cls.GFF_REF)
        cls.primer_bed = os.path.join(cls.test_data_dir, cls.PRIMERS)
        cls.target_bed = os.path.join(cls.test_data_dir, cls.TARGETS)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".test.sam") as test_sam, \
                tempfile.NamedTemporaryFile(mode="w", suffix=".test.sam") as test_invalid_sam:

            test_sam.write(TEST_SAM)
            fu.flush_files((test_sam,))
            cls.test_bam = su.sort_and_index(test_sam.name)

            test_invalid_sam.write(TEST_INVALID_SAM)
            fu.flush_files((test_invalid_sam,))
            cls.test_invalid_bam = su.sort_and_index(test_invalid_sam.name)

        # Select mate-concordant and discordant reads for testing variant calling
        with pysam.AlignmentFile(cls.test_bam, "rb") as test_af:
            for i, read in enumerate(test_af.fetch(until_eof=True)):

                # For some reason, samtools sort is not a stable sort
                if i == 0:
                    cls.test_align_seg_r2_positive_discordant = read
                if i == 1:
                    cls.test_align_seg_r1_negative_discordant = read
                if i == 2:
                    cls.test_align_seg_r1_positive_concordant = read
                if i == 4:
                    cls.test_align_seg_r2_negative_concordant = read

            test_af.reset()

        cls.vc = vc.VariantCaller(
            am=cls.test_bam, ref=cls.ref, trx_gff=cls.gff, gff_ref=cls.gff_ref,
            targets=cls.target_bed, primers=cls.primer_bed, output_dir=cls.tempdir)

        cls.vc_invalid = vc.VariantCaller(
            am=cls.test_invalid_bam, ref=cls.ref, trx_gff=cls.gff, gff_ref=cls.gff_ref,
            targets=cls.target_bed, primers=cls.primer_bed, output_dir=cls.tempdir)

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestVariantCaller."""

        fu.safe_remove((cls.tempdir, cls.test_bam, cls.test_invalid_bam,
                        fu.add_extension(cls.test_bam, su.BAM_INDEX_SUFFIX),
                        fu.add_extension(cls.test_invalid_bam, su.BAM_INDEX_SUFFIX)), force_remove=True)

    def test_is_indel_no_indel(self):
        """Tests that InDels are not detected from a pysam aligned pair."""

        observed = self.vc._is_indel((4, 10, "G"))
        self.assertFalse(observed)

    def test_is_indel_ins(self):
        """Tests that InDels are detected for an insertion."""

        observed = self.vc._is_indel((4, None, None))
        self.assertTrue(observed)

    def test_is_indel_del(self):
        """Tests that InDels are detected for a deletion."""

        observed = self.vc._is_indel((None, 10, "G"))
        self.assertTrue(observed)

    def test_get_haplotype_ref(self):
        """Tests that a proper REF field is constructed from mismatch tuples."""

        expected = ("ACGCTGGGCAG", [0, 10])

        observed = self.vc._get_haplotype_ref(
            contig="CBS_pEZY3", first_pos=2456, last_pos=2466, mm_pos_set={2456, 2466})

        self.assertEqual(expected, observed)

    def test_get_haplotype_dict_within_window(self):
        """Tests that we can generate a dictionary containing positions of mismatches involved in a haplotype."""

        expected_dict_key = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="ACG", alt="GCC", refs=None, alts=None, positions=None)

        expected_dict_val = {2456, 2458}

        # This mismatch exists in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)

        # This mismatch does not exist in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        # but we create it for testing
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2458, "G", "C", 40, 62)

        filt_r_mms = [mm1, mm2]
        observed_haplo_dict, observed_pos_blacklist = self.vc._get_haplotype_dict(filt_r_mms, 3)

        test_1 = list(observed_haplo_dict.keys())[0] == expected_dict_key
        test_2 = len(list(observed_haplo_dict.values())[0] - expected_dict_val) == 0
        test_3 = len(observed_pos_blacklist - expected_dict_val) == 0
        test_res = (test_1, test_2, test_3)
        self.assertTrue(all(test_res))

    def test_get_haplotype_dict_outside_window(self):
        """Tests that we can generate an empty dictionary for mismatches at a distance greater than the window."""

        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2458, "G", "C", 40, 62)
        filt_r_mms = [mm1, mm2]
        observed_haplo_dict, observed_pos_blacklist = self.vc._get_haplotype_dict(filt_r_mms, 2)

        # The dictionary and set should be empty as the haplotype span is greater than the width
        test_1 = len(observed_haplo_dict) == 0
        test_2 = len(observed_pos_blacklist) == 0

        self.assertTrue(all((test_1, test_2)))

    def test_get_haplotype_dict_snp_dint_mnp(self):
        """Tests that we can generate a haplotype dictionary for a di-nt MNP downstream of a SNP."""

        expected_dict_key = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="ACG", alt="GCC", refs=None, alts=None, positions=None)

        expected_dict_val = {2456, 2458}

        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2450, "A", "G", 39, 54)

        # This mismatch exists in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)

        # This mismatch does not exist in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        # but we create it for testing
        mm3 = vc.MM_TUPLE("CBS_pEZY3", 2458, "G", "C", 40, 62)

        filt_r_mms = [mm1, mm2, mm3]
        observed_haplo_dict, observed_pos_blacklist = self.vc._get_haplotype_dict(filt_r_mms, 3)

        test_1 = list(observed_haplo_dict.keys())[0] == expected_dict_key
        test_2 = len(list(observed_haplo_dict.values())[0] - expected_dict_val) == 0
        test_3 = len(observed_pos_blacklist - expected_dict_val) == 0
        test_res = (test_1, test_2, test_3)
        self.assertTrue(all(test_res))

    def test_get_haplotype_dict_snp_trint_mnp(self):
        """Tests that we can generate a haplotype dictionary for a tri-nt MNP downstream of a SNP."""

        expected_dict_key = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="ACG", alt="GAC", refs=None, alts=None, positions=None)

        expected_dict_val = {2456, 2457, 2458}

        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2450, "A", "G", 39, 54)

        # This mismatch exists in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)

        # These mismatches do not exist in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        # but we create it for testing
        mm3 = vc.MM_TUPLE("CBS_pEZY3", 2457, "C", "A", 40, 61)
        mm4 = vc.MM_TUPLE("CBS_pEZY3", 2458, "G", "C", 40, 62)

        filt_r_mms = [mm1, mm2, mm3, mm4]
        observed_haplo_dict, observed_pos_blacklist = self.vc._get_haplotype_dict(filt_r_mms, 3)

        test_1 = list(observed_haplo_dict.keys())[0] == expected_dict_key
        test_2 = len(list(observed_haplo_dict.values())[0] - expected_dict_val) == 0
        test_3 = len(observed_pos_blacklist - expected_dict_val) == 0
        test_res = (test_1, test_2, test_3)
        self.assertTrue(all(test_res))

    def test_get_haplotype_dict_two_dint_mnp(self):
        """Tests that we can generate a haplotype dictionary for two di-nt MNPs not within the window."""

        expected_dict_keys = {
            vc.CALL_TUPLE(contig="CBS_pEZY3", pos=2456, ref="AC", alt="GA", refs=None, alts=None, positions=None),
            vc.CALL_TUPLE(contig="CBS_pEZY3", pos=2461, ref="GG", alt="AA", refs=None, alts=None, positions=None)
        }

        expected_dict_val = {2456, 2457, 2461, 2462}

        # This mismatch exists in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)

        # These mismatches do not exist in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        # but we create it for testing
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2457, "C", "A", 40, 61)
        mm3 = vc.MM_TUPLE("CBS_pEZY3", 2461, "G", "A", 40, 65)
        mm4 = vc.MM_TUPLE("CBS_pEZY3", 2462, "G", "A", 40, 66)

        filt_r_mms = [mm1, mm2, mm3, mm4]
        observed_haplo_dict, observed_pos_blacklist = self.vc._get_haplotype_dict(filt_r_mms, 3)

        test_1 = len(set(observed_haplo_dict.keys()) - expected_dict_keys) == 0
        test_2 = len(list(observed_haplo_dict.values())[0] - expected_dict_val) == 0
        test_3 = len(observed_pos_blacklist - expected_dict_val) == 0
        test_res = (test_1, test_2, test_3)
        self.assertTrue(all(test_res))

    def test_get_haplotype_dict_two_trint_mnp(self):
        """Tests that we can generate a haplotype dictionary for two tri-nt MNPs not within the window."""

        expected_dict_keys = {
            vc.CALL_TUPLE(contig="CBS_pEZY3", pos=2456, ref="ACG", alt="GAC", refs=None, alts=None, positions=None),
            vc.CALL_TUPLE(contig="CBS_pEZY3", pos=2461, ref="GGG", alt="CCC", refs=None, alts=None, positions=None)
        }

        expected_dict_val = {2456, 2457, 2458, 2461, 2462, 2463}

        # This mismatch exists in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)

        # These mismatches do not exist in MG01HS02:1483:HG7MTBCX3:1:1215:11984:60374_CGTTGATC
        # but we create it for testing
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2457, "C", "A", 40, 61)
        mm3 = vc.MM_TUPLE("CBS_pEZY3", 2458, "G", "C", 40, 62)

        mm4 = vc.MM_TUPLE("CBS_pEZY3", 2461, "G", "C", 40, 65)
        mm5 = vc.MM_TUPLE("CBS_pEZY3", 2462, "G", "C", 40, 66)
        mm6 = vc.MM_TUPLE("CBS_pEZY3", 2463, "G", "C", 40, 67)

        filt_r_mms = [mm1, mm2, mm3, mm4, mm5, mm6]
        observed_haplo_dict, observed_pos_blacklist = self.vc._get_haplotype_dict(filt_r_mms, 3)

        test_1 = len(set(observed_haplo_dict.keys()) - expected_dict_keys) == 0
        test_2 = len(list(observed_haplo_dict.values())[0] - expected_dict_val) == 0
        test_3 = len(observed_pos_blacklist - expected_dict_val) == 0
        test_res = (test_1, test_2, test_3)
        self.assertTrue(all(test_res))

    def test_call_haplotypes_one_mismatch(self):
        """Tests that no haplotypes are returned if there is only one mismatch in the mate pair."""

        filt_r_mms = [vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)]

        observed = self.vc._call_haplotypes(filt_r_mms)

        self.assertIsNone(observed)

    def test_call_haplotypes_two_mismatches_outside_window(self):
        """Tests that no haplotypes are returned two mismatches are at distance greater than the window."""

        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2466, "G", "C", 40, 70)
        filt_r_mms = [mm1, mm2]

        observed = self.vc._call_haplotypes(filt_r_mms)

        self.assertIsNone(observed)

    # No need to test the other logic path through _call_haplotypes() as we test _get_haplotype_dict()

    def test_call_snps_no_haplotypes(self):
        """Tests that we call SNPs for mismatches when there are no called MNPs/haplotypes."""

        expected_dict_keys = {
            vc.CALL_TUPLE(contig="CBS_pEZY3", pos=2456, ref="A", alt="G", refs=None, alts=None, positions=None),
            vc.CALL_TUPLE(contig="CBS_pEZY3", pos=2466, ref="G", alt="C", refs=None, alts=None, positions=None)
        }

        expected_dict_vals = {2456, 2466}

        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2466, "G", "C", 40, 70)
        filt_r_mms = [mm1, mm2]
        haplo_dict, pos_blacklist = self.vc._get_haplotype_dict(filt_r_mms, 3)

        # Because MNP span is less than the distance between mismatches, both mismatches should be called as SNPs
        observed = self.vc._call_snps(filt_r_mms, pos_blacklist, haplo_dict)

        test_1 = len(expected_dict_keys - set(observed.keys())) == 0

        # Just combine all positions from the dict values for testing
        test_2_res = set()
        for e in observed.values():
            test_2_res |= e

        test_2 = len(expected_dict_vals - test_2_res) == 0
        test_res = (test_1, test_2)
        self.assertTrue(all(test_res))

    def test_call_snps_mismatch_outside_haplotypes(self):
        """Tests that we call SNPs for mismatches not in an already-called MNP/haplotype."""

        expected_dict_key = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2470, ref="C", alt="T", refs=None, alts=None, positions=None)

        expected_dict_val = {2470}

        # mm1 and mm2 will be a part of a haplotype and mm3 will be the SNP
        # even the mm3 is within the window distance to mm2, because we join from left to right, it will
        # be excluded by the pos_blacklist
        mm1 = vc.MM_TUPLE("CBS_pEZY3", 2456, "A", "G", 39, 60)
        mm2 = vc.MM_TUPLE("CBS_pEZY3", 2466, "G", "C", 40, 70)
        mm3 = vc.MM_TUPLE("CBS_pEZY3", 2470, "C", "T", 40, 74)
        filt_r_mms = [mm1, mm2, mm3]
        haplo_dict, pos_blacklist = self.vc._get_haplotype_dict(filt_r_mms, 11)

        # Because MNP span is greater than the distance between mismatches, both mismatches should be called as SNPs
        observed = self.vc._call_snps(filt_r_mms, pos_blacklist, haplo_dict)

        test_1 = expected_dict_key == list(observed.keys())[0]
        test_2 = expected_dict_val == list(observed.values())[0]
        test_res = (test_1, test_2)
        self.assertTrue(all(test_res))

    def test_enumerate_mismatches(self):
        """Tests that we properly enumerate mismatches in a read."""

        # Mismatch tuples have a 1-based position and the read_pos is 1-based with respect to the 5' end of the read
        expected = [vc.MM_TUPLE(contig="CBS_pEZY3", pos=2292, ref="A", alt="T", bq=16, read_pos=103),
                    vc.MM_TUPLE(contig="CBS_pEZY3", pos=2337, ref="T", alt="G", bq=16, read_pos=58),
                    vc.MM_TUPLE(contig="CBS_pEZY3", pos=2361, ref="A", alt="T", bq=16, read_pos=34)]

        observed = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r1_negative_discordant, min_bq=10)

        self.assertEqual(expected, observed)

    def test_enumerate_mismatches_bq_filter(self):
        """Tests that we do not enumerate mismatches in a read if the BQs are lower than the min BQ threshold."""

        observed = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r1_negative_discordant, min_bq=20)
        self.assertEqual(0, len(observed))

    def test_intersect_mismatches_discordant(self):
        """Test that we find no intersected mismatches for a discordant pair."""

        expected = ([], [])

        # Here we have R1 mismatches but no R2 mismatches
        r1_mms = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r1_negative_discordant, min_bq=10)
        r2_mms = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r2_positive_discordant, min_bq=10)

        observed = self.vc._intersect_mismatches(r1_mms, r2_mms)

        self.assertEqual(expected, observed)

    def test_intersect_mismatches_concordant(self):
        """Test that we find intersected mismatches for a concordant pair."""

        expected = ([vc.MM_TUPLE(contig="CBS_pEZY3", pos=2456, ref="A", alt="G", bq=39, read_pos=60)],
                    [vc.MM_TUPLE(contig="CBS_pEZY3", pos=2456, ref="A", alt="G", bq=40, read_pos=103)])

        r1_mms = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r1_positive_concordant, min_bq=30)
        r2_mms = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r2_negative_concordant, min_bq=30)
        observed = self.vc._intersect_mismatches(r1_mms, r2_mms)

        self.assertEqual(expected, observed)

    def test_unpack_stats_snp(self):
        """Tests extraction of positions, reference bases, alternate bases and quality info from SNPs in a pair."""

        r1_mms = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r1_positive_concordant, min_bq=30)
        r2_mms = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r2_negative_concordant, min_bq=30)
        filt_r1_mms, filt_r2_mms = self.vc._intersect_mismatches(r1_mms, r2_mms)
        support_pos = {2456}

        refs, alts, positions, per_bp_stats = self.vc._unpack_stats(
            supporting_positions=support_pos, filt_r1_mms=filt_r1_mms, filt_r2_mms=filt_r2_mms)

        test_1 = refs == "A"
        test_2 = alts == "G"
        test_3 = positions == "2456"
        test_4 = per_bp_stats == vc.PER_BP_STATS(
            r1_bqs=np.array([39], dtype=np.int32), r2_bqs=np.array([40], dtype=np.int32),
            r1_read_pos=np.array([60], dtype=np.int32), r2_read_pos=np.array([103], dtype=np.int32))

        test_res = (test_1, test_2, test_3, test_4)
        self.assertTrue(all(test_res))

    def test_assign_stats_snp(self):
        """Tests assignment of SNP BQ and read position stats to the variant call dict."""

        call_tuple = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="A", alt="G", refs="A", alts="G", positions="2456")

        self.vc.variant_counts = collections.OrderedDict()
        self.vc.variant_counts[call_tuple] = [
            [0, [], [], []],  # R1, +
            [0, [], [], []],  # R1, -
            [0, [], [], []],  # R2, +
            [0, [], [], []]  # R2, -
        ]

        per_bp_stats = vc.PER_BP_STATS(
            r1_bqs=np.array([39], dtype=np.int32), r2_bqs=np.array([40], dtype=np.int32),
            r1_read_pos=np.array([60], dtype=np.int32), r2_read_pos=np.array([103], dtype=np.int32))

        expected = collections.OrderedDict()
        expected[call_tuple] = [
                [0, ["39"], ["60"], []],  # R1, +
                [0, [], [], []],  # R1, -
                [0, [], [], []],  # R2, +
                [0, ["40"], ["103"], []]   # R2, -
            ]

        # This method has a side-effect on the caller obj dict
        self.vc._assign_stats(r_index=self.vc.R1_PLUS_INDEX, call_tuple=call_tuple, per_bp_stats=per_bp_stats)
        self.vc._assign_stats(r_index=self.vc.R2_MINUS_INDEX, call_tuple=call_tuple, per_bp_stats=per_bp_stats)

        self.assertEqual(expected, self.vc.variant_counts)

    def test_assign_stats_mnp(self):
        """Tests assignment of MNP BQ and read position stats to the variant call dict."""

        call_tuple = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="AC", alt="GT", refs="A,C", alts="G,T", positions="2456,2457")

        self.vc.variant_counts = collections.OrderedDict()
        self.vc.variant_counts[call_tuple] = [
                [0, [], [], []],  # R1, +
                [0, [], [], []],  # R1, -
                [0, [], [], []],  # R2, +
                [0, [], [], []]   # R2, -
            ]

        per_bp_stats = vc.PER_BP_STATS(
            r1_bqs=np.array([39, 39], dtype=np.int32), r2_bqs=np.array([40, 39], dtype=np.int32),
            r1_read_pos=np.array([60, 61], dtype=np.int32), r2_read_pos=np.array([103, 102], dtype=np.int32))

        expected = collections.OrderedDict()
        expected[call_tuple] = [
                [0, ["39,39"], ["60,61"], []],  # R1, +
                [0, [], [], []],  # R1, -
                [0, [], [], []],  # R2, +
                [0, ["40,39"], ["103,102"], []]   # R2, -
            ]

        # This method has a side-effect on the caller obj dict
        self.vc._assign_stats(r_index=self.vc.R1_PLUS_INDEX, call_tuple=call_tuple, per_bp_stats=per_bp_stats)
        self.vc._assign_stats(r_index=self.vc.R2_MINUS_INDEX, call_tuple=call_tuple, per_bp_stats=per_bp_stats)

        self.assertEqual(expected, self.vc.variant_counts)

    def test_add_counts_and_stats_mnp(self):
        """Tests that counts, BQs, NMs, and read positions for an MNP are added to the variant counts dict."""

        self.vc.variant_counts = collections.OrderedDict()

        per_bp_stats = vc.PER_BP_STATS(
            r1_bqs=np.array([39, 39], dtype=np.int32), r2_bqs=np.array([40, 39], dtype=np.int32),
            r1_read_pos=np.array([60, 61], dtype=np.int32), r2_read_pos=np.array([103, 102], dtype=np.int32))

        call_tuple = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="AC", alt="GT", refs="A,C", alts="G,T", positions="2456,2457")

        expected = collections.OrderedDict()
        expected[call_tuple] = [
            [1, ["39,39"], ["60,61"], [2]],  # R1, +
            [0, [], [], []],  # R1, -
            [0, [], [], []],  # R2, +
            [1, ["40,39"], ["103,102"], [2]]  # R2, -
        ]

        self.vc._add_counts_and_stats(
            call_tuple=call_tuple, per_bp_stats=per_bp_stats, r1_nm=2, r2_nm=2,
            r1_strand=su.Strand.PLUS, r2_strand=su.Strand.MINUS)

        self.assertEqual(expected, self.vc.variant_counts)

    def test_update_counts(self):
        """Tests that the variant counts dict is updated starting from an initial pair-specific call dict."""

        self.vc.variant_counts = collections.OrderedDict()

        expected = collections.OrderedDict()

        call_tuple_final = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="A", alt="G", refs="A", alts="G", positions="2456")

        expected[call_tuple_final] = [
            [1, ["39"], ["60"], [2]],  # R1, +
            [0, [], [], []],  # R1, -
            [0, [], [], []],  # R2, +
            [1, ["40"], ["103"], [2]]  # R2, -
        ]

        collective_variants = collections.defaultdict(set)

        call_tuple_initial = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="A", alt="G", refs=None, alts=None, positions=None)

        collective_variants[call_tuple_initial] |= {2456}

        r1_mms = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r1_positive_concordant, min_bq=30)
        r2_mms = self.vc._enumerate_mismatches(align_seg=self.test_align_seg_r2_negative_concordant, min_bq=30)
        filt_r1_mms, filt_r2_mms = self.vc._intersect_mismatches(r1_mms, r2_mms)

        self.vc._update_counts(
            collective_variants=collective_variants, filt_r1_mms=filt_r1_mms, filt_r2_mms=filt_r2_mms,
            r1_nm=2, r2_nm=2, r1_strand=su.Strand.PLUS, r2_strand=su.Strand.MINUS)

        self.assertEqual(expected, self.vc.variant_counts)

    def test_get_unmasked_positions(self):
        """Tests that unmasked (non-zero) BQs indices are returned."""

        expected = set(self.test_align_seg_r1_positive_concordant.get_reference_positions()[15:])

        test_align_seg = copy.deepcopy(self.test_align_seg_r1_positive_concordant)
        test_align_seg_bqs = list(test_align_seg.query_qualities)
        test_align_seg_bqs[0:15] = [0] * 15
        test_align_seg.query_qualities = test_align_seg_bqs

        observed = self.vc._get_unmasked_positions(test_align_seg)

        self.assertEqual(0, len(expected - observed))

    def test_reads_overlap(self):
        """Tests that overlapping reads are detected as such."""

        observed = self.vc._reads_overlap(
            r1=self.test_align_seg_r1_negative_discordant, r2=self.test_align_seg_r2_positive_discordant)

        self.assertTrue(observed)

    def test_reads_overlap_none(self):
        """Tests that non-overlapping reads are detected as such."""

        observed = self.vc._reads_overlap(
            r1=self.test_align_seg_r1_negative_discordant, r2=self.test_align_seg_r2_negative_concordant)

        self.assertFalse(observed)

    def test_update_pos_dp_partial_masking(self):
        """Tests that the coverage DP dict is updated for reads with partial BQ masking."""

        # non-zero coverage will start at 2412 and end at 2559, 1-based
        expected = collections.defaultdict(int)
        for pos in range(2411, 2559):
            expected[vc.COORDINATE_KEY("CBS_pEZY3", pos + 1)] += 1

        # Make sure to re-init the dict as it us updated by side effect by various test methods
        self.vc.coordinate_counts = collections.defaultdict(int)

        # Include masked BQs which should drop out of the fragment coverage
        test_r1_align_seg = copy.deepcopy(self.test_align_seg_r1_positive_concordant)
        test_r1_align_seg_bqs = list(test_r1_align_seg.query_qualities)
        test_r1_align_seg_bqs[0:15] = [0] * 15
        test_r1_align_seg.query_qualities = test_r1_align_seg_bqs

        test_r2_align_seg = copy.deepcopy(self.test_align_seg_r2_negative_concordant)
        test_r2_align_seg_bqs = list(test_r2_align_seg.query_qualities)
        test_r2_align_seg_bqs[0:4] = [0] * 4
        test_r2_align_seg.query_qualities = test_r2_align_seg_bqs

        self.vc._update_pos_dp(r1=test_r1_align_seg, r2=test_r2_align_seg)

        # Test that the keys are in the expected range
        test_1 = len(set(expected.keys()) - set(self.vc.coordinate_counts.keys())) == 0

        # Test that 1 is set for all positions
        test_2 = len(set(self.vc.coordinate_counts.values()) - {1}) == 0

        test_res = (test_1, test_2)
        self.assertTrue(all(test_res))

    def test_update_pos_dp_all_masked(self):
        """Tests that None is returned for reads with complete BQ masking."""

        # Make sure to re-init the dict as it us updated by side effect by various test methods
        self.coordinate_counts = collections.defaultdict(int)

        test_r1_align_seg = copy.deepcopy(self.test_align_seg_r1_positive_concordant)
        test_r1_align_seg_bqs = [0] * test_r1_align_seg.query_length
        test_r1_align_seg.query_qualities = test_r1_align_seg_bqs

        test_r2_align_seg = copy.deepcopy(self.test_align_seg_r2_negative_concordant)
        test_r2_align_seg_bqs = [0] * test_r2_align_seg.query_length
        test_r2_align_seg.query_qualities = test_r2_align_seg_bqs

        observed = self.vc._update_pos_dp(r1=test_r1_align_seg, r2=self.test_align_seg_r2_negative_concordant)

        self.assertIsNone(observed)

    def test_iterate_over_reads_concordant(self):
        """Tests that only concordant variants are called."""

        expected = collections.OrderedDict()

        # There are only two concordant variants in the test data
        call_tuple_1 = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="A", alt="G", refs="A", alts="G", positions="2456")

        call_tuple_2 = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2489, ref="G", alt="T", refs="G", alts="T", positions="2489")

        expected[call_tuple_1] = [
            [1, ["39"], ["60"], [2]],  # R1, +
            [0, [], [], []],  # R1, -
            [0, [], [], []],  # R2, +
            [1, ["40"], ["103"], [2]]  # R2, -
        ]

        expected[call_tuple_2] = [
            [1, ["40"], ["92"], [2]],  # R1, +
            [0, [], [], []],  # R1, -
            [0, [], [], []],  # R2, +
            [1, ["39"], ["71"], [2]]  # R2, -
        ]

        self.vc.variant_counts = collections.OrderedDict()

        with pysam.AlignmentFile(self.vc.vc_preprocessor.r1_calling_bam, "rb", check_sq=False) as af1, \
                pysam.AlignmentFile(self.vc.vc_preprocessor.r2_calling_bam, "rb", check_sq=False) as af2:

            self.vc._iterate_over_reads(af1=af1, af2=af2, min_bq=30, max_nm=20, max_mnp_window=3)

        self.assertEqual(expected, self.vc.variant_counts)

    def test_iterate_over_reads_nm_filter(self):
        """Tests that the edit distance filter is applied."""

        self.vc.variant_counts = collections.OrderedDict()

        with pysam.AlignmentFile(self.vc.vc_preprocessor.r1_calling_bam, "rb", check_sq=False) as af1, \
                pysam.AlignmentFile(self.vc.vc_preprocessor.r2_calling_bam, "rb", check_sq=False) as af2:

            self.vc._iterate_over_reads(af1=af1, af2=af2, min_bq=30, max_nm=1, max_mnp_window=3)

        self.assertEqual(0, len(self.vc.variant_counts))

    def test_iterate_over_reads_qname_mismatch(self):
        """Tests that an exception is generated if we have a mismatch between qnames of a pair."""

        self.vc.variant_counts = collections.OrderedDict()

        with pysam.AlignmentFile(self.vc_invalid.vc_preprocessor.r1_calling_bam, "rb", check_sq=False) as af1, \
                pysam.AlignmentFile(self.vc_invalid.vc_preprocessor.r2_calling_bam, "rb", check_sq=False) as af2:

            with self.assertRaises(RuntimeError):
                self.vc._iterate_over_reads(af1=af1, af2=af2, min_bq=30, max_nm=1, max_mnp_window=3)

    def test_summarize_stats_bq(self):
        """Tests summarization of BQ stats for an MNP."""

        expected_1 = "38.0"
        expected_2 = "39.0"

        # An di-nt MNP with two supporting reads; normally we would have R2, - data but no need for testing
        var_list = [
            [2, ["39,39", "37,39"], ["60,61", "59,60"], [2, 2]],  # R1, +
            [0, [], [], []],  # R1, -
            [0, [], [], []],  # R2, +
            [0, [], [], []]  # R2, -
        ]

        observed_1 = self.vc._summarize_stats(
            var_list=var_list, read_index=self.vc.R1_PLUS_INDEX, stat_index=self.vc.R_BQ_INDEX, pos_index=0)

        observed_2 = self.vc._summarize_stats(
            var_list=var_list, read_index=self.vc.R1_PLUS_INDEX, stat_index=self.vc.R_BQ_INDEX, pos_index=1)

        test_1 = expected_1 == observed_1
        test_2 = expected_2 == observed_2
        test_res = (test_1, test_2)
        self.assertTrue(all(test_res))

    def test_summarize_stats_nm(self):
        """Tests summarization of read position stats for an MNP."""

        expected_1 = "59.5"
        expected_2 = "60.5"

        # An di-nt MNP with two supporting reads; normally we would have R2, - data but no need for testing
        var_list = [
            [2, ["39,39", "37,39"], ["60,61", "59,60"], [2, 2]],  # R1, +
            [0, [], [], []],  # R1, -
            [0, [], [], []],  # R2, +
            [0, [], [], []]  # R2, -
        ]

        observed_1 = self.vc._summarize_stats(
            var_list=var_list, read_index=self.vc.R1_PLUS_INDEX, stat_index=self.vc.R_RP_INDEX, pos_index=0)

        observed_2 = self.vc._summarize_stats(
            var_list=var_list, read_index=self.vc.R1_PLUS_INDEX, stat_index=self.vc.R_RP_INDEX, pos_index=1)

        test_1 = expected_1 == observed_1
        test_2 = expected_2 == observed_2
        test_res = (test_1, test_2)
        self.assertTrue(all(test_res))

    # Skip testing of _get_read_pos_stats(), _get_read_bq_stats() as they simply wrap _summarize_stats()

    def test_get_read_nm_stats(self):
        """Tests that NM stats are summarized when supporting reads are present."""

        expected = ("4.0", su.R_COMPAT_NA, su.R_COMPAT_NA, su.R_COMPAT_NA)

        var_list = [
            [2, ["39,39", "37,39"], ["60,61", "59,60"], [2, 6]],  # R1, +
            [0, [], [], []],  # R1, -
            [0, [], [], []],  # R2, +
            [0, [], [], []]  # R2, -
        ]

        observed = self.vc._get_read_nm_stats(var_list=var_list)

        self.assertEqual(expected, observed)

    def test_call_variants_exceeds_min_qnames(self):
        """Tests that variants are called and stats are reported for each component base in a call."""

        expected = collections.OrderedDict()

        # This is quite verbose but we need to check that we can assign stats for each component base in a MNP
        vcst_1 = vc.VARIANT_CALL_SUMMARY_TUPLE(
            POS_NT=2456, REF_NT="A", ALT_NT="G", UP_REF_NT="C", DOWN_REF_NT="C", DP=3, CAO=3,
            CAF=1.0, NORM_CAO=1000000.0, R1_PLUS_AO=2, R1_MINUS_AO=1, R2_PLUS_AO=0, R2_MINUS_AO=0,
            R1_PLUS_MED_RP="59.5", R1_MINUS_MED_RP="60", R2_PLUS_MED_RP=su.R_COMPAT_NA, R2_MINUS_MED_RP=su.R_COMPAT_NA,
            R1_PLUS_MED_BQ="37.5", R1_MINUS_MED_BQ="20", R2_PLUS_MED_BQ=su.R_COMPAT_NA, R2_MINUS_MED_BQ=su.R_COMPAT_NA,
            R1_PLUS_MED_NM="3.5", R1_MINUS_MED_NM="3", R2_PLUS_MED_NM=su.R_COMPAT_NA, R2_MINUS_MED_NM=su.R_COMPAT_NA
        )

        vcst_2 = vc.VARIANT_CALL_SUMMARY_TUPLE(
            POS_NT=2457, REF_NT="C", ALT_NT="T", UP_REF_NT="A", DOWN_REF_NT="G", DP=3, CAO=3,
            CAF=1.0, NORM_CAO=1000000.0, R1_PLUS_AO=2, R1_MINUS_AO=1, R2_PLUS_AO=0, R2_MINUS_AO=0,
            R1_PLUS_MED_RP="60.5", R1_MINUS_MED_RP="61", R2_PLUS_MED_RP=su.R_COMPAT_NA, R2_MINUS_MED_RP=su.R_COMPAT_NA,
            R1_PLUS_MED_BQ="39.0", R1_MINUS_MED_BQ="30", R2_PLUS_MED_BQ=su.R_COMPAT_NA, R2_MINUS_MED_BQ=su.R_COMPAT_NA,
            R1_PLUS_MED_NM="3.5", R1_MINUS_MED_NM="3", R2_PLUS_MED_NM=su.R_COMPAT_NA, R2_MINUS_MED_NM=su.R_COMPAT_NA
        )

        expected[vc.VARIANT_CALL_KEY_TUPLE(contig="CBS_pEZY3", pos=2456, ref="AC", alt="GT", index=0)] = vcst_1
        expected[vc.VARIANT_CALL_KEY_TUPLE(contig="CBS_pEZY3", pos=2456, ref="AC", alt="GT", index=1)] = vcst_2

        self.vc.variant_counts = collections.OrderedDict()

        call_tuple = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="AC", alt="GT", refs="A,C", alts="G,T", positions="2456,2457")

        # Normally R2 data would be present but for the purpose of testing omit
        var_list = [
            [2, ["38,39", "37,39"], ["60,61", "59,60"], [3, 4]],  # R1, +
            [1, ["20,30"], ["60,61"], [3]],  # R1, -
            [0, [], [], []],  # R2, +
            [0, [], [], []]  # R2, -
        ]

        self.vc.variant_counts[call_tuple] = var_list
        self.vc.coordinate_counts = collections.defaultdict(int)
        self.vc.coordinate_counts[vc.COORDINATE_KEY("CBS_pEZY3", 2456)] += 3
        self.vc.coordinate_counts[vc.COORDINATE_KEY("CBS_pEZY3", 2457)] += 3

        observed = self.vc._call_variants(min_supporting_qnames=1)

        self.assertEqual(expected, observed)

    def test_call_variants_min_qname_filter(self):
        """Tests that no variants are called if the CAO does not exceed min_supporting_qnames."""

        self.vc.variant_counts = collections.OrderedDict()

        call_tuple = vc.CALL_TUPLE(
            contig="CBS_pEZY3", pos=2456, ref="AC", alt="GT", refs="A,C", alts="G,T", positions="2456,2457")

        # Normally R2 data would be present but for the purpose of testing omit
        var_list = [
            [2, ["38,39", "37,39"], ["60,61", "59,60"], [3, 4]],  # R1, +
            [1, ["20,30"], ["60,61"], [3]],  # R1, -
            [0, [], [], []],  # R2, +
            [0, [], [], []]  # R2, -
        ]

        self.vc.variant_counts[call_tuple] = var_list
        self.coordinate_counts = collections.defaultdict(int)
        self.coordinate_counts[vc.COORDINATE_KEY("CBS_pEZY3", 2456)] += 3
        self.coordinate_counts[vc.COORDINATE_KEY("CBS_pEZY3", 2457)] += 3

        observed = self.vc._call_variants(min_supporting_qnames=4)

        self.assertEqual(0, len(observed))

    def test_workflow(self):
        """Tests that output files are generated from the whole workflow."""

        vcf, bed = self.vc.workflow(min_bq=30, max_nm=5, min_supporting_qnames=1, max_mnp_window=3,
                                    out_prefix=os.path.join(self.tempdir, "test"))
        self.assertTrue(all((os.path.exists(vcf), os.path.exists(bed),)))

    def test_workflow_primers_and_min_bq(self):
        """Tests that a NotImplementedError is raised if primers are provided and min_bq is 0."""

        with self.assertRaises(NotImplementedError):
            _, _ = self.vc.workflow(min_bq=0, max_nm=5, min_supporting_qnames=1, max_mnp_window=3,
                                    out_prefix=os.path.join(self.tempdir, "test"))

    def test_workflow_max_mnp_window(self):
        """Tests that a NotImplementedError is raised if max_mnp_window > 3."""

        with self.assertRaises(NotImplementedError):
            _, _ = self.vc.workflow(min_bq=30, max_nm=5, min_supporting_qnames=1, max_mnp_window=4,
                                    out_prefix=os.path.join(self.tempdir, "test"))
