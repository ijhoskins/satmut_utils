#!/usr/bin/env python3
"""Tests for analysis.seq_utils."""

import os
import pysam
import random
from shutil import copyfile
import tempfile
import unittest

import analysis.aligners as al
from analysis.references import index_reference
import analysis.seq_utils as su
import core_utils.file_utils as fu
from satmut_utils.definitions import DEFAULT_QUALITY_OFFSET

TEST_FASTA = """>gi|571026644|ref|NM_014453.3|:1-941 Homo sapiens charged multivesicular body protein 2A (CHMP2A), transcript variant 2, mRNA
ATGTTACGTCACGTCGTTTAATCCGGAAACGGCGGCGGCGGCGACAGGACCGAGGGGCCTTAGTTGGTGG
GCAAGTCGGGGATCCCAGAAAGAGAAGCGTGACCCGGAAGCGGAAACGGGTGTCCGTCCCAGCTCCGGCC
TGCCAGTGAGCTTCTACCATCATGGACCTATTGTTCGGGCGCCGGAAGACGCCAGAGGAGCTACTGCGGC
AGAACCAGAGGGCCCTGAACCGTGCCATGCGGGAGCTGGACCGCGAGCGACAGAAACTAGAGACCCAGGA
GAAGAAAATCATTGCAGACATTAAGAAGATGGCCAAGCAAGGCCAGATGGATGCTGTTCGCATCATGGCA
AAAGACTTGGTGCGCACCCGGCGCTATGTGCGCAAGTTTGTATTGATGCGGGCCAACATCCAGGCTGTGT
CCCTCAAGATCCAGACACTCAAGTCCAACAACTCGATGGCACAAGCCATGAAGGGTGTCACCAAGGCCAT
GGGCACCATGAACAGACAGCTGAAGTTGCCCCAGATCCAGAAGATCATGATGGAGTTTGAGCGGCAGGCA
GAGATCATGGATATGAAGGAGGAGATGATGAATGATGCCATTGATGATGCCATGGGTGATGAGGAAGATG
AAGAGGAGAGTGATGCTGTGGTGTCCCAGGTTCTGGATGAGCTGGGACTTAGCCTAACAGATGAGCTGTC
GAACCTCCCCTCAACTGGGGGCTCGCTTAGTGTGGCTGCTGGTGGGAAAAAAGCAGAGGCCGCAGCCTCA
GCCCTAGCTGATGCTGATGCAGACCTGGAGGAACGGCTTAAGAACCTGCGGAGGGACTGAGTGCCCCTGC
CACTCCGAGATAACCAGTGGATGCCCAGGATCTTTTACCACAACCCCTCTGTAATAAAAGAGATTTGACA
CTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""
TEST_FASTA_LEN = 941

TEST_BED = """chr19	59066354	59066491	gi|571026644|ref|NM_014453.3|:1-137	0	-
chr19	59065411	59065603	gi|571026644|ref|NM_014453.3|:138-329	100	-
chr19	59063625	59063805	gi|571026644|ref|NM_014453.3|:330-509	200	-
chr19	59063421	59063552	gi|571026644|ref|NM_014453.3|:510-640	300	-
chr19	59063272	59063334	gi|571026644|ref|NM_014453.3|:641-702	400	-
chr19	59062932	59063143	gi|571026644|ref|NM_014453.3|:703-913	500	-
"""

TEST_SAM = """@HD	VN:1.0	SO:coordinate
@SQ	SN:chr19	LN:59128983
@RG	ID:20Feb2019.08:11PM.read_gen.dna
chr19:59063248-59063359_657LTC3T6CFQ	163	chr19	59063256	44	103M1S	=	59063256	-104	GGATTTGGAGCACTCACTCGACAGCTCATCTGTTAGGCTAAGTCCCAGCTCATCCAGAACCTGGGACACCACAGCATCACTAGGGAAGAGAGAGAACTCAGTGT	IRLRQKLQSSMQKIJOKLOSKSPQISJKQSJMJJIOPMKQISLSPIKMPSQOMLQMMIORJQQLOPLIOSIQIOPNSQKMMSSJQNMPOSRIMMKJRRSOLOOL	AS:i:206	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:103	YS:i:200	YT:Z:CP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59063248-59063359_657LTC3T6CFQ	83	chr19	59063256	44	104M	=	59063256	-104	GGATTTGGAGCACTCACTCGACAGCTCATCTGTTAGGCTAAGTCCCAGCTCTTCCAGAACCTGGGACACCACAGCATCACTAGGGAAGAGAGAGAACTCAGTGA	IRKJKORRNKRMRIMQSKINRPJMOLMSLROOQSMJKNSQKSIOLQKJNMJOOPLPKKLMSPMOMQROQJKKRPSPLNRLNMNSRMOSNNMQJQJOMSIRLKPR	AS:i:200	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:51A52	YS:i:206	YT:Z:CP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59063248-59063359_KOQWXDQ80WB2	99	chr19	59063268	44	77M	=	59063268	-77	CTCACTCGACAGCTCATCGGTTAGGCTAAGTCCCAGCTCATCCAGAACCTGGGACACCACAGCATCACTAGGGAAGA	MILPPRLKINQJPQLMRLRMJOPLNPPJSKPPIIRJRSRJKOIQQMRLRKKSPNKMNIJPLSSRRMLOQLORSNKMN	AS:i:146	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:18T58	YS:i:146	YT:Z:CP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59063248-59063359_KOQWXDQ80WB2	147	chr19	59063268	44	77M	=	59063268	-77	CTCACTCGACAGCTCATCTGTTAGGCTAAGTCCCAGCTCATTCAGAACCTGGGACACCACAGCATCACTAGGGAAGA	QOMQLLLLIMPMINRISPRISPIJJIJMNKRIRQOQSJKQSJORSJLKSIPKKNNRLLLSSRSPMQNKQKNMMRNLJ	AS:i:146	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:41C35	YS:i:146	YT:Z:CP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59063397-59063577_W8X67QWSMZTL	145	chr19	59063413	44	142M	=	59063414	-142	TCCCCATACCTCTCCTCTTCATCTTCCTCATCACCCATGGCATCATCAATGGCATCATTCATCATCTCCTCCTTCATATCCATGATCTCTGCCTGCCGCTCAAACTCCATCATTATCTTCTGGATCTGGGGCAACTTCAGCT	LSKJMKRJKNQILPMPNLNMNLQMNQLNRRSRKLQJLNLIJRNMPMJIIOOOMPLRKNMPOSMMNPLNSSIJOJLMNIJRNSMOMOIPNLOPORSQMNMMRMOLKRLNJLPRNRLRIIRPPNPPJJLJKKIPMPISOIPIRP	AS:i:276	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:113G28	YS:i:282	YT:Z:DP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59063397-59063577_W8X67QWSMZTL	97	chr19	59063414	44	1S141M	=	59063413	-142	ACCCCATACCTCTCCTCTTCATCTTCCTCATCACCCATGGCATCATCAATGGCATCATTCATCATCTCCTCCTTCATATCCATGATCTCTGCCTGCCGCTCAAACTCCATCATGATCTTCTGGATCTGGGGCAACTTCAGCT	LNLNSKNPLLLQJRINNNMNPLIKNSSISJNMJOKLPPRKILRNMRKMIRJOOLINSQIOSNJMJSLMPRKQLNQNSMSLLOSLOKNKIKNPQKQMOJMLIPROKNPNKKOJJLKLQOOOJLNORRKPIQJPKPLRJQKQOR	AS:i:282	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:141	YS:i:276	YT:Z:DP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59063397-59063577_GDBKER4FJ023	163	chr19	59063414	44	150M	=	59063415	151	CCCCATACCTCTCCTTTTCATCTTCCTCATCACCCATGGCATCATCAATGGCATCATTCATCATCTCCTCCTTCATATCCATGATCTCTGCCTGCCGCTCAAACTCCATCATGATCTTCTGGATCTGGGGCAACTTCAGCTGGACGTGAC	QKSQLOPKQNLKKIPKJLLJLRRSQRSINMPLINIISIILOMOKJLLPLLIPLMSSPNLSRSJLLKNIOQJRNQOKOIKKKNPPSPIPQRKILOQPMLOSSSKJMINSOSKJSNMPQOPJJNMJKQQOSNPIIOMQLLMJKSJIOILRON	AS:i:284	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:15C128A5	YS:i:284	YT:Z:CP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59063397-59063577_GDBKER4FJ023	83	chr19	59063415	44	150M	=	59063414	-151	CCCATACCTCTCCTCTTCATCTTCCTCATCACCCATGGCATCATCAATGGCATCATTCATCCTCTCCTCCTTCATATCCATGATCTCTGCCTGCCGCTCAAACTCCATCATGATCTTCTTGATCTGGGGCAACTTCAGCTGGAAGTGACA	MOIRKKSJINJOQMIQNIRKNIMNMPSPOSLQQIINMLJQKKIKMNNKRJRLNMMKPRNPQOMQJSLSOQNPNIKKKIJPRMIOJMQPPJIPNMKRKQOOJNNLRIOQSLLMIRMQPPOLKLQNPPKJNRPMQJRNQKQPQJMSLKPQNP	AS:i:284	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:61A57G30	YS:i:284	YT:Z:CP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59066330-59066516_THRZC67832XD	99	chr19	59066332	44	150M	=	59066343	161	ACGGTCCAGCAGGGGCCACTCACCGGAGCTGGGACGGACACCCGTTTCCGCTTCCGGGTCACGCTTATCTTTCTGGGATCCCCGACTTGCCCACCAACTAAGGCCCCTCGGTCCTGTCGCCGCCGCCGCCGTTTCCGGATCAAACGACGT	OOOKIRJNSSONNONLSLQIMMRSKIOPLISMPQKQNILRLRQMORNKIQORQKSOJRJOMMQKSMJLQOMRKLIKJNJSONINILRJLQNOLMMNQMQMNSJMMSKKQJKROOMKRJJQPPSJMSQJMNJSPIMLRMOPPPNONRQOKP	AS:i:284	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:66C73T9	YS:i:284	YT:Z:CP	RG:Z:20Feb2019.08:11PM.read_gen.dna
chr19:59066330-59066516_THRZC67832XD	147	chr19	59066343	44	150M	=	59066332	-161	GGGGCCACTCACCGGAGCTGGGACGGACACCCGTGTCCGCTTCCGGGTCACGCTTCTCTTTCTGGGATCCCCGACTTGCCCACCAACTAAGGCCCCTCGGTCCTGTCGCCGCCGCCGCCGTTTCCGGATTAAACGACGTCACGTAACATG	IRIRKQMOJSLIMKPONSJINNMIRONOPKJPQSSKNOKMKNJMJSIQMSOJMNRJJPSNNNISRISPQMNMJMMKRPSOLJRQKQSJIRKMSQSOSJKOPLQQRIRNORJQOLRPIPSQJIONJKQMRLJRNPKJIJQLQOLORSQRNK	AS:i:284	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:34T104G10	YS:i:284	YT:Z:CP	RG:Z:20Feb2019.08:11PM.read_gen.dna
"""


class TestSeqUtils(unittest.TestCase):
    """Tests for TestSeqUtils."""

    APPRIS_REF = "GRCh38.chr21.fa.gz"
    CBS_REF = "CBS_pEZY3.fa"

    @classmethod
    def setUpClass(cls):
        """Setup for TestSeqUtilsSetup."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.abspath(os.path.join(cls.test_dir, "..", "test_data"))
        cls.cbs_ref = os.path.join(cls.test_data_dir, cls.CBS_REF)
        cls.cbs_ref_copy = os.path.join(cls.tempdir, cls.CBS_REF)
        copyfile(cls.cbs_ref, cls.cbs_ref_copy)
        index_reference(cls.cbs_ref_copy)
        cls.seq = "ATCGATTACG"

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestSeqUtilsSetup."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

    def test_reverse_complement(self):
        """Tests for proper reverse complement of string."""

        observed = su.reverse_complement(self.seq)
        expected = "CGTAATCGAT"
        self.assertEqual(observed, expected)

    def test_introduce_error(self):
        """Test for proper introduction of error into a sequence."""

        observed_seq = su.introduce_error(self.seq, 0.1)
        expected_num = 1
        observed_num = 0

        for (bp1, bp2) in zip(observed_seq, self.seq):
            if bp1 != bp2:
                observed_num += 1

        self.assertEqual(observed_num, expected_num)

    def test_introduce_indels(self):
        """Test for introduction of InDels in a read."""

        random.seed(9)
        wt_seq = "GGAGAAGGGCTTCGACCAGGCGCCCGTGGTGGATGAGGCGGGGGTAATCCTGGGAATGGTGACGCTTGGGAACATGCTCTCGTCCCTGCTTGCCGGGAAGGTGCA"
        indel_seq = su.introduce_indels(wt_seq, 0.03)
        fasta_lines = [">WT", wt_seq, ">InDels", indel_seq]
        has_indels = False

        with tempfile.NamedTemporaryFile(mode="w", suffix=".test.fa", delete=False, dir=self.tempdir) as test_fasta:
            test_fasta.write(fu.FILE_NEWLINE.join(fasta_lines) + fu.FILE_NEWLINE)
            test_fasta_fn = test_fasta.name

        ba = al.Bowtie2(
            f1=test_fasta_fn, output_dir=self.tempdir, config=al.BowtieConfig(
                self.cbs_ref_copy, False, 1, DEFAULT_QUALITY_OFFSET, "f"))

        with pysam.AlignmentFile(ba.output_bam, "rb") as af:
            for align_seg in af.fetch():
                if align_seg.query_name == "InDels":
                    has_indels = "I" in align_seg.cigarstring or "D" in align_seg.cigarstring

        fu.safe_remove((ba.output_bam, fu.add_extension(ba.output_bam, su.BAM_INDEX_SUFFIX),))
        self.assertTrue(has_indels)

    def test_extract_seq(self):
        """Test that we properly extract sequences from the genome."""

        observed = su.extract_seq("CBS_pEZY3", 3200, 3204, self.cbs_ref_copy).upper()
        expected = "TTCAC"
        self.assertEqual(observed, expected)


class TestSamtools(unittest.TestCase):
    """Tests for pysam wrappers of samtools commands."""

    @classmethod
    def setUpClass(cls):
        """Setup for TestSamtools."""

        cls.tempdir = tempfile.mkdtemp()

        with tempfile.NamedTemporaryFile(mode="w", suffix=".test.sam", delete=False, dir=cls.tempdir) as test_sam:
            test_sam.write(TEST_SAM)
            cls.test_sam_name = test_sam.name

        random.seed(9)

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestSamtools."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

    def test_sam_view(self):
        """Tests the pysam.view call."""

        # Just check that we can extract the body
        res = su.sam_view(self.test_sam_name, output_format="SAM")

        with open(res, "r") as out_sam:
            observed = out_sam.read()

        expected = fu.FILE_NEWLINE.join(
            [line for line in TEST_SAM.splitlines() if not line.startswith(su.SAM_HEADER_CHAR)]) + fu.FILE_NEWLINE

        fu.safe_remove((res,))
        self.assertEqual(expected, observed)

    def test_sam_view_with_flag(self):
        """Tests the pysam.view call with a flag."""

        # Check that we can extract the header
        res = su.sam_view(self.test_sam_name, None, "SAM", 0, "H")

        with open(res, "r") as out_sam:
            # Remove the newly added program header
            observed = out_sam.readlines()[0:3]
            observed = "".join(observed)

        # Do not copy the last header line which will have a new line added by samtools
        expected = fu.FILE_NEWLINE.join([
            line for line in TEST_SAM.splitlines() if line.startswith(su.SAM_HEADER_CHAR)]) + fu.FILE_NEWLINE

        fu.safe_remove((res,))
        self.assertEqual(expected, observed)

    def test_sort_bam(self):
        """Test that we properly sort a BAM."""

        # Randomly shuffle the lines
        with tempfile.NamedTemporaryFile(
                "w+", suffix=".unsorted.test.sam", delete=False, dir=self.tempdir) as test_unsorted_sam:

            sam_lines = TEST_SAM.splitlines()
            header_lines = []
            body_lines = []

            for sam_line in sam_lines:
                if sam_line.startswith(su.SAM_HEADER_CHAR):
                    header_lines.append(sam_line)
                else:
                    body_lines.append(sam_line)

            body_lines.reverse()
            unsorted_lines = header_lines + body_lines
            test_unsorted_sam.write(fu.FILE_NEWLINE.join(unsorted_lines) + fu.FILE_NEWLINE)
            out_sam = test_unsorted_sam.name

        # Convert to BAM
        out_bam = su.sam_view(out_sam)

        # Sort the BAM
        res = su.sort_bam(out_bam, output_format="SAM")

        with open(res, "r") as outsam:
            for i, line in enumerate(outsam):
                if i == 5:
                    start_pos = int(line.split(fu.FILE_DELIM)[3])
                if i == 14:
                    end_pos = int(line.split(fu.FILE_DELIM)[3])

        fu.safe_remove((out_sam, out_bam, res,))
        self.assertTrue(start_pos < end_pos)

    def test_index_bam(self):
        """Test we can create an index file on a sorted BAM."""

        res = su.sam_view(self.test_sam_name)
        su.index_bam(res)
        test_res = os.path.exists(fu.add_extension(res, su.BAM_INDEX_SUFFIX))
        fu.safe_remove((res, fu.add_extension(res, su.BAM_INDEX_SUFFIX),))
        self.assertTrue(test_res)

    def test_sort_and_index(self):
        """Test for proper creation and indexing of a coordinate-sorted BAM."""

        res = su.sort_and_index(self.test_sam_name)
        test_res = res.endswith(su.BAM_SUFFIX) and os.path.exists(fu.add_extension(res, su.BAM_INDEX_SUFFIX))
        fu.safe_remove((res, fu.add_extension(res, su.BAM_INDEX_SUFFIX),))
        self.assertTrue(test_res)


class TestFastaToFastq(unittest.TestCase):
    """Tests for analysis.seq_utils.fasta_to_fastq."""

    @classmethod
    def setUpClass(cls):
        """Tests for TestFastaToFastq."""

        cls.tempdir = tempfile.mkdtemp()

        with tempfile.NamedTemporaryFile(mode="w", suffix=".test.fa", delete=False, dir=cls.tempdir) as test_fasta:
            test_fasta.write(TEST_FASTA)
            cls.test_fasta = test_fasta.name

        with tempfile.NamedTemporaryFile(suffix=".out.fq", delete=False, dir=cls.tempdir) as test_fastq_file:
            output_fn = test_fastq_file.name
            su.fasta_to_fastq(cls.test_fasta, output_fn, max_len=None, snp_prob=0.1, letters=("A", "G",))

        with open(output_fn, "r") as output_fh:

            observed_seq = "".join([line.rstrip("\r\n") for num, line in enumerate(output_fh.readlines()) if num == 1])
            expected_seq = "ATGTTACGTCACGTCGTTTAATCCGGAAACGGCGGCGGCGGCGACAGGACCGAGGGGCCTTAGTTGGTGGGCAAGTCGGGGATCCCAGAAAGAGAAGCGTGACCCGGAAGCGGAAACGGGTGTCCGTCCCAGCTCCGGCCTGCCAGTGAGCTTCTACCATCATGGACCTATTGTTCGGGCGCCGGAAGACGCCAGAGGAGCTACTGCGGCAGAACCAGAGGGCCCTGAACCGTGCCATGCGGGAGCTGGACCGCGAGCGACAGAAACTAGAGACCCAGGAGAAGAAAATCATTGCAGACATTAAGAAGATGGCCAAGCAAGGCCAGATGGATGCTGTTCGCATCATGGCAAAAGACTTGGTGCGCACCCGGCGCTATGTGCGCAAGTTTGTATTGATGCGGGCCAACATCCAGGCTGTGTCCCTCAAGATCCAGACACTCAAGTCCAACAACTCGATGGCACAAGCCATGAAGGGTGTCACCAAGGCCATGGGCACCATGAACAGACAGCTGAAGTTGCCCCAGATCCAGAAGATCATGATGGAGTTTGAGCGGCAGGCAGAGATCATGGATATGAAGGAGGAGATGATGAATGATGCCATTGATGATGCCATGGGTGATGAGGAAGATGAAGAGGAGAGTGATGCTGTGGTGTCCCAGGTTCTGGATGAGCTGGGACTTAGCCTAACAGATGAGCTGTCGAACCTCCCCTCAACTGGGGGCTCGCTTAGTGTGGCTGCTGGTGGGAAAAAAGCAGAGGCCGCAGCCTCAGCCCTAGCTGATGCTGATGCAGACCTGGAGGAACGGCTTAAGAACCTGCGGAGGGACTGAGTGCCCCTGCCACTCCGAGATAACCAGTGGATGCCCAGGATCTTTTACCACAACCCCTCTGTAATAAAAGAGATTTGACACTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            cls.zipped_seqs = list(zip(observed_seq, expected_seq))

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestFastaToFastq."""

        fu.safe_remove((cls.tempdir,), force_remove=True)

    def test_fasta_to_fastq_proportion_error(self):
        """Tests correct error proportion."""

        observed_num = 0
        expected_num = 94

        for (bp1, bp2) in self.zipped_seqs:
            if bp1 != bp2:
                observed_num += 1

        self.assertEqual(observed_num, expected_num)

    def test_fasta_to_fastq_letters(self):
        """Tests introduction of specific bases as errors."""

        expected_num = 0
        observed_num = 0

        for (bp1, bp2) in self.zipped_seqs:
            # Check to ensure no complement of the letters has been introduced
            if bp1 != bp2 and bp1 in {"T", "C"}:
                observed_num += 1

        self.assertEqual(observed_num, expected_num)


class TestEnums(unittest.TestCase):
    """Tests for Enums."""
    
    def test_Strand(self):
        """Test strand equality with various inputs."""

        # Comes from AlignedSegment.is_reverse flag
        self.assertTrue(su.Strand("+") == su.Strand(False))
        
    def test_HumanContig_strings(self):
        """Test contig equality with string inputs."""

        self.assertTrue(su.HumanContig("chrX") == su.HumanContig("NC_000023"))

    def test_HumanContig_string_int(self):
        """Test contig equality with string and int inputs."""

        self.assertTrue(su.HumanContig("chr1") == su.HumanContig(1))

