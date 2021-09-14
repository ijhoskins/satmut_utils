#!/usr/bin/env/python
"""Tests for analysis.seq_utils."""

import os
import pysam
import random
import re
import tempfile
import unittest

import analysis.aligners as al
import analysis.seq_utils as su
import core_utils.file_utils as fu

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
@SQ	SN:chr1	LN:249250621
@SQ	SN:chr2	LN:243199373
@SQ	SN:chr3	LN:198022430
@SQ	SN:chr4	LN:191154276
@SQ	SN:chr5	LN:180915260
@SQ	SN:chr6	LN:171115067
@SQ	SN:chr7	LN:159138663
@SQ	SN:chrX	LN:155270560
@SQ	SN:chr8	LN:146364022
@SQ	SN:chr9	LN:141213431
@SQ	SN:chr10	LN:135534747
@SQ	SN:chr11	LN:135006516
@SQ	SN:chr12	LN:133851895
@SQ	SN:chr13	LN:115169878
@SQ	SN:chr14	LN:107349540
@SQ	SN:chr15	LN:102531392
@SQ	SN:chr16	LN:90354753
@SQ	SN:chr17	LN:81195210
@SQ	SN:chr18	LN:78077248
@SQ	SN:chr20	LN:63025520
@SQ	SN:chrY	LN:59373566
@SQ	SN:chr19	LN:59128983
@SQ	SN:chr22	LN:51304566
@SQ	SN:chr21	LN:48129895
@RG	ID:20Feb2019.08:11PM.read_gen.dna
@PG	ID:bowtie2	PN:bowtie2	VN:2.3.4.3	CL:"/Users/ianhoskins/.conda/envs/Cenik_lab/bin/bowtie2-align-s --wrapper basic-0 --local --rg-id 20Feb2019.08:11PM.read_gen.dna -x /usr/local/bin/hg19.fa -1 read_gen_results/20Feb2019.08:11PM.read_gen.dna.r1.fq -2 read_gen_results/20Feb2019.08:11PM.read_gen.dna.r2.fq"
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


class TestSeqUtilsSetup(object):

    @classmethod
    def setUpClass(cls):
        """Setup for TestSeqUtils."""

        cls.seq = "ATCGATTACG"

        # We'll need another BED file to check intersection, etc.
        cls.test_bed_b_str = "".join(TEST_BED.split("\n")[0]) + fu.FILE_NEWLINE

        with tempfile.NamedTemporaryFile("w", suffix=".test.fa", delete=False) as test_fasta, \
                tempfile.NamedTemporaryFile("w", suffix=".test.a.bed", delete=False) as test_bed_a, \
                tempfile.NamedTemporaryFile("w", suffix=".test.b.bed", delete=False) as test_bed_b, \
                tempfile.NamedTemporaryFile("w", suffix=".test.sam", delete=False) as test_sam:

            test_fasta.write(TEST_FASTA)
            test_bed_a.write(TEST_BED)
            test_bed_b.write(cls.test_bed_b_str)
            test_sam.write(TEST_SAM)

            cls.test_fasta = test_fasta.name
            cls.test_bed_a = test_bed_a.name
            cls.test_bed_b = test_bed_b.name
            cls.test_sam = test_sam.name

    @classmethod
    def tearDownClass(cls):
        test_files = [fn for fn in os.listdir(tempfile.gettempdir()) if re.search("test", fn)]
        fu.safe_remove(tuple(test_files))


# Note multiple inheritance may not work with more than one non-unittest.TestCase superclass:
# https://stackoverflow.com/questions/18896877/python-multiple-inheritance-unittest

class TestSeqUtils(TestSeqUtilsSetup, unittest.TestCase):
    """Tests for TestSeqUtils."""

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

        wt_seq = "GACGGACACCCGTTTCCGCTTCCGGGTCACGCTTCTCTTTCTGGGATCCCCGACTTGCCCACCAACTAAGGCCCCTCGGTCCTGTCGCCGCCGCCGCCGTTTCCGGATTAAACGACGTGACGTAACATGCCCCGCCCGCACCCGGAACGT"
        indel_seq = su.introduce_indels(wt_seq, 0.02)
        fasta_lines = [">WT", wt_seq, ">InDels", indel_seq]

        has_indels = False

        with tempfile.NamedTemporaryFile("w", suffix=".test.fa") as test_fasta, \
                tempfile.NamedTemporaryFile("w+b", suffix=".test.bam") as test_bam:

            test_fasta.write(fu.FILE_NEWLINE.join(fasta_lines) + fu.FILE_NEWLINE)
            fu.flush_files((test_fasta,))
            
            al.Bowtie2(f1=test_fasta.name, output_bam=test_bam.name, config=al.BowtieConfig(su.GENOMIC_FASTA, False, "f"))
            fu.flush_files((test_bam,))

            with pysam.AlignmentFile(test_bam.name, "rb") as af:
                for align_seg in af.fetch(until_eof=True):

                    if align_seg.is_unmapped:
                        continue

                    if align_seg.query_name == "InDels":
                        has_indels = "I" in align_seg.cigarstring or "D" in align_seg.cigarstring

        self.assertTrue(has_indels)

    def test_bq_int_to_ascii(self):
        """Test for proper conversion between Illumina-offset int to ASCII format."""

        observed = su.bq_int_to_ascii(0)
        expected = "!"
        self.assertEqual(observed, expected)

    def test_bq_ascii_to_int(self):
        """Test for proper conversion between Illumina-offset ASCII to int."""

        observed = su.bq_ascii_to_int("!")
        expected = 0
        self.assertEqual(observed, expected)

    def test_extract_seq(self):
        """Test that we properly extract sequences from the genome."""

        observed = su.extract_seq("chr19", 59065436, 59065475)
        expected = "CTTAATGTCTGCAATGATTTTCTTCTCCTGGGTCTCTAGT"
        self.assertEqual(observed, expected)


class TestSamtools(TestSeqUtilsSetup, unittest.TestCase):
    """Tests for pysam wrappers of samtools commands."""

    @classmethod
    def setUpClass(cls):
        """Setup for TestSamtools."""

        with tempfile.NamedTemporaryFile("w+", suffix=".test.sam", delete=False) as test_sam:
            test_sam.write(TEST_SAM)
            cls.test_sam_name = test_sam.name

        random.seed(9)

    def test_sam_view(self):
        """Tests the pysam.view call."""

        # Just check that we can extract the body
        res = su.sam_view(self.test_sam_name, output_format="SAM")

        with open(res, "r") as out_sam:
            obs = out_sam.read()

        exp = fu.FILE_NEWLINE.join(
            [line for line in TEST_SAM.splitlines() if not line.startswith(su.SAM_HEADER_CHAR)]) + fu.FILE_NEWLINE

        self.assertEqual(obs, exp)

    def test_sam_view_with_flag(self):
        """Tests the pysam.view call with a flag."""

        # Check that we can extract the header
        res = su.sam_view(self.test_sam_name, None, "SAM", 0, "H")

        with open(res, "r") as out_sam:
            obs = out_sam.read()

        exp = fu.FILE_NEWLINE.join(
            [line for line in TEST_SAM.splitlines() if line.startswith(su.SAM_HEADER_CHAR)]) + fu.FILE_NEWLINE

        self.assertEqual(obs, exp)

    def test_sort_bam(self):
        """Test that we properly sort a BAM."""

        # Randomly shuffle the lines
        with tempfile.NamedTemporaryFile("w+", suffix=".unsorted.test.sam", delete=False) as test_unsorted_sam:

            sam_lines = TEST_SAM.splitlines()
            header_lines = []
            body_lines = []

            for sam_line in sam_lines:
                if sam_line.startswith(su.SAM_HEADER_CHAR):
                    header_lines.append(sam_line)
                else:
                    body_lines.append(sam_line)

            shuffled_body_lines = random.sample(body_lines, len(body_lines))
            reassembled_lines = header_lines + shuffled_body_lines

            test_unsorted_sam.write(fu.FILE_NEWLINE.join(reassembled_lines) + fu.FILE_NEWLINE)
            out_sam = test_unsorted_sam.name

        # Convert to BAM
        out_bam = su.sam_view(out_sam)

        # Sort the BAM
        res = su.sort_bam(out_bam, output_format="SAM")

        with open(res, "r") as outsam:
            obs = outsam.read()

        self.assertEqual(obs, TEST_SAM)

    def test_index_bam(self):
        """Test we can create an index file on a sorted BAM."""

        res = su.sam_view(self.test_sam_name)

        su.index_bam(res)

        self.assertTrue(os.path.exists("{}.{}".format(res, "bai")))

    def test_sort_and_index(self):
        """Test for proper creation and indexing of a coordinate-sorted BAM."""

        res = su.sort_and_index(self.test_sam_name)

        self.assertTrue(res.endswith(su.BAM_SUFFIX) and os.path.exists("{}.{}".format(res, "bai")))

    def test_merge_ams(self):
        """Test that we can merge BAMs."""
        # TODO
        # Not essential at this time
        pass


class TestFastaToFastq(TestSeqUtilsSetup, unittest.TestCase):
    """ Tests for su.fasta_to_fastq."""

    @classmethod
    def setUpClass(cls):
        """Tests for TestFastaToFastq."""

        super(TestFastaToFastq, cls).setUpClass()

        with tempfile.NamedTemporaryFile(suffix=".out.fq", delete=False) as test_fastq_file:
            output_fn = test_fastq_file.name
            su.fasta_to_fastq(cls.test_fasta, output_fn, max_len=None, snp_prob=0.1, letters=("A", "G",))

        with open(output_fn, "r") as output_fh:

            observed_seq = "".join([line.rstrip("\r\n") for num, line in enumerate(output_fh.readlines()) if num == 1])
            expected_seq = "ATGTTACGTCACGTCGTTTAATCCGGAAACGGCGGCGGCGGCGACAGGACCGAGGGGCCTTAGTTGGTGGGCAAGTCGGGGATCCCAGAAAGAGAAGCGTGACCCGGAAGCGGAAACGGGTGTCCGTCCCAGCTCCGGCCTGCCAGTGAGCTTCTACCATCATGGACCTATTGTTCGGGCGCCGGAAGACGCCAGAGGAGCTACTGCGGCAGAACCAGAGGGCCCTGAACCGTGCCATGCGGGAGCTGGACCGCGAGCGACAGAAACTAGAGACCCAGGAGAAGAAAATCATTGCAGACATTAAGAAGATGGCCAAGCAAGGCCAGATGGATGCTGTTCGCATCATGGCAAAAGACTTGGTGCGCACCCGGCGCTATGTGCGCAAGTTTGTATTGATGCGGGCCAACATCCAGGCTGTGTCCCTCAAGATCCAGACACTCAAGTCCAACAACTCGATGGCACAAGCCATGAAGGGTGTCACCAAGGCCATGGGCACCATGAACAGACAGCTGAAGTTGCCCCAGATCCAGAAGATCATGATGGAGTTTGAGCGGCAGGCAGAGATCATGGATATGAAGGAGGAGATGATGAATGATGCCATTGATGATGCCATGGGTGATGAGGAAGATGAAGAGGAGAGTGATGCTGTGGTGTCCCAGGTTCTGGATGAGCTGGGACTTAGCCTAACAGATGAGCTGTCGAACCTCCCCTCAACTGGGGGCTCGCTTAGTGTGGCTGCTGGTGGGAAAAAAGCAGAGGCCGCAGCCTCAGCCCTAGCTGATGCTGATGCAGACCTGGAGGAACGGCTTAAGAACCTGCGGAGGGACTGAGTGCCCCTGCCACTCCGAGATAACCAGTGGATGCCCAGGATCTTTTACCACAACCCCTCTGTAATAAAAGAGATTTGACACTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            cls.zipped_seqs = list(zip(observed_seq, expected_seq))

    def test_fasta_to_fastq_proportion_error(self):
        """Test for the correct error proportion."""

        observed_num = 0
        expected_num = 94

        for (bp1, bp2) in self.zipped_seqs:
            if bp1 != bp2:
                observed_num += 1

        self.assertEqual(observed_num, expected_num)

    def test_fasta_to_fastq_letters(self):
        """Test that we can introduce specific bases as errors."""

        expected_num = 0
        observed_num = 0

        for (bp1, bp2) in self.zipped_seqs:
            # Check to ensure no complement of the letters has been introduced
            if bp1 != bp2 and bp1 in {"T", "C"}:
                observed_num += 1

        self.assertEqual(observed_num, expected_num)


class TestAddErrorToFasta(TestSeqUtilsSetup, unittest.TestCase):
    """Tests for su.add_error_to_fasta."""
    # TODO
    # Not a priority right now as we already test the error functions by themselves
    pass


class TestEnums(unittest.TestCase):
    """Tests for Enums."""
    
    def test_Strand(self):
        """Test strand equality with various inputs."""
        
        self.assertTrue(su.Strand("+") == su.Strand(False))
        
    def test_HumanContig(self):
        """Test contig equality with various inputs."""

        self.assertTrue(su.HumanContig("chrX") == su.HumanContig("NC_000023"))

