#!/usr/bin/env python3
"""Tests for analysis.read_editor"""

import collections
import copy
import gzip
import pysam
import shutil
import tempfile
import unittest

import analysis.read_editor as ed
from analysis.references import index_reference
from analysis.seq_utils import sort_and_index, COORD_FORMAT, DEFAULT_MAPQ, ReadMate, MASKED_BQ, SAM_EDITED_TAG, FASTQ_QNAME_CHAR, BAM_INDEX_SUFFIX
import core_utils.file_utils as fu
from core_utils.vcf_utils import get_variant_type
from satmut_utils.definitions import *

tempfile.tempdir = DEFAULT_TEMPDIR

TEST_VCF = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|>
##INFO=<ID=AAM_AA_CHANGE,Number=.,Type=String,Description="Can_Cenik_lab.Ian_Hoskins.bioinfo.analysis.coordinate_mapper.AminoAcidMapper comma-delimited amino acid changes.">
##INFO=<ID=VARTYPE,Number=.,Type=String,Description="Variant type. MNPs are confined to a codon. Haplotypes span codons.">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency in range (0,1).">
##INFO=<ID=IE,Number=1,Type=String,Description="Introduce sequencing error to a specific strand, either + or -.">
##INFO=<ID=IR,Number=1,Type=String,Description="Introduce error (sample strand chemical change) to a specific read mate, either R1 or R2.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	795	.	GAC	ATG	.	.	AAM_AA_CHANGE=p.D179M;VARTYPE=tri_nt_MNP;AF=1.0
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	804	.	C	G	.	.	AAM_AA_CHANGE=p.R182G;VARTYPE=snp;AF=0.25
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	804	.	C	T	.	.	AAM_AA_CHANGE=p.R182W;VARTYPE=snp;AF=0.25
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	812	.	GG	AT	.	.	AAM_AA_CHANGE=p.G185M;VARTYPE=di_nt_MNP;AF=0.25
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	955	.	ACT	A	.	.	VARTYPE=del;AF=1.0
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1094	.	T	TG	.	.	VARTYPE=ins;AF=1.0
"""

TEST_SAM = """@HD	VN:1.0	SO:coordinate
@SQ	SN:ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	LN:2605
0000005815	163	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	241	42	148M	=	241	-148	ACCATCTGTCCGGTCCCAGCATGCCTTCTGAGACCCCCCAGGCAGAAGTGGGGCCCACAGGCTGCCCCCACCGCTCAGGGCCACACTCGGCGAAGGGGAGCCTGGAGAAGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTG	SSSSSSSSSSLOKROLRPMKIQQPJRIQJSNLIINKMSJMQNIINKSPQNNLNQNPNPJPRSOPSKPLOPOIPMSSPLOINSJQMMQMKLSNKQLIPIQRNROKPRKMJKIKRLIMSJMRRRPRRRPQNLOIKSJMIRSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:148	YS:i:0	YT:Z:CP
0000005815	83	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	241	42	148M	=	241	-148	ACCATCTGTCCGGTCCCAGCATGCCTTCTGAGACCCCCCAGGCAGAAGTGGGGCCCACAGGCTGCCCCCACCGCTCAGGGCCACACTCGGCGAAGGGGAGCCTGGAGAAGGGGTCCCCAGAGGATAAGGAAGCCAAGGAGCCCCTGTG	SSSSSSSSSSINKKNIRJRIJONJRNJIRQISRQIMIJMSLQJMIRRMISPKPQJJMIOQKQOQJKPIMPKJPMSNNSILPILKNRRLIJPQKMMKQIOLPSMIMKLRSILNIMNPISOPPRLRQSNRPMQPJJMRRRSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:148	YS:i:0	YT:Z:CP
0000253933	163	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	566	42	150M	=	569	153	GTGTGAGCTCTTGGCCAAGTGTGAGTTCTTCAACGCGGGCGGGAGCGTGAAGGACCGCATCAGCCTGCGGATGATTGAGGATGCTGAGCGCGACGGGACGCTGAAGCCCGGGGACACGATTATCGAGCCGACATCCGGGAACACCGGGAT	SSSSSSSSSSNJMMKKPJQMSSSMINMMKOSQKPKINROPPORPQLKMMIJNSSORRIIPRKRIMIPOMOPMJSNLIQQLORJRJILRJMLPSILQSKMLKKSLJPNMMQIMLKRNQMNSPKPOPQIIINIPMOSRIKROPSLSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:0	YT:Z:CP
0000253933	83	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	569	42	150M	=	566	-153	TGAGCTCTTGGCCAAGTGTGAGTTCTTCAACGCGGGCGGGAGCGTGAAGGACCGCATCAGCCTGCGGATGATTGAGGATGCTGAGCGCGACGGGACGCTGAAGCCCGGGGACACGATTATCGAGCCGACATCCGGGAACACCGGGATCGG	SSSSSSSQIKKSQINKPSPSLRILIONLIMNKJNQPRINPOJOKONRPSKSQININOKOPONIPINNOIRSPOPJSSQPMPPLOOSKNOSJRINIKJNOMRLSSKJLLJLONQOJSMNPMJLMIMMPLSIIMMLLLONMSSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:0	YT:Z:CP
0000224875	163	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	680	42	150M	=	683	153	CACGATTATCGAGCCGACATCCGGGAACACCGGGATCGGGCTGGCCCTGGCTGCGGCAGTGAGGGGCTATCGCTGCATCATCGTGATGCCAGAGAAGATGAGCTCCGAGAAGGTGGACGTGCTGTGGGCACTGGGGGCTGAGATTGTGAG	SSSSSSSSSSOKROSNOMKISLLKKLNPLIJQOSJJMRLQIJJQRSRNMKLLOLJSJRNJRSPPPLMRQKQMIRQRLJOPSOJMKLRRMJIJOOJQLSPLKKMJSNNQLOMJSOJKRRMPPLOOIQPNRQJSLJSOKOKMMMMSSSSSSS	AS:i:-4	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:124C25	YS:i:0	YT:Z:CP
0000224875	83	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	683	42	150M	=	680	-153	GATTATCGAGCCGACATCCGGGAACACCGGGATCGGGCTGGCCCTGGCTGCGGCAGTGAGGGGCTATCGCTGCATCATCGTGATGCCAGAGAAGATGAGCTCCGAGAAGGTGGACGTGCTGCGGGCACTGGGGGCTGAGATTGTGAGGAC	SSSSSSSNSRKQNMLQIKNLNJPKLIINNNQRNIROIQJMJLLKOQSQIIKSOMOOQQIRRRNNMLIOLKSJQQJKKLQRSIQPLJLRMLPKOPMNPLIRLLIJNQLKISPOSKKIMIIKSMKMJKSKKORKISNPROOKSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:-4	YT:Z:CP
0000001369	163	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	794	42	150M	=	797	153	GGACGTGCTGAGGGCACTGGGGGCTGAGATTGTGAGGACGCCCACCAATGCCAGGTTCGACTCCCCGGAGTCACACGTGGGGGTGGCCTGGCGGCTGAAGAACGAAATCCCCAATTCTCACATCCTAGACCAGTACCGCAACGCCAGCAA	SSSSSSSSSSMOINJNKNNLLRNLPRRNNKLLNMRNNLJPMLSOINMSILMLMRLSLPNKLLRNJMSKJJIPQMJNRLQOJQLOSRORNQSQRPOPOSIPMPILKNNLPMNIQJPOQPILOOMORLPILSONMKLRLPPRIQPSSSSSSS	AS:i:-4	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:10C139	YS:i:0	YT:Z:CP
0000051899	163	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	794	42	17M2D131M	=	797	153	GGACGTGCTGCGGGCACGGGGCTGAGATTGTGAGGACGCCCACCAATGCCAGGTTCGACTCCCCGGAGTCACACGTGGGGGTGGCCTGGCGGCTGAAGAACGAAATCCCCAATTCTCACATCCTAGACCAGTACCGCAACGCCAGCAA	SSSSSSSSSSKPLJJQLSLOMIPSIPMILSORJONIQLPKJSMJRKNIIJRPSSIOJSOOMSKQKOQMNLIOQRSQIPILKSRQPMNIRMSNKQOOJQSLMQRNJJPJOJLKRSOKSQIPSKKQPSPSLJJKSPRMNJNJLSSSSSSS	AS:i:-14	XN:i:0	XM:i:0	XO:i:1	XG:i:2	NM:i:2	MD:Z:17^TG131	YS:i:0	YT:Z:CP
0000225689	99	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	794	42	17M1I133M	=	797	153	GGACGTGCTGCGGGCACTTGGGGGCTGAGATTGTGAGGACGCCCACCAATGCCAGGTTCGACTCCCCGGAGTCACACGTGGGGGTGGCCTGGCGGCTGAAGAACGAAATCCCCAATTCTCACATCCTAGACCAGTACCGCAACGCCAGCAA	SSSSSSSSSSPNIIJKQLRRPOSSJKSMKOQKPIQKONRJNQRRPJNOQLOMPSPMSSNQRJIRPOOJIPLNSKNLMKPNSRLJLRMMJKIRPJIRKSJKMLSQOQQSJLPMSIPRLLQMIIROKOPNROSOPSKJJSPPSQIQSSSSSSS	AS:i:-10	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:150	YS:i:0	YT:Z:CP
0000001369	83	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	797	42	150M	=	794	-153	CGTGCTGCGGGCACTGGGGGCTGAGATTGTGAGGACGCCCACCAATGCCAGGTTCGACTCCCCGGAGTCACACGTGGGGGTGGCCTGGCGGCTGAAGAACGAAATCCCCAATTCTCACATCCTAGACCAGTACCGCAACGCCAGCAACCC	SSSSSSSIRMRPPMKORMQQKPOQRJLQPLOKSQJOIMJJILKPKJQKNRIJJLJPRMNIRKJMJKKOJQPSKSSNJIILLLSRMPPNLIJOOJNNMMLPIOKLPSMQSQOQKPOSSQKISLQSLLRRKNSMRORRNJKOSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:-4	YT:Z:CP
0000051899	83	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	797	42	150M	=	794	-153	CGTGCTGCGGGCACTGGGGGCTGAGATTGTGAGGACGCCCACCAATGCCAGGTTCGACTCCCCGGAGTCACACGTGGGGGTGGCCTGGCGGCTGAAGAACGAAATCCCCAATTCTCACATCCTAGACCAGTACCGCAACGCCAGCAACCC	SSSSSSSNMRMMNRNPRQKRNRIOMSLPIIISLOQIOKKOOOMLOLQINPNMNNQOIMNJOQRSMPQISIIQINSJORMRROLSJQMLISKMLPMJMMSQINOOQOOSQMJMIJJLRSKPNSSSOMQPNSNKQJLIJJKKSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:-14	YT:Z:CP
0000225689	147	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	797	42	150M	=	794	-153	CGTGCTGCGGGCACTGGGGGCTGAGATTGTGAGGACGCCCACCAATGCCAGGTTCGACTCCCCGGAGTCACACGTGGGGGTGGCCTGGCGGCTGAAGAACGAAATCCCCAATTCTCACATCCTAGACCAGTACCGCAACGCCAGCAACCC	SSSSSSSNKMQOSOKSKJRQJJPKLKPRPOISPKNLQSSKLMLIKKPOSOQONSMILKSQISIJJMQOMNKQPQMQRIQRQQQSSJISSMSSQJILPIPMLOLQNMIMIRLQJSNKLKQJOSPRNMSPJKQJIQMKSIPSSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:-10	YT:Z:CP
0000003129	163	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	908	42	147M	=	908	-147	TTCTCACATCCTAGACCAGTACCGCAACGCCAGCAACCCCCTGGCTCACTACGACACCACCGCTGATGAGATCCTGCAGCAGTGTGATGGGAAGCTGGACATGCTGGTGGCTTCAGTGGGCACGGGCGGCACCATCACGGGCATTGC	SSSSSSSSSSOLPSKLNQOJLLNMPKLJJJNKQKOMMOKQPNLNOOLMJOMSQKIIQLKNLIPRQKNQIQRQIMJQSIQONKSIMONLJQSLMQOLNJSONKMNNLNNPKLQKPPIOMPQQQKOKPMSRPRKONMNNSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:147	YS:i:0	YT:Z:CP
0000003129	83	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	908	42	147M	=	908	-147	TTCTCACATCCTAGACCAGTACCGCAACGCCAGCAACCCCCTGGCTCACTACGACACCACCGCTGATGAGATCCTGCAGCAGTGTGATGGGAAGCTGGACATGCTGGTGGCTTCAGTGGGCACGGGCGGCACCATCACGGGCATTGC	SSSSSSSSSSIJJMPQKJKKINNNMKNPOQPJJNMLIRNKIRQPJNRPRSOQQMIMOIKNPLILONPKONSOOQJRSIQMIPNQSQSSNSKSOONKJJSQRIPOJQONPOQJRJSKRNMRQKRLNSQILNSQKOOIRSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:147	YS:i:0	YT:Z:CP
0000323491	99	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1016	42	150M	=	1019	153	GGCTTCAGTGGGCACGGGCGGCACCATCACGGGCATTGCCAGGAAGCTGAAGGAGAAGTGTCCTGGATGCAGGATCATTGGGGTGGATCCCGAAGGGTCCATCCTCGCAGAGCCGGAGGAGCTGAACCAGACGGAGCAGACAACCTACGA	SSSSSSSSSSNIILOQONPRQOPJOIOLNJJIOPNKQNPOKQRRSRJMPKSSPQQNJOMIJMKNIQRQPRSRIMNQNPONMRPKROKNORRRQQIISIQRONILSJKKRPSSIPKILKIMKQPNNOLJNKRIORKIMPMIPIKSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:0	YT:Z:CP
0000323491	147	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1019	42	150M	=	1016	-153	TTCAGTGGGCACGGGCGGCACCATCACGGGCATTGCCAGGAAGCTGAAGGAGAAGTGTCCTGGATGCAGGATCATTGGGGTGGATCCCGAAGGGTCCATCCTCGCAGAGCCGGAGGAGCTGAACCAGACGGAGCAGACAACCTACGAGGT	SSSSSSSIRJLKSRNOKILLMKJMQKINJQPNSJIQNIJMLJMRPPRQIRRIOQLRKKNQMPISMOJNJOOKKPJOSMNMSKSOONJNKQRSLQPOSMLQMRLRRSOIPLKIJMNKLMJJMKQQKPPOKJLNLSJNLKSMSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:0	YT:Z:CP
0000214705	163	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1358	42	150M	=	1361	153	GGGCCAGCGCTGCGTGGTCATTCTGCCCGACTCAGTGCGGAACTACATGACCAAGTTCCTGAGCGACAGGTGGATGCTGCAGAAGGGCTTTCTGAAGGAGGAGGACCTCACGGAGAAGAAGCCCTGGTGGTGGCACCTCCGTGTTCAGGA	SSSSSSSSSSILOKJNNLJKRORJMSOLMLQMNRSOQQIPKJRKOQKMKMRNKIMNIIRIJPNQIIPNKQRNQMIMRRJLPSIMJRNNNNPKNPPIJJSKRRRKRNOSLQQOOSPNKMOSJLNOSONNJRQONORIRPMSOMLSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:0	YT:Z:CP
0000214705	83	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1361	42	150M	=	1358	-153	CCAGCGCTGCGTGGTCATTCTGCCCGACTCAGTGCGGAACTACATGACCAAGTTCCTGAGCGACAGGTGGATGCTGCAGAAGGGCTTTCTGAAGGAGGAGGACCTCACGGAGAAGAAGCCCTGGTGGTGGCACCTCCGTGTTCAGGAGCT	SSSSSSSKKMRPKISSNKMISKISRQPOJKJJQLNRKRIQPMJLPMPJLRKJIJQSKJJPSIQLQLRMQONOMLKNSQLKKOOSPLNJRIJSIIIQQRQLRJJPKRNSOKPRJNMSLSKSLMRLRLPQIQIIPIPKSOIKSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:150	YS:i:0	YT:Z:CP
0000099776	99	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1805	42	132M	=	1805	-132	GATCCAGTACCACAGCACCGGGAAGTCCAGTCAGCGGCAGATGGTGTTCGGGGTGGTCACCGCCATTGACTTGCTGAACTTCGTGGCCGCCCAGGAGCGGGACCAGAAGTGAAGTCCGGAGCGCTGGGCGGT	SSSSSSSSSSONQJOQQSQKNJSKRMNQKJLLOSJQJKLMQPLKKLSLRQLPPIPSLJSSJJMMOJINMRJJLSJMQLMPNNPSJQNNILPJKJNLNKPNLRJRSROMQRKKSRPRKIPLIQSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:132	YS:i:0	YT:Z:CP
0000099776	147	ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1805	42	132M	=	1805	-132	GATCCAGTACCACAGCACCGGGAAGTCCAGTCAGCGGCAGATGGTGTTCGGGGTGGTCACCGCCATTGACTTGCTGAACTTCGTGGCCGCCCAGGAGCGGGACCAGAAGTGAAGTCCGGAGCGCTGGGCGGT	SSSSSSSSSSINNIKMSNPPOONNOMMQIQOMRRSNKPLKLILMPKPPNMQNSLPRMRILNMIMIJNSLNJMOOQOPKRJNILIQPKQRNMOSNOQNJRKNONPSJLRRQIIMKOQKPOJMKSSSSSSSSSS	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:132	YS:i:0	YT:Z:CP
"""

TEST_PRIMERS = """ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	235	250	CBS_target	.	+
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	378	393	CBS_target	.	-
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	560	575	CBS_target	.	+
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	674	689	CBS_target	.	+
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	708	723	CBS_target	.	-
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	788	803	CBS_target	.	+
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	822	837	CBS_target	.	-
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	902	917	CBS_target	.	+
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	936	951	CBS_target	.	-
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1010	1025	CBS_target	.	+
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1044	1059	CBS_target	.	-
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1158	1173	CBS_target	.	-
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1352	1367	CBS_target	.	+
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1500	1515	CBS_target	.	-
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1799	1814	CBS_target	.	+
ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|	1926	1941	CBS_target	.	-
"""


class TestReadEditor(unittest.TestCase):
    """Tests for ReadEditor."""

    CBS_REF = "CBS.fa"

    @classmethod
    def setUpClass(cls):
        """Set up for TestReadEditor."""

        cls.tempdir = tempfile.mkdtemp()
        cls.test_dir = os.path.dirname(__file__)
        cls.test_data_dir = os.path.abspath(os.path.join(cls.test_dir, "..", "test_data"))
        cls.cbs_ref = os.path.join(cls.test_data_dir, cls.CBS_REF)

        # Copy the reference to the tempdir and index it
        cls.cbs_ref_copy = os.path.join(cls.tempdir, cls.CBS_REF)
        shutil.copyfile(cls.cbs_ref, cls.cbs_ref_copy)
        index_reference(cls.cbs_ref_copy)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".test.sam", dir=cls.tempdir) as test_sam, \
                tempfile.NamedTemporaryFile(
                    mode="w", suffix=".test.vcf", delete=False, dir=cls.tempdir) as test_vcf, \
                tempfile.NamedTemporaryFile(
                    mode="w", suffix=".test.primers.bed", delete=False, dir=cls.tempdir) as test_primers:

            test_sam.write(TEST_SAM)
            fu.flush_files((test_sam,))
            cls.test_bam = sort_and_index(test_sam.name)

            test_vcf.write(TEST_VCF)
            cls.test_vcf = test_vcf.name

            test_primers.write(TEST_PRIMERS)
            cls.test_primers = test_primers.name

        # Store some aligned segments for testing editing objects
        with pysam.AlignmentFile(cls.test_bam, "rb") as test_bam:
            for i, align_seg in enumerate(test_bam.fetch(until_eof=True)):
                if i == 0:
                    cls.test_align_seg = align_seg
                break

        # Force edit of all variants despite their frequency sum exceeding 1
        cls.ed = ed.ReadEditor(
            cls.test_bam, variants=cls.test_vcf, ref=cls.cbs_ref_copy, primers=cls.test_primers,
            output_dir=cls.tempdir, output_prefix="test_editor", buffer=3, force_edit=True)

        cls.observed_edit_configs = cls.ed._get_edit_configs()

        # The qnames at POS 804 are 0000001369, 000224875, 000225689, 0000051899
        cls.test_variant_config = ed.VARIANT_CONFIG_TUPLE(
            type=get_variant_type("C", "G"),
            contig="ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|",
            pos=804, ref="C", alt="G", af=0.5)

        # The qname at POS 795 is 0000224875
        cls.test_variant_config_af1 = ed.VARIANT_CONFIG_TUPLE(
            type=get_variant_type("GAC", "ATG"),
            contig="ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|",
            pos=795, ref="GAC", alt="ATG", af=1.0)

        cls.contig = "ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|"

        # Single-stranded errors introduced into 000224875 and 0000001369 R2s at position 804 to check REF base
        # A del was introduced into 0000051899 R2, and an ins into 0000225689 R1, to check InDels

    @classmethod
    def tearDownClass(cls):
        """Tear down for TestReadEditor."""

        fu.safe_remove((cls.tempdir, cls.test_bam, fu.add_extension(cls.test_bam, BAM_INDEX_SUFFIX)), force_remove=True)

    def test_get_variant_configs(self):
        """Test that variant configurations are extracted from the input VCF."""

        expected = [
            ed.VARIANT_CONFIG_TUPLE(
                type=get_variant_type("GAC", "ATG"), contig=self.contig, pos=795, ref="GAC", alt="ATG", af=1.0),
            ed.VARIANT_CONFIG_TUPLE(
                type=get_variant_type("C", "G"), contig=self.contig, pos=804, ref="C", alt="G", af=0.25),
            ed.VARIANT_CONFIG_TUPLE(
                type=get_variant_type("C", "T"), contig=self.contig, pos=804, ref="C", alt="T", af=0.25),
            ed.VARIANT_CONFIG_TUPLE(
                type=get_variant_type("GG", "AT"), contig=self.contig, pos=812, ref="GG", alt="AT", af=0.25),
            ed.VARIANT_CONFIG_TUPLE(
                type=get_variant_type("ACT", "A"), contig=self.contig, pos=955, ref="ACT", alt="A", af=1.0),
            ed.VARIANT_CONFIG_TUPLE(
                type=get_variant_type("T", "TG"), contig=self.contig, pos=1094, ref="T", alt="TG", af=1.0),
        ]

        observed = self.ed._get_variant_configs()

        self.assertEqual(expected, observed)

    def test_get_edit_qnames_enough_amenable(self):
        """Tests editable read names are returned for proper allele frequencies."""

        expected = 2

        amenable_qnames = {0, 1, 2}
        total_amenable_qnames = 4
        qnames_to_edit = self.ed._get_edit_qnames(amenable_qnames, self.test_variant_config, total_amenable_qnames)
        observed = len(qnames_to_edit)

        self.assertEqual(expected, observed)

    def test_get_edit_qnames_none_amenable(self):
        """Tests that no read names are returned if none are editable at the position."""

        expected = 0
        amenable_qnames = set()
        total_amenable_qnames = 4
        qnames_to_edit = self.ed._get_edit_qnames(amenable_qnames, self.test_variant_config, total_amenable_qnames)
        observed = len(qnames_to_edit)

        self.assertEqual(expected, observed)

    def test_get_edit_qnames_bump(self):
        """Tests editable read names are returned."""

        expected = 1
        amenable_qnames = {0, 1, 2}
        total_amenable_qnames = 4

        variant_config = ed.VARIANT_CONFIG_TUPLE(
            type=get_variant_type("C", "G"),
            contig="ENST00000398165.7|ENSG00000160200.17|OTTHUMG00000086834.7|OTTHUMT00000195525.1|CBS-204|CBS|2605|protein_coding|",
            pos=804, ref="C", alt="G", af=0.1)

        qnames_to_edit = self.ed._get_edit_qnames(amenable_qnames, variant_config, total_amenable_qnames)
        observed = len(qnames_to_edit)

        self.assertEqual(expected, observed)

    def test_iterate_over_pileup_reads_caf(self):
        """Tests that edit configs are appended to the edit dictionary and the expected truth frequency is returned."""

        expected = 1.0

        observed_edit_configs = dict()
        with pysam.AlignmentFile(self.ed.editor_preprocessor.edit_background, "rb") as edited_background_af:

            # For purpose of testing change just the input coordinate
            target = COORD_FORMAT.format(self.contig, 795, 795)

            qname_blacklist = set()

            for pc in edited_background_af.pileup(
                    region=target, truncate=True, max_depth=self.ed.MAX_DP, stepper="all",
                    ignore_overlaps=False, ignore_orphans=False,
                    min_base_quality=self.ed.MIN_BQ, min_mapping_quality=DEFAULT_MAPQ):

                amenable_qnames = []
                for pileup_read in pc.pileups:
                    qname_alias = self.ed._get_qname_alias(pileup_read.alignment.query_name)
                    amenable_qnames.append(qname_alias)

                amenable_qname_counter = collections.Counter(amenable_qnames)
                amenable_qnames = {qname for (qname, qname_count) in amenable_qname_counter.items() if qname_count == 2}
                total_amenable_qnames = len(amenable_qnames)

                # Now we can test the method
                _, observed = self.ed._iterate_over_pileup_reads(
                    pc, self.test_variant_config_af1, observed_edit_configs, amenable_qnames,
                    total_amenable_qnames, qname_blacklist)

                edited_background_af.reset()
                self.assertEqual(expected, observed)

    def test_iterate_over_pileup_reads_edit_configs(self):
        """Tests that edit configs are appended to the edit dictionary and the expected truth frequency is returned."""

        observed_edit_configs = dict()

        with pysam.AlignmentFile(self.ed.editor_preprocessor.edit_background, "rb") as edited_background_af:

            # For purpose of testing change just the input coordinate
            target = COORD_FORMAT.format(self.contig, 795, 795)

            qname_blacklist = set()

            for pc in edited_background_af.pileup(
                    region=target, truncate=True, max_depth=self.ed.MAX_DP, stepper="all",
                    ignore_overlaps=False, ignore_orphans=False,
                    min_base_quality=self.ed.MIN_BQ, min_mapping_quality=DEFAULT_MAPQ):

                amenable_qnames = []
                for pileup_read in pc.pileups:
                    qname_alias = self.ed._get_qname_alias(pileup_read.alignment.query_name)
                    amenable_qnames.append(qname_alias)

                amenable_qname_counter = collections.Counter(amenable_qnames)
                amenable_qnames = {qname for (qname, qname_count) in amenable_qname_counter.items() if qname_count == 2}
                total_amenable_qnames = len(amenable_qnames)

                # Now we can test the method
                _, _ = self.ed._iterate_over_pileup_reads(
                    pc, self.test_variant_config_af1, observed_edit_configs, amenable_qnames, total_amenable_qnames, qname_blacklist)

                qname_alias = self.ed.qname_lookup["0000224875"]
                edit_key_r1 = ed.EDIT_KEY_TUPLE(qname=qname_alias, mate=ReadMate("R1"))
                edit_config_r1 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=795, ref="GAC", alt="ATG", read_pos=112)

                edit_key_r2 = ed.EDIT_KEY_TUPLE(qname=qname_alias, mate=ReadMate("R2"))
                edit_config_r2 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=795, ref="GAC", alt="ATG", read_pos=115)

                test_1 = len({edit_key_r1, edit_key_r2} - set(observed_edit_configs.keys())) == 0
                test_2 = observed_edit_configs[edit_key_r1] == edit_config_r1
                test_3 = observed_edit_configs[edit_key_r2] == edit_config_r2
                test_res = (test_1, test_2, test_3)

                edited_background_af.reset()
                self.assertTrue(all(test_res))

    def test_get_window_indices(self):
        """Tests that indices for a window about a position are properly returned."""

        expected = (2, 11)
        test_pos = 5
        test_len = 3
        observed = self.ed._get_window_indices(test_pos, test_len)
        self.assertEqual(expected, observed)

    def test_ref_matches_window_snp(self):
        """Test that a read segment about a SNP matches the corresponding reference segment."""

        mock_read = "CTGCGGG"

        observed = self.ed._ref_matches_window(
            contig=self.contig, query_seq=mock_read, query_pos=3, ref_len=1, ref_pos=804)

        self.assertTrue(observed)

    def test_ref_matches_window_snp_nomatch(self):
        """Test that a read segment about a SNP does not match the corresponding reference segment."""

        mock_read = "CTGCGAG"

        observed = self.ed._ref_matches_window(
            contig=self.contig, query_seq=mock_read, query_pos=3, ref_len=1, ref_pos=804)

        self.assertFalse(observed)

    def test_ref_matches_window_mnp(self):
        """Test that a read segment about a MNP matches the corresponding reference segment."""

        mock_read = "GTGGACGTG"

        observed = self.ed._ref_matches_window(
            contig=self.contig, query_seq=mock_read, query_pos=3, ref_len=3, ref_pos=795)

        self.assertTrue(observed)

    def test_get_edit_configs_trint_mnp_and_masking_detection(self):
        """Tests update of the edit config dictionary for a tri-nt MNP and implicitly tests BQ masking detection."""

        # Based on creation of the test dataset, specific reads are expected to be configured for editing
        # through a combination of the frequency and which reads have pre-existing errors

        # tri-nt MNP; presence of 1 editable read validates BQ masking detection
        qname_alias_0000224875 = self.ed.qname_lookup["0000224875"]
        edit_key_r1_0000224875 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000224875, mate=ReadMate("R1"))
        edit_config_r1_0000224875 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=795, ref="GAC", alt="ATG", read_pos=112)
        edit_key_r2_0000224875 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000224875, mate=ReadMate("R2"))
        edit_config_r2_0000224875 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=795, ref="GAC", alt="ATG", read_pos=115)

        test_1 = self.observed_edit_configs[edit_key_r1_0000224875] == edit_config_r1_0000224875
        test_2 = self.observed_edit_configs[edit_key_r2_0000224875] == edit_config_r2_0000224875
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_get_edit_configs_snp_ref_base_match_check(self):
        """Tests update of the edit config dictionary for a SNP and implicitly tests REF base check."""

        # SNP A (C>G) will be in 0000051899
        qname_alias_0000051899 = self.ed.qname_lookup["0000051899"]
        edit_key_r1_0000051899 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000051899, mate=ReadMate("R1"))
        edit_key_r2_0000051899 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000051899, mate=ReadMate("R2"))
        edit_config_r1_0000051899 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=804, ref="C", alt="T", read_pos=7)
        edit_config_r2_0000051899 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=804, ref="C", alt="T", read_pos=10)

        test_1 = self.observed_edit_configs[edit_key_r1_0000051899] == edit_config_r1_0000051899
        test_2 = self.observed_edit_configs[edit_key_r2_0000051899] == edit_config_r2_0000051899
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_get_edit_configs_multiple_snps(self):
        """Tests update of the edit config dictionary for a second SNP at the same position as another."""

        # SNP B (C>T) will be in 0000225689
        qname_alias_0000225689 = self.ed.qname_lookup["0000225689"]
        edit_key_r1_0000225689 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000225689, mate=ReadMate("R1"))
        edit_key_r2_0000225689 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000225689, mate=ReadMate("R2"))
        edit_config_r1_0000225689 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=804, ref="C", alt="G", read_pos=10)
        edit_config_r2_0000225689 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=804, ref="C", alt="G", read_pos=7)

        test_1 = self.observed_edit_configs[edit_key_r1_0000225689] == edit_config_r1_0000225689
        test_2 = self.observed_edit_configs[edit_key_r2_0000225689] == edit_config_r2_0000225689
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_get_edit_configs_dint_mnp_and_indel_detection(self):
        """Tests update of the edit config dictionary for a second SNP at the same position as another."""

        # di-nt MNP will be in 0000001369
        qname_alias_0000001369 = self.ed.qname_lookup["0000001369"]
        edit_key_r1_0000001369 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000001369, mate=ReadMate("R1"))
        edit_key_r2_0000001369 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000001369, mate=ReadMate("R2"))
        edit_config_r1_0000001369 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=812, ref="GG", alt="AT", read_pos=15)
        edit_config_r2_0000001369 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=812, ref="GG", alt="AT", read_pos=18)

        test_1 = self.observed_edit_configs[edit_key_r1_0000001369] == edit_config_r1_0000001369
        test_2 = self.observed_edit_configs[edit_key_r2_0000001369] == edit_config_r2_0000001369
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_get_edit_configs_del(self):
        """Tests update of the edit config dictionary for a deletion."""

        # 2-bp del will be in 0000003129
        qname_alias_0000003129 = self.ed.qname_lookup["0000003129"]
        edit_key_r1_0000003129 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000003129, mate=ReadMate("R1"))
        edit_key_r2_0000003129 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000003129, mate=ReadMate("R2"))
        edit_config_r1_0000003129 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=955, ref="ACT", alt="A", read_pos=47)
        edit_config_r2_0000003129 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=955, ref="ACT", alt="A", read_pos=47)

        test_1 = self.observed_edit_configs[edit_key_r1_0000003129] == edit_config_r1_0000003129
        test_2 = self.observed_edit_configs[edit_key_r2_0000003129] == edit_config_r2_0000003129
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_get_edit_configs_ins(self):
        """Tests update of the edit config dictionary for an insertion."""

        # 1-bp ins will be in 0000323491
        qname_alias_0000323491 = self.ed.qname_lookup["0000323491"]
        edit_key_r1_0000323491 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000323491, mate=ReadMate("R1"))
        edit_key_r2_0000323491 = ed.EDIT_KEY_TUPLE(qname=qname_alias_0000323491, mate=ReadMate("R2"))
        edit_config_r1_0000323491 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=1094, ref="T", alt="TG", read_pos=78)
        edit_config_r2_0000323491 = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=1094, ref="T", alt="TG", read_pos=75)

        test_1 = self.observed_edit_configs[edit_key_r1_0000323491] == edit_config_r1_0000323491
        test_2 = self.observed_edit_configs[edit_key_r2_0000323491] == edit_config_r2_0000323491
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_unmask_quals(self):
        """Tests for proper unmasking of qualities prior to write of reads."""

        observed = self.ed._unmask_quals([MASKED_BQ] * 10)
        self.assertTrue(MASKED_BQ not in set(observed))

    def test_edit_snp(self):
        """Tests edit of a SNP in a read object."""

        test_align_seg = copy.deepcopy(self.test_align_seg)
        ect = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=251, ref="C", alt="A", read_pos=10)

        self.ed._edit(align_seg=test_align_seg, variant=ect)

        expected_1 = "A"
        expected_2 = self.ed.VAR_TAG_DELIM.join([self.contig, str(251), "C", "A"])

        test_1 = test_align_seg.query_sequence[10] == expected_1
        test_2 = test_align_seg.get_tag(SAM_EDITED_TAG) == expected_2
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_edit_mnp(self):
        """Tests edit of a MNP in a read object."""

        test_align_seg = copy.deepcopy(self.test_align_seg)
        ect = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=252, ref="GGT", alt="AAC", read_pos=11)

        self.ed._edit(align_seg=test_align_seg, variant=ect)

        expected_1 = "AAC"
        expected_2 = self.ed.VAR_TAG_DELIM.join([self.contig, str(252), "GGT", "AAC"])

        test_1 = test_align_seg.query_sequence[11:14] == expected_1
        test_2 = test_align_seg.get_tag(SAM_EDITED_TAG) == expected_2
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_edit_ins(self):
        """Tests edit of an insertion in a read object."""

        test_align_seg = copy.deepcopy(self.test_align_seg)
        ect = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=251, ref="C", alt="A", read_pos=10)

        self.ed._edit(align_seg=test_align_seg, variant=ect)

        expected_1 = "A"
        expected_2 = self.ed.VAR_TAG_DELIM.join([self.contig, str(251), "C", "A"])

        test_1 = test_align_seg.query_sequence[10] == expected_1
        test_2 = test_align_seg.get_tag(SAM_EDITED_TAG) == expected_2
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_edit_del(self):
        """Tests edit of a deletion in a read object."""

        test_align_seg = copy.deepcopy(self.test_align_seg)
        ect = ed.EDIT_CONFIG_TUPLE(contig=self.contig, pos=251, ref="C", alt="A", read_pos=10)

        self.ed._edit(align_seg=test_align_seg, variant=ect)

        expected_1 = "A"
        expected_2 = self.ed.VAR_TAG_DELIM.join([self.contig, str(251), "C", "A"])

        test_1 = test_align_seg.query_sequence[10] == expected_1
        test_2 = test_align_seg.get_tag(SAM_EDITED_TAG) == expected_2
        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_iterate_over_reads(self):
        """Tests that multiple reads are edited."""

        self.ed._iterate_over_reads(self.observed_edit_configs)

        tests = []
        # Test that we actually edited the reads with a specific ALT at the expected read position
        with pysam.AlignmentFile(self.ed.temp_edit_bam, "rb") as af:
            for align_seg in af.fetch(until_eof=True):
                if align_seg.query_name == "0000224875":
                    if align_seg.is_read1:
                        tests.append(align_seg.query_sequence[112:115] == "ATG")
                    else:
                        tests.append(align_seg.query_sequence[115:118] == "ATG")
                if align_seg.query_name == "0000051899":
                    if align_seg.is_read1:
                        tests.append(align_seg.query_sequence[7] == "T")
                    else:
                        tests.append(align_seg.query_sequence[10] == "T")
                if align_seg.query_name == "0000225689":
                    if align_seg.is_read1:
                        tests.append(align_seg.query_sequence[10] == "G")
                    else:
                        tests.append(align_seg.query_sequence[7] == "G")
                if align_seg.query_name == "0000001369":
                    if align_seg.is_read1:
                        tests.append(align_seg.query_sequence[15:17] == "AT")
                    else:
                        tests.append(align_seg.query_sequence[18:20] == "AT")
                if align_seg.query_name == "0000003129":
                    if align_seg.is_read1:
                        tests.append(align_seg.query_sequence[47:50] == "AAC")
                    else:
                        tests.append(align_seg.query_sequence[47:50] == "AAC")
                if align_seg.query_name == "0000323491":
                    if align_seg.is_read1:
                        tests.append(align_seg.query_sequence[78:80] == "TG")
                    else:
                        tests.append(align_seg.query_sequence[75:77] == "TG")

        self.assertTrue(all(tests))

    def test_workflow(self):
        """Smoke test that the workflow runs, outputs files, and that the FASTQs are similarly ordered."""

        # Check that the output files exist and that the FASTQs have ordered qnames
        edit_bam, r1_fq, r2_fq = self.ed.workflow()
        test_1 = [os.path.exists(e) for e in (edit_bam, r1_fq, r2_fq,)]

        with gzip.open(r1_fq) as r1, gzip.open(r2_fq) as r2:
            r1_qnames = [e for e in r1 if e.decode().startswith(FASTQ_QNAME_CHAR)]
            r2_qnames = [e for e in r2 if e.decode().startswith(FASTQ_QNAME_CHAR)]
            test_2 = r1_qnames == r2_qnames

        test_res = (test_1, test_2,)
        self.assertTrue(all(test_res))

    def test_force_edit(self):
        """Test that a InvalidVariantConfig exception is raised with invalid variant frequency configurations."""

        with self.assertRaises(ed.InvalidVariantConfig):
            _ = ed.ReadEditor(
                self.test_bam, variants=self.test_vcf, ref=self.cbs_ref_copy, primers=self.test_primers,
                output_dir=self.tempdir, output_prefix="test_editor", force_edit=False)
