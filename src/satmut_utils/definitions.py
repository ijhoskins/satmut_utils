#!/usr/bin/env python3
"""Collection of project definitions."""

import logging
import os

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
PROJECT_LAB = "Can_Cenik_lab"
PROJECT_AUTHOR = "Ian_Hoskins"

APPRIS_CONTIG_IDS = "appris_human_v1_actual_regions_contigs.txt"
APPRIS_TRX_IDS = "appris_human_v1_actual_regions_trx_ids.txt"
APPRIS_GENE_IDS = "appris_human_v1_actual_regions_gene_ids.txt"
APPRIS_TRX_FASTA = "appris_human_24_01_2019_selected.fa"
GENCODE_TRX_GFF = "gencode.v29.annotation.gtf"
GRCH38_FASTA = "GRCh38.fa.gz"

NEB_ADAPTER_P5 = "AGACGTGTGCTCTTCCGATCT"  # 5' tail of tile-seq R primers
NEB_ADAPTER_P5_RC = "AGATCGGAAGAGCACACGTCT"  # to be trimmed from 3' end of R1 in tile-seq with NEBNext adapters
NEB_ADAPTER_P7 = "TACACGACGCTCTTCCGATCT"  # 5' tail of tile-seq F primers
NEB_ADAPTER_P7_RC = "AGATCGGAAGAGCGTCGTGTA"   # to be trimmed from 3' end of R2 in tile-seq with NEBNext adapters

PEZY3_ATTB2_P5 = "ACCACTTTGTACAAGAAAGTTGG"  # exists only on the terminal R primer (for last tile)
PEZY3_ATTB2_P5_RC = "CCAACTTTCTTGTACAAAGTGGT"
PEZY3_ATTB1_P7 = "CAAGTTTGTACAAAAAAGTTGGC"  # exists only on the terminal F primer (for tile 1)
PEZY3_ATTB1_P7_RC = "GCCAACTTTTTTGTACAAACTTG"

PEZY3_ATTB1_EXT_P7 = "TGCTGGAATTGATTAATA"  # 5' Extension for tile_1_F primer

TILESEQ_UMI_LEN = 8
TILESEQ_NM_ALLOWANCE = 2

TILESEQ_UMI_REGEX = "(?P<umi_1>[ATCG]{%i})" % TILESEQ_UMI_LEN + \
                    "(?P<discard_1>%s){e<=%i}" % (PEZY3_ATTB1_EXT_P7, TILESEQ_NM_ALLOWANCE)

# CR=12 bp adapter duplex "common region" plus T from A tail ligation (13 bp)
AMP_CR = "AACCGCCAGGAGT"
AMP_CR_RC = "ACTCCTGGCGGTT"
AMP_GSP2_TAIL = "GTTCAGACGTGTGCTCTTCCGATCT"
AMP_GSP2_TAIL_RC = "AGATCGGAAGAGCACACGTCTGAAC"

AMP_UMI_LEN = 8
AMP_UMI_BUFFER = 6
AMP_CR_NM_ALLOWANCE = 2
AMP_UMI_NM_ALLOWANCE = 1

# The final UMI will be composed of the 8 bp MBC, and the first bases of the start site determined by the buffer
# This uses the python regex library syntax with named capture groups and fuzzy matching
AMP_UMI_REGEX = "(?P<umi_1>[ATCG]{%i})" % AMP_UMI_LEN + \
                "(?P<discard_1>%s){e<=%i}" % (AMP_CR, AMP_CR_NM_ALLOWANCE)

# For primer masking, some input read names (qnames) may need a different sort key than these two
ILLUMINA_QNAME_DELIM = ":"
ILLUMINA_QNAME_NFIELDS = 7
ILLUMINA_FORMAT_INDEX = 0
ILLUMINA_QNAME_INDEX_DELIM = " "
INT_FORMAT_INDEX = 1
ILLUMINA_QNAME_SORT = ("sort", "-t:", "-k1,1", "-k2,2n", "-k3,3", "-k4,4n", "-k5,5n", "-k6,6n", "-k7,7n")
INT_QNAME_SORT = ("sort", "-k1,1n")
QNAME_SORTS = (ILLUMINA_QNAME_SORT, INT_QNAME_SORT)

DEFAULT_MUT_SIG = "NNN"
VALID_MUT_SIGS = {"NNN", "NNK", "NNS"}
KEEP_INTERMEDIATES = False
DEFAULT_QUALITY_OFFSET = 33
PRE_V1p8_QUALITY_OFFSET = 64

LOG_FORMATTER = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

DEFAULT_TEMPDIR = os.getenv("SCRATCH", "/tmp")
