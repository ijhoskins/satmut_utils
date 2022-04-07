#!/usr/bin/env python3
"""Collection of sequence manipulation utilities."""

import aenum
import logging
import os
import pysam
import random
import subprocess
import tempfile

from core_utils.file_utils import flush_files, add_extension, remove_extension, FILE_NEWLINE
from satmut_utils.definitions import DEFAULT_TEMPDIR

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

tempfile.tempdir = DEFAULT_TEMPDIR
logger = logging.getLogger(__name__)

DEFAULT_ERROR_RATE = 0.01
DEFAULT_INDEL_RATE = 0.005
MASKED_BQ = 0
DEFAULT_MIN_BQ = 35
DEFAULT_MAX_BQ = 40
DEFAULT_MAPQ = 0
DEFAULT_READ_LEN = 150
MIN_ALIGN_LEN = 30
MAX_INDELS = 7
DEFAULT_FRAG_LEN = 250
ILLUMINA_BQ_OFFSET = 33
MIN_FRAG_LENGTH = 20
SD_FROM_MEAN_FACTOR = 0.25  # controls uniformity of fragment length and read counts for in silico read generation
SAM_FLAG_R1 = 64
SAM_FLAG_R2 = 128
SAM_FLAG_MINUS = 16
SAM_FLAG_PLUS = 32
SAM_FLAG_PAIRED = 1
SAM_FLAG_PROPER_PAIR = 2
SAM_FLAG_UNMAP = 4
SAM_FLAG_MUNMAP = 8
SAM_FLAG_SUPPL = 2048
SAM_FLAG_SECONDARY = 256
SAM_FLAG_REVERSE = 16

SAM_HD_TAG = "HD"
SAM_SQ_TAG = "SQ"
SAM_RG_TAG = "RG"
SAM_PG_TAG = "PG"
SAM_HD_SO_TAG = "SO"
SAM_HD_SO_QNAME_VAL = "queryname"

DNA_BASES = ("A", "C", "G", "T")
RNA_BASES = ("A", "C", "G", "U")
STOP_CODONS = ("TGA", "TAA", "TAG")
UNKNOWN_BASE = "N"
FASTA_HEADER_CHAR = ">"
FASTA_CONTIG_DELIM = "|"
FASTA_CONTIG_ID_FIELD = 0
FASTQ_QNAME_CHAR = "@"
FASTQ_SPACER_CHAR = "+"
FASTQ_SUFFIX = "fq"

FASTA_FILETYPE = "fa"
FASTA_INDEX_SUFFIX = "fai"
FASTA_FILETYPES = {"fa", "fasta", "FA", "FASTA"}
FASTQ_FILETYPES = {"fq", "fastq", "FQ", "FASTQ"}
GFF_FILETYPES = {"gff", "gtf", "GFF", "GTF"}
GFF_DEFAULT_EXT = "gff"

SAM_QNAME_INDEX = 0
SAM_BAM_SUFFIX = "am"
SAM_SUFFIX = "sam"
BAM_SUFFIX = "bam"
BAM_INDEX_SUFFIX = "bai"
SAM_HEADER_CHAR = "@"
SAM_HD_HEADER = "HD"
SAM_SQ_HEADER = "SQ"
SAM_CO_HEADER = "CO"
SAM_EDIT_DIST_TAG = "NM"
SAM_MD_TAG = "MD"
SAM_EDITED_TAG = "IN"
SAM_X0_TAG = "X0"
SAM_CIGAR_SOFTCLIP = "S"
SAM_CIGAR_HARDCLIP = "H"
SAM_CIGAR_INS = "I"
SAM_CIGAR_DEL = "D"
SAM_NULL_VAL = "."
SAM_MULTIMAP_TAG = "X0"

REVCOMP_TRANSL = str.maketrans("ACGTNacgtnRYSWKMBVDHryswkmbvdh", "TGCANtgcanYRSWMKVBHDyrswmkvbhd")
PYSAM_CIGARTUPLES_MATCH = 0
PYSAM_CIGARTUPLES_INS = 1
PYSAM_CIGARTUPLES_DEL = 2
PYSAM_CIGARTUPLES_SOFTCLIP = 4
PYSAM_CIGARTUPLES_HARDCLIP = 5
COORD_FORMAT = "{}:{}-{}"
COORD_FORMAT_STRAND = "{}:{}-{}:{}"
R_COMPAT_NA = "NA"

PYSAM_ALIGNED_PAIRS_QPOS_INDEX = 0
PYSAM_ALIGNED_PAIRS_RPOS_INDEX = 1
PYSAM_ALIGNED_PAIRS_REFB_INDEX = 2

MD_MATCH = "M"
MD_MISMATCH = "0"  # precedes base and alternate
MD_DEL = "^"  # precedes del length


class Strand(aenum.MultiValueEnum):
    """Enum for strand representation."""
    PLUS = "+", False  # alternate values for pysam.AlignedSegment.is_reverse
    MINUS = "-", True
    UNKNOWN = "", "."


class ReadMate(aenum.MultiValueEnum):
    """Enum for read mate representation."""
    R1 = "R1", True, "1"  # last values for pysam.AlignedSegment.is_read1
    R2 = "R2", False, "2"


class HumanContig(aenum.MultiValueEnum):
    """Enum for human genomic contig representation."""

    CHR1 = "1", "chr1", "NC_000001", "CM000663.2", 1
    CHR2 = "2", "chr2", "NC_000002", "CM000664.2", 2
    CHR3 = "3", "chr3", "NC_000003", "CM000665.2", 3
    CHR4 = "4", "chr4", "NC_000004", "CM000666.2", 4
    CHR5 = "5", "chr5", "NC_000005", "CM000667.2", 5
    CHR6 = "6", "chr6", "NC_000006", "CM000668.2", 6
    CHR7 = "7", "chr7", "NC_000007", "CM000669.2", 7
    CHR8 = "8", "chr8", "NC_000008", "CM000670.2", 8
    CHR9 = "9", "chr9", "NC_000009", "CM000671.2", 9
    CHR10 = "10", "chr10", "NC_000010", "CM000672.2", 10
    CHR11 = "11", "chr11", "NC_000011", "CM000673.2", 11
    CHR12 = "12", "chr12", "NC_000012", "CM000674.2", 12
    CHR13 = "13", "chr13", "NC_000013", "CM000675.2", 13
    CHR14 = "14", "chr14", "NC_000014", "CM000676.2", 14
    CHR15 = "15", "chr15", "NC_000015", "CM000677.2", 15
    CHR16 = "16", "chr16", "NC_000016", "CM000678.2", 16
    CHR17 = "17", "chr17", "NC_000017", "CM000679.2", 17
    CHR18 = "18", "chr18", "NC_000018", "CM000680.2", 18
    CHR19 = "19", "chr19", "NC_000019", "CM000681.2", 19
    CHR20 = "20", "chr20", "NC_000020", "CM000682.2", 20
    CHR21 = "21", "chr21", "NC_000021", "CM000683.2", 21
    CHR22 = "22", "chr22", "NC_000022", "CM000684.2", 22
    CHRX = "X", "chrX", "NC_000023", "CM000685.2", "23", 23
    CHRY = "Y", "chrY", "NC_000024", "CM000686.2", "24", 24
    CHRM = "MT", "chrM", "NC_012920", "J01415.2", "M"


def get_contig_lookup(transcriptome_reference):
    """Gets a lookup of short_name: full_name of reference transcript contigs.

    :param str transcriptome_reference: path of the transcriptome reference
    :return dict: mapping of short Ensembl transcript IDs to the full contig name
    """

    contig_lookup = {}

    grep_call = ("grep", "^>", transcriptome_reference)

    with tempfile.NamedTemporaryFile(suffix=".grep.res", mode="r+") as res_file:

        subprocess.run(grep_call, stdout=res_file)
        flush_files((res_file,))

        for line in res_file:
            full_name = line.strip(FASTA_HEADER_CHAR + FILE_NEWLINE)
            short_name = full_name.split(FASTA_CONTIG_DELIM)[FASTA_CONTIG_ID_FIELD]
            contig_lookup[short_name] = full_name

    return contig_lookup


def dna_to_rna(seq):
    """Converts Ts to Us.

    :param str seq: input sequence
    :return str: RNA sequence
    """

    dna_to_rna_trans = str.maketrans("Tt", "Uu")
    res = seq.translate(dna_to_rna_trans)
    return res


def reverse_complement(seq):
    """Reverse complements a sequence.

    :param str seq: input sequence
    :return str: reverse complemented sequence
    """

    rc = seq.translate(REVCOMP_TRANSL)[::-1]
    return rc


def extract_seq(contig, start, stop, ref):
    """Extracts a sequence from a reference given coordinates. Reference should have a samtools faidx index.

    :param str contig: name of contig/chrom
    :param int start: 1-based start coordinate
    :param int stop: stop coordinate
    :param str ref: indexed reference genome FASTA
    :return str: sequence
    """

    coord_str = COORD_FORMAT.format(*list(map(str, (contig, start, stop,))))

    # Use pysam.faidx for fast random access of seqs
    with pysam.FastaFile(ref) as fasta:
        seq = fasta.fetch(region=coord_str).upper()

    #seq_fa = pysam.faidx(ref, coord_str)
    #seq could wrap on many lines, but first is always header
    #seq = "".join([line for line in seq_fa.split(FILE_NEWLINE)[1:]])

    return seq


def sam_view(am, output_am=None, output_format="BAM", nthreads=0, *args, **kwargs):
    """samtools view an alignment file.

    :param str am: alignment file
    :param str | None output_am: optional output name
    :param str output_format: One of SAM, BAM, or CRAM
    :param int nthreads: number additional threads to use
    :param sequence args: single flags to pass to samtools view
    :param sequence kwargs: key-value pairs to pass to view
    :return str: output file
    """

    call_args = []

    file_mode = "w"
    if str(output_format).lower() == BAM_SUFFIX:
        file_mode = "w+b"
        call_args.extend(["-b"])

    outname = output_am
    if output_am is None:
        outname = tempfile.NamedTemporaryFile(file_mode, suffix="." + output_format.lower(), delete=False).name

    for f in args:
        call_args.extend(["-" + f])

    for k, v in kwargs.items():
        call_args.extend(["-{}".format(k), str(v)])

    call_args.extend(["-o", outname, "-O", output_format, "-@", str(nthreads)])

    call_args += [am]

    pysam.view(*call_args, catch_stdout=False)

    return outname


def sort_bam(bam, output_am=None, output_format="BAM", by_qname=False, nthreads=0, *args, **kwargs):
    """samtools sort an alignment file.

    :param str bam: alignment file
    :param str | None output_am: optional output name
    :param str output_format: One of SAM, BAM, or CRAM
    :param bool by_qname: sort by qname/read name? Default False. Sort by coordinate.
    :param int nthreads: number additional threads to use
    :param sequence args: single flags to pass to samtools view
    :param sequence kwargs: key-value pairs to pass to view
    :return str: output file
    :raises NotImplementedError: if an unrecognized samtools option is provided

    WARNINGS: --input-fmt-option, --output-fmt, and --output-fmt-option are currently not supported.
    """

    single_args = {"write-index"}
    single_char_kwargs = {"l", "m", "K", "t", "o", "T", "O"}
    multi_char_kwargs = {"no-PG", "reference", "verbosity"}

    file_mode = "w"
    if str(output_format).lower() == BAM_SUFFIX:
        file_mode = "w+b"

    outname = output_am
    if output_am is None:
        outname = tempfile.NamedTemporaryFile(file_mode, suffix="." + output_format.lower(), delete=False).name

    call_args = ["samtools", "sort", "-o", outname, "-O", output_format, "-@", str(nthreads)]

    if by_qname:
        call_args += ["-n"]

    for f in args:
        if f in single_args:
            call_args.extend(["--" + f])
        else:
            call_args.extend(["-" + f])

    for k, v in kwargs.items():
        if k in single_char_kwargs:
            call_args.extend(["-{}".format(k), v])
        elif k in multi_char_kwargs:
            call_args.extend(["--{}".format(k), v])
        else:
            raise NotImplementedError("The option %s is not recognized." % k)

    call_args += [bam]

    subprocess.run(call_args)

    # This caused issues in the context of running many TestCases at once
    # pysam.sort(*call_args, catch_stdout=False)

    return outname


def index_bam(bam):
    """Indexes a BAM.

    :param str bam: BAM file name
    """

    pysam.index(bam)


def sort_and_index(am, output_am=None, nthreads=0):
    """Sorts and indexes a SAM/BAM file.

    :param str am: alignment file
    :param str output_am: optional output name
    :param int nthreads: number additional threads to use
    :return str: output BAM file
    """

    outname = output_am
    if output_am is None:
        outname = tempfile.NamedTemporaryFile("w+b", suffix=".bam", delete=False).name

    sorted_bam = sort_bam(bam=am, output_am=outname, nthreads=nthreads)
    index_bam(sorted_bam)

    return sorted_bam


def cat_bams(bams, output_bam=None):
    """Concatenates BAM files.

    :param tuple bams: BAM files to merge
    :param str | None output_bam: optional output BAM name
    :return str: output BAM file
    """

    outname = output_bam
    if output_bam is None:
        outname = tempfile.NamedTemporaryFile(suffix=".cat.bam", delete=False).name

    file_list = list(bams)
    cat_args = ["samtools", "cat", "-o", outname] + file_list
    subprocess.call(cat_args)

    return outname


def bam_fixmate(bam, output_bam=None, nthreads=0):
    """Fixes the read SAM flags in a paired BAM.

    :param str bam: Paired BAM file to update
    :param str output_bam: optional output BAM name
    :param int nthreads: number additional threads to use
    :return str: output BAM file
    """

    outname = output_bam
    if output_bam is None:
        outname = tempfile.NamedTemporaryFile(suffix=".fixmate.sorted.bam", delete=False).name

    with tempfile.NamedTemporaryFile(suffix=".qname.sorted.bam") as qname_sorted, \
            tempfile.NamedTemporaryFile(suffix=".fixmate.bam") as fixmate_bam:

        sort_bam(bam=bam, output_am=qname_sorted.name, by_qname=True, nthreads=nthreads)
        flush_files((qname_sorted,))

        fixmate_args = ("-c", "-m", "-@", nthreads, qname_sorted.name, fixmate_bam.name)
        pysam.fixmate(*fixmate_args)
        flush_files((fixmate_bam,))

        sort_and_index(am=fixmate_bam.name, output_am=outname, nthreads=nthreads)

    return outname


def bam_to_fastq(bam, out_prefix=None, is_paired=True, nthreads=0, *args, **kwargs):
    """Converts BAM to FASTQ.

    :param str bam: BAM file to convert
    :param str | None out_prefix: output prefix to use for FASTQs, containing the directory path
    :param bool is_paired: does the BAM consist of paired reads? Default True
    :param int nthreads: number additional threads to use
    :param sequence args: additional flags to pass in as str, no - prefix
    :param dict kwargs: key-value args to pass, no --prefix
    :return tuple: (str, str | None) paths of the R1 and R2 (if present) FASTQ files
    """

    single_char_opts = {"o", "f", "F", "G", "s", "T", "v", "c"}
    multi_char_opts = {
        "i1", "i2", "barcode-tag", "quality-tag", "index-format", "reference", "verbosity"}

    outfile_prefix = out_prefix
    if out_prefix is None:
        outfile_prefix = remove_extension(os.path.abspath(bam))

    call_args = ["-n", "-@", str(nthreads)]
    for f in args:
        call_args.extend(["-" + f])

    for k, v in kwargs.items():
        if k in single_char_opts:
            call_args.extend(["-{}".format(k), v])
        elif k in multi_char_opts:
            call_args.extend(["--{}".format(k), v])
        else:
            raise NotImplementedError("The option %s is not recognized." % k)

    r1_fastq = add_extension(add_extension(outfile_prefix, "R1"), FASTQ_SUFFIX)
    r2_fastq = None
    call_args.extend(["-1", r1_fastq])

    if is_paired:
        r2_fastq = add_extension(add_extension(outfile_prefix, "R2"), FASTQ_SUFFIX)
        call_args.extend(["-2", r2_fastq])

    call_args.extend([bam])

    pysam.fastq(*call_args)

    return r1_fastq, r2_fastq


def get_edit_distance(align_seg):
    """Gets the edit distance of an alignment.

    :param pysam.AlignedSegment align_seg: read object
    :return int: edit distance
    :raises RuntimeError: if the read object has no NM or MD tags.
    """

    if align_seg.has_tag(SAM_EDIT_DIST_TAG):
        return align_seg.get_tag(SAM_EDIT_DIST_TAG)

    # Following requires the MD tag
    if not align_seg.has_tag(SAM_MD_TAG):
        raise RuntimeError("Cannot determine edit distance of read with no NM or MD tags.")

    nm = 0
    last_query_pos = 0
    last_ref_pos = 0

    for query_pos, ref_pos, b in align_seg.get_aligned_pairs(with_seq=True):

        # Mismatch
        if str(b).islower():
            nm += 1

        # Paradigm below is to count transitions from a InDel to a match
        # Insertion
        if last_query_pos is None and query_pos is not None:
            nm += 1

        # Deletion
        if last_ref_pos is None and ref_pos is not None:
            nm += 1

        last_query_pos = query_pos
        last_ref_pos = ref_pos

    return nm


def introduce_error(seq, prob=DEFAULT_ERROR_RATE, letters=DNA_BASES):
    """Introduce sequencing error in a sequence.

    :param str seq: sequence string
    :param float prob: proportion of bases to mutate
    :param tuple letters: available letters to mutate to
    :return str: sequence with error introduced
    """

    num_bases_to_mutate = round(len(seq) * prob)
    indices_to_mutate = set(random.sample(range(len(seq)), num_bases_to_mutate))
    letter_set = list(set(letters))

    mutated_list = []
    for i, b in enumerate(seq):

        if i in indices_to_mutate:
            new_base = random.choice(letter_set)
            while b == new_base:
                new_base = random.choice(letter_set)
            mutated_list.append(new_base)
        else:
            mutated_list.append(b)

    mutated_seq = "".join(mutated_list)
    return mutated_seq


def introduce_indels(seq, prob=DEFAULT_INDEL_RATE, letters=DNA_BASES):
    """Introduce InDels in a sequence.

    :param str seq: sequence string
    :param float prob: proportion of bases to insert into/delete
    :param tuple letters: available letters to mutate to
    :return str: sequence with InDels introduced
    """

    seq_len = len(seq)
    num_bases_to_mutate = round(seq_len * prob)
    indices_to_mutate = set(random.sample(range(seq_len), num_bases_to_mutate))
    letter_set = list(set(letters))

    # Length list (of 100 elements); controls weighting of InDel lengths
    len_list = [1] * 80 + [2] * 14 + [3] * 3 + [4] * 2 + [5] * 1

    mutated_list = []
    for i, b in enumerate(seq):

        if i in indices_to_mutate:

            # Determine a length of the InDel
            indel_len = random.choice(len_list)

            # Every other base to mutate will have equal probability of being a ins or del
            if random.random() >= 0.5:
                # Make an insertion
                mutated_list.append(b)  # first retain the original base
                ins = random.choices(letter_set, k=indel_len)
                mutated_list.extend(ins)
            else:
                # Make a deletion to the left of the chosen index to mutate;
                # this could potentially delete other induced errors, but is the simplest implementation
                mutated_list = mutated_list[:-indel_len]
                mutated_list.append(b)
        else:
            mutated_list.append(b)

    if len(mutated_list) > seq_len:
        mutated_list = mutated_list[:seq_len]

    mutated_seq = "".join(mutated_list)
    return mutated_seq


def bq_int_to_ascii(bq, offset=ILLUMINA_BQ_OFFSET):
    """Converts a integer base quality to an ASCII.

    :param int bq: integer BQ
    :param int offset: encoding offset
    :return str: ASCII BQ
    """

    return str(chr(bq + offset))


def bq_ascii_to_int(bq, offset=ILLUMINA_BQ_OFFSET):
    """Converts a integer base quality to an ASCII.

    :param str bq: ASCII BQ
    :param int offset: encoding offset
    :return int: integer BQ
    """

    return ord(bq) - offset


def fasta_to_fastq(in_fasta, out_fastq=None, min_bq=DEFAULT_MIN_BQ, max_bq=DEFAULT_MAX_BQ, max_len=DEFAULT_READ_LEN,
                   add_error_snps=True, add_error_indels=False, snp_prob=DEFAULT_ERROR_RATE,
                   indel_prob=DEFAULT_INDEL_RATE, letters=DNA_BASES):
    r"""Adds random quality scores to a FASTA to generate a FASTQ, optionally trims sequences to a max read length, \
    and adds sequencing error to simulate real sequencing reads.

    :param str in_fasta: input FASTQ
    :param str | None out_fastq: Optional output FASTQ
    :param int min_bq: Min BQ to assign
    :param int max_bq: Max BQ to assign
    :param int | None max_len: Max read length to trim reads to if larger than this value; if None, will not trim seqs
    :param bool add_error_snps: should sequencing errors in the form of SNPs be added? Default True.
    :param bool add_error_indels: should sequencing errors in the form of InDels be added? Default False.
    :param float snp_prob: proportion of bases to mutate to SNPs
    :param float indel_prob: proportion of bases to mutate to InDels
    :param tuple letters: available letters to mutate to
    :return str: name of the output FASTQ
    """

    def _create_fastq_lines(out_fh):
        """Writes a FASTQ record from a FASTA sequence.

        :param file out_fh: output file handle to write to
        """

        seq_len = len(seq)

        new_seq = seq
        if max_len is not None and seq_len > max_len:
            new_seq = seq[:max_len]

        # Add sequencing error and InDels to challenge variant callers
        if add_error_snps:
            new_seq = introduce_error(new_seq, prob=snp_prob, letters=letters)

        if add_error_indels:
            new_seq = introduce_indels(new_seq, prob=indel_prob, letters=letters)

        new_bqs = "".join([bq_int_to_ascii(random.randint(min_bq, max_bq)) for _ in range(len(new_seq))])

        out_fh.write(new_seq + FILE_NEWLINE)
        out_fh.write(FASTQ_SPACER_CHAR + FILE_NEWLINE)
        out_fh.write(new_bqs + FILE_NEWLINE)

    out_fq = out_fastq
    if out_fastq is None:
        out_fq = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False).name

    with open(in_fasta, "r") as infile, \
            open(out_fq, "w") as outfile:

        # Need to potentially concatenate multi-line FASTA sequences which makes for somewhat more complicated parsing
        seq = ""
        for (rec_num, line) in enumerate(infile):
            if line.startswith(FASTA_HEADER_CHAR):
                if rec_num == 0:
                    outfile.write(FASTQ_QNAME_CHAR + line[1:])
                else:
                    # Write the last entry's seq and quals before a new one
                    _create_fastq_lines(outfile)

                    # Reset the last seq then write the new qname
                    seq = ""
                    outfile.write(FASTQ_QNAME_CHAR + line[1:])
            else:
                seq += line.rstrip("\r\n")
        else:
            # This allows us to write single-record FASTAs and the last records of a multi-line FASTA
            _create_fastq_lines(outfile)

        return out_fq

