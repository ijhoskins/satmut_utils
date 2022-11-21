#!/usr/bin/env python3
"""Objects for read pre-processing."""

import collections
import gzip
import logging
import os
import numpy as np
import pybedtools
import pysam
import regex
import subprocess
import tempfile

import analysis.seq_utils as su
import core_utils.feature_file_utils as ffu
import core_utils.file_utils as fu
import core_utils.vcf_utils as vu
from satmut_utils.definitions import *
from scripts.run_bowtie2_aligner import workflow as baw

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

tempfile.tempdir = DEFAULT_TEMPDIR
logger = logging.getLogger(__name__)

MATE_STRAND_POS_TUPLE = collections.namedtuple("MATE_STRAND_POS_TUPLE", "mate, strand, pos, ref")
CONSENSUS_STATS_TUPLE = collections.namedtuple("CONSENSUS_STATS_TUPLE", "base, bq, nm")
N_INDEX_TUPLE = collections.namedtuple("N_INDEX_TUPLE", "start, len")


DEDUP_FLAG = False
CDEDUP_FLAG = False
UMITOOLS_UG_TAG = "UG"  # unique group ID

# UMI_DELIM is used as the delimiter for primer appending; UMI_SEP is used by umitools to add the canonical UMI
UMI_SEP = "_"
UMI_DELIM = "."


class QnameVerification(object):
    """Class for checking and validating read name formats, which may affect processing."""

    def __init__(self, fastq=None, bam=None):
        """Constructor for QnameVerification.

        :param str | None fastq: FASTQ file, may be gzipped.
        :param str | None bam: BAM file
        """

        self.fastq = fastq
        self.bam = bam
        self.illumina_regex = "^[a-zA-Z]"
        self.error_msg = """The QNAME format is incompatible with the primer masking workflow. A work-around is to \
        replace qnames (read names) with unique integers and provide."""
        self.compatible, self.format_index = self.workflow()

    def verify_qname_format(self, qname):
        """Determine if the qname is compatible with primer masking (Illumina standard format or a single integer).

        :param str qname: example read name
        :return tuple: (bool, int | None) whether or not the qname matches a compatible format, and the index of the format
        """

        qname_split = qname.split(ILLUMINA_QNAME_DELIM)
        qname_split_len = len(qname_split)

        # Illumina format has 7 fields and the first field normally has a non-numeric character
        if qname_split_len == ILLUMINA_QNAME_NFIELDS and regex.match(self.illumina_regex, qname_split[0]):
            return True, ILLUMINA_FORMAT_INDEX

        if qname_split_len == 1:
            try:
                _ = int(qname_split[0])
                return True, INT_FORMAT_INDEX
            except ValueError:
                return False, None

    def _verify_fastq(self):
        """Verifies qname format for a FASTQ.

        :return tuple: (bool, int | None) whether or not the qname matches a compatible format, and the index of the format
        :raises NotImplementedError: if the read name is not a recognizable format.
        """

        open_func = open
        file_mode = "r"
        if self.fastq.endswith(fu.GZ_EXTENSION):
            open_func = gzip.open
            file_mode = "rt"

        with open_func(self.fastq, mode=file_mode) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    compatible, format_index = self.verify_qname_format(
                        line.strip(su.SAM_HEADER_CHAR + fu.FILE_NEWLINE).split(ILLUMINA_QNAME_INDEX_DELIM)[0])

                    if not compatible:
                        fh.seek(0)
                        raise NotImplementedError(self.error_msg)

                    fh.seek(0)
                    return compatible, format_index

    def _verify_bam(self):
        """Verifies qname format for a BAM.

        :return tuple: (bool, int | None) whether or not the qname is a compatible format, and the index of the format
        :raises NotImplementedError: if the read name is not a recognizable format.
        """

        with pysam.AlignmentFile(self.bam, "rb") as af:
            for i, align_seg in enumerate(af.fetch(until_eof=True)):
                if i == 0:
                    compatible, format_index = self.verify_qname_format(align_seg.query_name)

                    if not compatible:
                        af.reset()
                        raise NotImplementedError(self.error_msg)

                    af.reset()
                    return compatible, format_index

    def workflow(self):
        """ Runs the qname verification workflow.

        :return tuple: (bool, int | None) whether or not the qname matches a compatible format, and the index of the format
        :raises NotImplementedError: if the read name is not a recognizable format.
        """

        logger.info("Started qname verification workflow.")

        if self.fastq is not None:
            compatible, format_index = self._verify_fastq()
            logger.info("Completed qname verification workflow.")
            return compatible, format_index

        if self.bam is not None:
            compatible, format_index = self._verify_bam()
            logger.info("Completed qname verification workflow.")
            return compatible, format_index

        logger.info("qname verification failed.")
        return False, None


class FastqPreprocessor(object):
    """Class for pre-processing FASTQs prior to alignment."""

    DEFAULT_ADAPTER = None
    OUTDIR = "."
    TRIM_QUALITY = 15
    OVERLAP_LEN = 8
    NCORES = 0
    NTRIMMED = 3
    MIN_LENGTH = 16
    TRIM_EXT = "trimmed.fq"
    TRIM_FLAG = False
    ADAPTER_DELIM = ","
    LOG_SUFFIX = ".cutadapt.json"

    def __init__(self, f1, f2, r1_fiveprime_adapters=DEFAULT_ADAPTER, r1_threeprime_adapters=DEFAULT_ADAPTER,
                 r2_fiveprime_adapters=DEFAULT_ADAPTER, r2_threeprime_adapters=DEFAULT_ADAPTER,
                 outdir=".", ncores=NCORES, trim_bq=TRIM_QUALITY, ntrimmed=NTRIMMED,
                 overlap_len=OVERLAP_LEN, no_trim=TRIM_FLAG):
        """Constructor for FastqPreprocessor.

        :param str f1: path of the R1 FASTQ
        :param str f2: path of the R2 FASTQ
        :param str | None r1_fiveprime_adapters: 5' adapters to be trimmed from R1. Default None.
        :param str | None r1_threeprime_adapters: 3' adapters to be trimmed from R1. Default None.
        :param str | None r2_fiveprime_adapters: 5' adapters to be trimmed from R2. Default None.
        :param str | None r2_threeprime_adapters: 3' adapters to be trimmed from R2. Default None.
        :param str outdir: Output directory to write preprocessed FASTQs to
        :param int ncores: Number of CPU cores to use in trimming. Default 0, autodetect.
        :param int trim_bq: quality score for quality trimming at the 3' end. Default 15.
        :param int ntrimmed: Max number of adapters to trim from each read. Default 3.
        :param int overlap_len: number of bases to match in read to trim. Default 8.
        :param bool no_trim: flag to turn off adapter and 3' base quality trimming. Default False.
        """

        self.f1 = f1
        self.f2 = f2
        self.r1_fiveprime_adapters = None
        self.r1_threeprime_adapters = None
        self.r2_fiveprime_adapters = None
        self.r2_threeprime_adapters = None

        if not no_trim:

            # Need to determine adapters for each read separately as we cannot assume the rev comp of each
            # adapter exists on the mate; this is common in standard Illumina library prep kits
            if r1_fiveprime_adapters is not None:
                self.r1_fiveprime_adapters = tuple(
                    [e.strip(fu.FILE_SPACE) for e in r1_fiveprime_adapters.split(self.ADAPTER_DELIM)])

            if r1_threeprime_adapters is not None:
                self.r1_threeprime_adapters = tuple(
                    [e.strip(fu.FILE_SPACE) for e in r1_threeprime_adapters.split(self.ADAPTER_DELIM)])

            if r2_fiveprime_adapters is not None:
                self.r2_fiveprime_adapters = tuple(
                    [e.strip(fu.FILE_SPACE) for e in r2_fiveprime_adapters.split(self.ADAPTER_DELIM)])

            if r2_threeprime_adapters is not None:
                self.r2_threeprime_adapters = tuple(
                    [e.strip(fu.FILE_SPACE) for e in r2_threeprime_adapters.split(self.ADAPTER_DELIM)])

        self.outdir = outdir
        self.ncores = ncores
        self.trim_bq = trim_bq
        self.ntrimmed = ntrimmed
        self.overlap_len = overlap_len
        self.no_trim = no_trim

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        self.log_file = os.path.join(outdir, os.path.basename(os.path.commonprefix((f1, f2,))) + self.LOG_SUFFIX)

        # Ensure trimmed FASTQs contain a gz extension if the input FASTQs were also gzipped,
        # as cutadapt outputs same compression as input FASTQs
        self.trimmed_f1 = os.path.join(outdir, fu.replace_extension(os.path.basename(f1), self.TRIM_EXT))
        if fu.is_gzipped(f1):
            self.trimmed_f1 = fu.add_extension(self.trimmed_f1, fu.get_extension(f1))

        self.trimmed_f2 = os.path.join(outdir, fu.replace_extension(os.path.basename(f2), self.TRIM_EXT))
        if fu.is_gzipped(f2):
            self.trimmed_f2 = fu.add_extension(self.trimmed_f2, fu.get_extension(f2))

        if not self.no_trim:
            self.workflow()
        else:
            self.trimmed_f1 = f1
            self.trimmed_f2 = f2

    def run_cutadapt(self):
        """Trims adapter sequence, poly-A sequences.

        Options:
        -a is 3' adapter
        -g is 5' adapter
        -e controls error tolerance; 10% by default
        -u/--cut enables hard clipping prior to adapter trimming
        """

        logger.info("Started cutadapt trimming.")

        common_call_args = ["cutadapt", "--quiet", "--quality-base", str(DEFAULT_QUALITY_OFFSET),
                            "--json", self.log_file, "-j", str(self.ncores),  "-n", str(self.ntrimmed),
                            "-q", str(self.trim_bq), "-m", str(self.MIN_LENGTH), "-O", str(self.overlap_len)]

        if self.r1_fiveprime_adapters is not None:
            r1_fp_args = ["-g {} ".format(e) for e in self.r1_fiveprime_adapters]
            common_call_args += r1_fp_args

        if self.r1_threeprime_adapters is not None:
            r1_tp_args = ["-a {} ".format(e) for e in self.r1_threeprime_adapters]
            common_call_args += r1_tp_args

        if self.r2_fiveprime_adapters is not None:
            r2_fp_args = ["-G {} ".format(e) for e in self.r2_fiveprime_adapters]
            common_call_args += r2_fp_args

        if self.r2_threeprime_adapters is not None:
            r2_tp_args = ["-A {} ".format(e) for e in self.r2_threeprime_adapters]
            common_call_args += r2_tp_args

        fq_args = ["-o", self.trimmed_f1, "-p", self.trimmed_f2, self.f1, self.f2]

        trim_call = common_call_args + fq_args
        subprocess.call(trim_call)

        logger.info("Completed cutadapt trimming.")

    def workflow(self):
        """Runs the FastqPreprocessor workflow."""

        logger.info("Started FASTQ preprocessor workflow.")
        self.run_cutadapt()
        logger.info("Completed FASTQ preprocessor workflow.")


class UMIExtractor(object):
    """Class for extracting UMIs from reads, and optionally appending primer tags for RACE-like (e.g. AMP) data."""

    DEFAULT_OUTDIR = "."
    UMI_FQ_SUFFIX = "umi.fq"
    R1_PRIMER_SUFFIX = ".r1.gsp2.fq"
    R2_PRIMER_SUFFIX = ".r2.gsp2.fq"
    STDERR_SUFFIX = "umitools_extract.stderr"
    PRIMER_SEQ_LEN = 16
    PRIMER_NM_ALLOW = 3  # NM allowance for string matching the primer sequence at the beginning of R2
    PRIMER_FASTA = None
    UNKNOWN_PRIMER_CHAR = "X"

    def __init__(self, r1_fastq, r2_fastq, umi_regex, primer_fasta=PRIMER_FASTA, primer_nm_allow=PRIMER_NM_ALLOW,
                 outdir=DEFAULT_OUTDIR):
        """Constructor for UMIExtractor.

        :param str r1_fastq: R1 FASTQ
        :param str r2_fastq: R2 FASTQ
        :param str umi_regex: regex for matching the UMIs (see umi_tools for docs)
        :param str | None primer_fasta: primer FASTA
        :param int primer_nm_allow: Edit distance allowance for primer matching Default 3.
        :param str outdir: Optional output directory. Default current working directory.
        """

        self.r1_fastq = r1_fastq
        self.r2_fastq = r2_fastq
        self.common_basename = os.path.basename(os.path.commonprefix((self.r1_fastq, self.r2_fastq)))
        self.umi_regex = umi_regex
        self.prepend_primer = True if primer_fasta is not None else False
        self.primer_fasta = primer_fasta
        self.primer_nm_allow = primer_nm_allow
        self.outdir = outdir

        self.fastq_basename = os.path.commonprefix([os.path.basename(r1_fastq), os.path.basename(r2_fastq)])
        self.extract_stderr = os.path.join(outdir, fu.add_extension(self.fastq_basename, self.STDERR_SUFFIX))

        self.r1_out_fastq = os.path.join(outdir, fu.replace_extension(os.path.basename(r1_fastq), self.UMI_FQ_SUFFIX))
        self.r2_out_fastq = os.path.join(outdir, fu.replace_extension(os.path.basename(r2_fastq), self.UMI_FQ_SUFFIX))

        self.primer_patterns = None
        if self.prepend_primer:
            self.primer_patterns = self._get_primer_patterns()

        self.workflow()

    def _get_primer_patterns(self):
        """Gets the GSP2 sequences for matching reads.

        :return dict: gsp2_name and sequence
        """

        with pysam.FastxFile(self.primer_fasta) as primer_fa:
            primer_dict = {rec.name: (
                len(rec.sequence), str(rec.sequence).upper(),
                regex.compile("(%s){e<=%i}" % (str(rec.sequence).upper(), self.primer_nm_allow))
            ) for rec in primer_fa}

        return primer_dict

    def get_orig_r2_primer(self, r2_seq):
        """Gets the originating primer sequence for a R2 based on fuzzy matching.

        :param str r2_seq: R2 sequence
        :return str: primer sequence
        """

        for primer_name, (primer_len, primer_seq, primer_re) in self.primer_patterns.items():
            exp_primer = r2_seq[:primer_len]

            # This will find the first primer from the 5' end
            # This logic might be improved by using coordinate information as well, to solve cases with high error
            # It was designed this way to be flexible to non-reference primers (vector expression and alignment to
            # insert-only reference).
            if primer_re.search(exp_primer.upper()):
                # We ought to return a sequence ID as this is more complex than primer names, which could differ
                # by as few as 1 character
                primer_id = primer_seq[-self.PRIMER_SEQ_LEN:]
                return primer_id
        else:
            # It'd be nice to just append None, but umitools extract
            # requires the UMIs to be the same length
            no_primer = "".join([self.UNKNOWN_PRIMER_CHAR] * self.PRIMER_SEQ_LEN)
            return no_primer

    def _append_primer_name(self):
        """Appends the originating R2 primer to the read names.

        :return tuple: output FASTQ filenames
        """

        with pysam.FastxFile(self.r1_fastq) as r1_ff, \
                pysam.FastxFile(self.r2_fastq) as r2_ff:

            zipped_reads = zip(r1_ff, r2_ff)

            with tempfile.NamedTemporaryFile(suffix=self.R1_PRIMER_SUFFIX, mode="w", delete=False) as r1_out, \
                    tempfile.NamedTemporaryFile(suffix=self.R2_PRIMER_SUFFIX, mode="w", delete=False) as r2_out:

                for r1, r2 in zipped_reads:
                    primer_seq = self.get_orig_r2_primer(r2.sequence.upper())
                    new_qname = r2.name + UMI_DELIM + primer_seq

                    # Paired FASTQs should always have the same names
                    r1.name = new_qname
                    r2.name = new_qname
                    r1_out.write(str(r1) + fu.FILE_NEWLINE)
                    r2_out.write(str(r2) + fu.FILE_NEWLINE)

                return r1_out.name, r2_out.name

    def _umitools_extract(self, r1_fastq, r2_fastq):
        """Extracts the molecular barcode from R1 and appends to the qname.

        :param str r1_fastq: R1 FASTQ
        :param str r2_fastq: R2 FASTQ
        """

        extract_call = ["umi_tools", "extract", "--extract-method=regex", "-p", self.umi_regex,
                        "-I", r1_fastq, "--read2-in=%s" % r2_fastq,
                        "-S", self.r1_out_fastq, "--read2-out=%s" % self.r2_out_fastq,
                        "-E", self.extract_stderr, "--log2stderr", "--temp-dir=%s" % DEFAULT_TEMPDIR]

        subprocess.call(extract_call)

    def workflow(self):
        """Runs the UMI extraction workflow."""

        logger.info("Started UMI extraction workflow.")

        r1_fq = self.r1_fastq
        r2_fq = self.r2_fastq

        # For RACE-like (AMP) data, potentially append the primer sequence (assumed start of R2) to each qname,
        # because multiple separate R2 start sites may have the same UMI/position through R1.
        if self.prepend_primer:
            logger.info("Finding originating primers for read pairs.")
            r1_fq, r2_fq = self._append_primer_name()

        logger.info("Extracting UMIs for %s." % self.common_basename)
        self._umitools_extract(r1_fastq=r1_fq, r2_fastq=r2_fq)

        if self.prepend_primer:
            fu.safe_remove((r1_fq, r2_fq,))

        logger.info("Completed UMI extraction workflow.")


class ReadGrouper(object):
    """Class for grouping UMIs in a BAM by addition of alignment tags."""

    DEFAULT_OUTDIR = "."
    GROUP_BAM_SUFFIX = "group.bam"
    STATS_SUFFIX = "group.txt"
    STDERR_SUFFIX = "umitools_group.stderr"
    UMI_NM_ALLOW = 1

    def __init__(self, in_bam, outdir=DEFAULT_OUTDIR):
        """Constructor for ReadGrouper.

        :param str in_bam: input BAM with extracted UMIs
        :param str outdir: Optional output directory. Default current working directory.
        """

        self.in_bam = in_bam
        self.outdir = outdir
        self.group_bam = os.path.join(outdir, fu.replace_extension(os.path.basename(in_bam), self.GROUP_BAM_SUFFIX))
        self.stderr = os.path.join(outdir, fu.replace_extension(os.path.basename(in_bam), self.STDERR_SUFFIX))
        self._workflow()

    def _umitools_group(self):
        """Runs umi_tools group."""

        # "--umi-group-tag", self.group_tag actually sets the tag containing the UMI, not the group
        # No option to set tag containing group ID
        group_call = ("umi_tools", "group", "-I", self.in_bam, "--paired", "--no-sort-output",
                      "--output-bam", "-S", self.group_bam,
                      "--edit-distance-threshold=%i" % self.UMI_NM_ALLOW,
                      "--umi-separator=%s" % UMI_SEP, "--ignore-tlen",
                      "--unpaired-reads=discard", "--unmapped-reads=discard",
                      "--multimapping-detection-method", su.SAM_MULTIMAP_TAG,
                      "-E", self.stderr, "--log2stderr",
                      "--temp-dir=%s" % DEFAULT_TEMPDIR)

        subprocess.call(group_call)

    def _workflow(self):
        """Runs the group workflow."""

        logger.info("Starting umi_tools group workflow.")
        self._umitools_group()
        logger.info("Completed umi_tools group workflow.")


class ReadDeduplicator(object):
    """Class for deduplicating reads generated from UMI-grouped BAMs without basic UMI-tools dedup."""

    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    DEDUP_BAM_SUFFIX = "dedup.bam"
    DEDUP_STDERR_SUFFIX = "umitools_dedup.stderr"

    def __init__(self, group_bam, outdir=DEFAULT_OUTDIR, nthreads=DEFAULT_NTHREADS):
        """Constructor for ReadDeduplicator.

        :param str group_bam: grouped input BAM
        :param str ref: reference FASTA
        :param str outdir: Optional output directory. Default current working directory.
        :param int nthreads: number threads for sort operations
        """

        self.group_bam = group_bam
        self.outdir = outdir
        self.nthreads = nthreads
        self.dedup_stderr = os.path.join(outdir, fu.add_extension(self.group_bam, self.DEDUP_STDERR_SUFFIX))
        self.dedup_bam = os.path.join(outdir, fu.add_extension(self.group_bam, self.DEDUP_BAM_SUFFIX))
        self._workflow()

    def _umitools_dedup(self):
        """Runs deduplication on a grouped BAM.

        :return str: name of the temp deduplicated BAM.
        """

        # We do not want singletons to remain in the output since we will not use them for variant calling anyways
        # Included unmapped reads, otherwise a mate of an unmapped read may remain to generate a singleton
        # "--multimapping-detection-method", su.SAM_MULTIMAP_TAG,
        with tempfile.NamedTemporaryFile(suffix=".dedup.bam", delete=False) as dedup_temp:

            dedup_call = ["umi_tools", "dedup", "-I", self.group_bam, "--paired", "-S", dedup_temp.name,
                          "--edit-distance-threshold=%i" % ReadGrouper.UMI_NM_ALLOW, "--umi-separator=%s" % UMI_SEP,
                          "--unpaired-reads=discard", "--unmapped-reads=discard", "-E", self.dedup_stderr,
                          "--output-stats=%s" % self.dedup_bam, "--log2stderr",
                          "--temp-dir=%s" % DEFAULT_TEMPDIR]

            subprocess.call(dedup_call)

            return dedup_temp.name

    def _workflow(self):
        """Runs the extraction workflow."""

        logger.info("Started deduplicating alignments for %s." % self.group_bam)
        dedup_temp = self._umitools_dedup()

        # We must filter the resulting deduplicated alignments so that we have no singletons
        # (as a result of one mate not mapping). --unmapped-reads=discard does not discard _pairs_
        # where one mate is unmapped.
        logger.info("Filtering deduplicating alignments for mapped pairs.")
        su.sam_view(am=dedup_temp, output_am=self.dedup_bam, F=su.SAM_FLAG_UNMAP + su.SAM_FLAG_MUNMAP,
                    nthreads=self.nthreads)

        su.index_bam(self.dedup_bam)
        logger.info("Completed deduplicating alignments for %s." % self.group_bam)


class ConsensusDeduplicatorPreprocessor(object):
    """Class for preprocessing reads prior to the ConsensusDeduplicator."""

    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    PREPROC_BAM_SUFFIX = "preprocess.bam"

    def __init__(self, group_bam, group_tag=UMITOOLS_UG_TAG, outdir=DEFAULT_OUTDIR, nthreads=DEFAULT_NTHREADS):
        """Constructor for ConsensusDeduplicatorPreprocessor.

        :param str group_bam: grouped input BAM
        :param str group_tag: BAM tag to store the group ID
        :param str outdir: Optional output directory. Default current working directory.
        :param int nthreads: number threads for sort operations
        """

        self.group_bam = group_bam
        self.group_tag = group_tag
        self.outdir = outdir
        self.nthreads = nthreads
        self.preprocess_bam = os.path.join(outdir, fu.add_extension(self.group_bam, self.PREPROC_BAM_SUFFIX))
        self._workflow()

    # The following class methods are exposed for convenience, in case the user wants any intermediate files
    @classmethod
    def update_tags(cls, qname_sorted, group_tag=UMITOOLS_UG_TAG, nthreads=DEFAULT_NTHREADS):
        """Updates R2 tags to contain the unique UMI group/network ID.

        :param str qname_sorted: filepath of a qname-sorted BAM
        :param str group_tag: BAM tag for the group ID. Default MI.
        :param int nthreads: number of threads to use for sorting operations
        :return str: name of an output BAM

        This method is needed because umi_tools group command outputs UMI ID tags to the R1s only.
        """

        r1_bam = su.sam_view(am=qname_sorted, nthreads=nthreads,
                             f=su.SAM_FLAG_R1, F=su.SAM_FLAG_UNMAP + su.SAM_FLAG_MUNMAP)

        r2_bam = su.sam_view(am=qname_sorted, nthreads=nthreads,
                             f=su.SAM_FLAG_R2, F=su.SAM_FLAG_UNMAP + su.SAM_FLAG_MUNMAP)

        # Iterate over the pairs, add the UG tags, then write out to a BAM to be sorted
        with tempfile.NamedTemporaryFile("wb", suffix=".tag.bam", delete=False) as tag_bam, \
                pysam.AlignmentFile(r1_bam, "rb") as af1, \
                pysam.AlignmentFile(r2_bam, "rb") as af2, \
                pysam.AlignmentFile(tag_bam, "wb", header=af1.header) as tag_af:

            for r1, r2 in zip(af1.fetch(until_eof=True), af2.fetch(until_eof=True)):

                r1_tag = r1.get_tag(group_tag)

                # Need to append the mate so the mates sort together with use of samtools sort -t
                r1_tag_final = "{}_{}".format(r1_tag, su.ReadMate.R1.value)
                r2_tag_final = "{}_{}".format(r1_tag, su.ReadMate.R2.value)

                r1.set_tag(group_tag, r1_tag_final)
                r2.set_tag(group_tag, r2_tag_final)

                tag_af.write(r1)
                tag_af.write(r2)

            tag_bam_name = tag_bam.name

        fu.safe_remove((r1_bam, r2_bam,))
        return tag_bam_name

    @classmethod
    def update_tags_from_grouped_bam(cls, in_bam, group_tag=UMITOOLS_UG_TAG, nthreads=ReadDeduplicator.DEFAULT_NTHREADS):
        """Updates the group ID tags with mate ID.

        :param str in_bam: grouped input BAM
        :param str group_tag: BAM tag for the group ID. Default MI.
        :param int nthreads: number of threads to use for sorting operations
        :return str: temp output BAM
        """

        logger.info("Sorting by read name and splitting into R1s and R2s.")
        qname_sorted = su.sort_bam(bam=in_bam, by_qname=True, nthreads=nthreads)

        logger.info("Updating group ID tags.")
        updated_bam = cls.update_tags(qname_sorted, group_tag, nthreads)

        fu.safe_remove((qname_sorted,))
        return updated_bam

    def _workflow(self):
        """Runs the preprocessing workflow for consensus deduplication.

        :return str: name of the grouped BAM (a BAM with UMI prepended to R1 UG tags)
        """

        logger.info("Started preprocessing workflow for consensus deduplication.")

        logger.info("Sorting by read name.")
        qname_sorted = su.sort_bam(bam=self.group_bam, by_qname=True, nthreads=self.nthreads)

        logger.info("Updating group ID tags for R2s.")
        updated_bam = self.update_tags(qname_sorted, self.group_tag, self.nthreads)

        logger.info("Sorting in preparation for consensus generation.")
        # Sort by the group ID tag (the mate was added to the group tag so no need for -n)
        su.sort_bam(bam=updated_bam, output_am=self.preprocess_bam, nthreads=self.nthreads, t=self.group_tag)

        fu.safe_remove((qname_sorted, updated_bam,))

        logger.info("Completed preprocessing workflow for consensus deduplication.")


class ConsensusDeduplicator(object):
    """Class for consensus deduplication starting from a UMI-grouped BAM."""

    DEFAULT_NTHREADS = 0
    DEFAULT_BAM = None
    DEFAULT_OUTDIR = "."
    UMI_NM_ALLOW = 1
    DEDUP_BAM_SUFFIX = ReadDeduplicator.DEDUP_BAM_SUFFIX
    DEFAULT_MAPQ = 40
    N_DUPLICATES_TAG = "ND"
    CONTIG_DEL_THRESH = 10  # candidate dels must be less than or equal to this value

    def __init__(self, in_bam, ref, group_tag=UMITOOLS_UG_TAG, outdir=DEFAULT_OUTDIR, out_bam=DEFAULT_BAM,
                 nthreads=DEFAULT_NTHREADS, contig_del_thresh=CONTIG_DEL_THRESH):
        """Constructor for ConsensusDeduplicator.

        :param str in_bam: input alignments with UMI network/group ID in alignment tag
        :param str ref: reference FASTA
        :param str group_tag: BAM tag for the group ID. Default UG.
        :param str outdir: Optional output directory. Default current working directory.
        :param str | None out_bam: Optional filepath of the output BAM.
        :param int nthreads: number of threads to use for alignment
        :param int contig_del_thresh: max deletion length for which del/N gaps in the merged R2 contig are called

        Note: a del/N gap refers to one of two cases:
        1) a true deletion in the alignment
        2) an unknown (N) segment between two or more merged R2s contributing to the R2 contig
        """

        self.in_bam = in_bam
        self.ref = ref
        self.group_tag = group_tag
        self.outdir = outdir
        self.nthreads = nthreads
        self.contig_del_thresh = contig_del_thresh

        outbam = out_bam
        if out_bam is None:
            outbam = os.path.join(outdir, fu.replace_extension(os.path.basename(in_bam), self.DEDUP_BAM_SUFFIX))

        self.out_bam = outbam

    def _extract_umi_network(self, align_seg):
        """Gets the UMI network/group tag from the read name.

        :param pysam.AlignedSegment align_seg: read object
        :return str: UMI network ID
        """

        res = str(align_seg.get_tag(self.group_tag)).split("_")[0]
        return res

    @staticmethod
    def _base_from_index(array_index):
        """Returns the DNA base for the corresponding array index.

        :param int array_index: index of the consensus array
        :return str: base call
        """

        if array_index == 0:
            return "A"
        elif array_index == 1:
            return "C"
        elif array_index == 2:
            return "G"
        elif array_index == 3:
            return "T"
        else:
            return su.UNKNOWN_BASE

    @staticmethod
    def _index_from_base(base):
        """Returns the array index for a given base.

        :param str base: one of {A, C, G, T, N)
        :return int | None: index into the consensus array, or None if base is not in set
        """

        if base == "A":
            array_idx = 0
        elif base == "C":
            array_idx = 1
        elif base == "G":
            array_idx = 2
        elif base == "T":
            array_idx = 3
        elif base == su.UNKNOWN_BASE:
            array_idx = 4
        else:
            return

        return array_idx

    def _update_array(self, array_ref, alt_base, bq):
        """Updates the consensus array for each position.

        :param numpy.array array_ref: reference to the array
        :param str alt_base: alternate base, usually in lowercase
        :param int bq: base quality
        :return if no base matches the reference
        """

        array_idx = self._index_from_base(alt_base)
        array_ref[array_idx] += 1
        array_ref[array_idx + 5] += bq

    @staticmethod
    def _init_positions(align_seg, mate_strand, consensus_dict):
        """Initializes new keys for consecutive and contiguous bases that the read spans in the reference

        :param pysam.AlignedSegment align_seg: read object
        :param collections.namedtuple mate_strand: mate and strand of the read
        :param collections.OrderedDict: consensus_dict
        """

        # 0-based positions
        ref_pos = align_seg.get_reference_positions()
        min_pos = min(ref_pos)
        max_pos = max(ref_pos)

        for i in range(min_pos, max_pos + 1):

            msp_tuple = MATE_STRAND_POS_TUPLE(
                mate=mate_strand.mate, strand=mate_strand.strand, pos=i, ref=mate_strand.ref)

            if msp_tuple not in consensus_dict:
                consensus_dict[msp_tuple] = np.zeros(shape=10, dtype=np.int32)

    def _get_missing_base_indices(self, consensus_quals):
        """Determines the indices of unknown bases between merged R2s in a R2 contig.

        :param list consensus_quals: BQs for the consensus
        :return set: indices in the quals list that should be converted to N
        """

        unknown_index_tuples = []
        unknown_indices = set()

        # First get a list of tuples indicating start index and contiguous length of None
        last_bq = 0
        len_counter = 0

        for i, bq in enumerate(consensus_quals):

            if bq is not None and last_bq is not None:
                last_bq = bq
                continue

            if last_bq is not None:
                start_index = i
                len_counter += 1
            elif bq is None and last_bq is None:
                len_counter += 1
            else:
                unknown_index_tuples.append(N_INDEX_TUPLE(start=start_index, len=len_counter))
                len_counter = 0

            last_bq = bq

        # At the end of the loop, we should have a list of tuples identifying the start and length of contiguous Ns
        # Finally iterate through this list and determine which contiguous stretches exceed the threshold
        for nit in unknown_index_tuples:
            if nit.len > self.contig_del_thresh:
                unknown_indices |= set(range(nit.start, nit.start + nit.len))

        return unknown_indices

    def _set_missing_bases(self, consensus_seq, consensus_quals):
        """Sets unknown base to region between merged R2s in a R2 contig.

        :param list consensus_seq: base calls for the consensus
        :param list consensus_quals: BQs for the consensus, with None at del or unknown positions
        :return tuple: updated consensus sequence and qualities, with None in qualities for True deletions
        """

        unknown_indices = self._get_missing_base_indices(consensus_quals)
        consensus_seq_update = []
        consensus_quals_update = []

        for i, (base, bq) in enumerate(zip(consensus_seq, consensus_quals)):

            if i in unknown_indices:
                consensus_seq_update.append(su.UNKNOWN_BASE)
                consensus_quals_update.append(su.DEFAULT_MAX_BQ)
            else:
                consensus_seq_update.append(base)
                consensus_quals_update.append(bq)

        return consensus_seq_update, consensus_quals_update

    def _update_consensus_dict(self, align_seg, consensus_dict, pos_list):
        """Updates the consensus and position dictionaries.

        :param pysam.AlignedSegment align_seg: read object
        :param collections.OrderedDict consensus_dict: ordered dict keeping aligned bases for each read-strand-position
        :param list pos_list: read start positions for duplicates
        """

        mate_strand = MATE_STRAND_POS_TUPLE(
            mate=su.ReadMate(align_seg.is_read1), strand=su.Strand(align_seg.is_reverse),
            pos=None, ref=align_seg.reference_name)

        # Keep track of the aligned start positions of the duplicates
        pos_list.append(align_seg.reference_start)

        # This is needed to deal with deletions in the first duplicate read,
        # whereas other duplicates may be full length. We need continuity in the dict keys.
        self._init_positions(align_seg, mate_strand, consensus_dict)

        # Need to extract the aligned pairs manually as get_aligned_pairs has unexpected behavior based on
        # matches_only kwarg. See https://github.com/pysam-developers/pysam/issues/357
        read_ap = align_seg.get_aligned_pairs(with_seq=True)
        for query_pos, ref_pos, ref_base in read_ap[align_seg.query_alignment_start:align_seg.query_alignment_end + 1]:

            # If we do not have an InDel to the reference, update the aligned base at each column
            # Consensus generation of deletions is handled after all duplicates aligned positions are queried
            if ref_pos is not None and query_pos is not None:

                mate_strand_pos = MATE_STRAND_POS_TUPLE(
                    mate=su.ReadMate(align_seg.is_read1), strand=su.Strand(align_seg.is_reverse),
                    pos=ref_pos, ref=align_seg.reference_name)

                # Enumerate the base for the match or mismatch
                # Consider adding the NM info: align_seg.get_tag(su.SAM_EDIT_DIST_TAG)
                self._update_array(
                    array_ref=consensus_dict[mate_strand_pos],
                    alt_base=align_seg.query_sequence[query_pos],
                    bq=align_seg.query_qualities[query_pos])

    def _get_consensus(self, consensus_key, consensus_array):
        """Determines the consensus base for all duplicates grouped by position.

        :param collections.namedtuple consensus_key: key containing the reference/position
        :param numpy.array consensus_array: array of counts for each base and sum of their BQs
        :return tuple: consensus base, BQ at position; ("", None) for dels
        """

        bases = consensus_array[:5]
        bqs = consensus_array[5:]
        sum_bases = sum(bases)

        # This enables facile handling of deletions or regions where reads do not cover a fragment
        # e.g. ----------->  ------------> (two non-overlapping reads on the same fragment)
        if sum_bases == 0:
            return "", None

        # Here we have duplicates
        if sum_bases == 2:
            if any(bases == 2):
                # If we have a consensus base call use it
                opt_base_index = bases.tolist().index(2)
            else:
                # Otherwise test if one of them matches the reference base and use it preferentially
                ref_base = su.extract_seq(
                    contig=consensus_key.ref, start=consensus_key.pos + 1, stop=consensus_key.pos + 1, ref=self.ref)

                opt_base_index = self._index_from_base(ref_base)

                # If neither matches the reference base, choose the one with the higher BQ (or random)
                if bases[opt_base_index] == 0:
                    opt_base_index = int(np.argmax(bqs))
        else:
            # For all other cases, use the mode to determine the consensus base
            # In the future also leverage BQs and REF:ALT information in a Bayesian framework
            opt_base_index = int(np.argmax(bases))

        opt_bq = int(bqs[opt_base_index] / bases[opt_base_index])
        cbase = self._base_from_index(opt_base_index)

        return cbase, opt_bq

    @staticmethod
    def _get_consensus_read_attrs(read_umi_network, curr_mate_strand):
        """Constructs a new read ID and SAM flag based on the consensus read data.

        :param str read_umi_network: UMI network ID
        :param collections.namedtuple curr_mate_strand: read mate and strand information
        :return tuple: (new qname, SAM flag)
        """

        sam_flag = su.SAM_FLAG_PAIRED + su.SAM_FLAG_PROPER_PAIR

        if curr_mate_strand.mate == su.ReadMate.R1:
            sam_flag += su.SAM_FLAG_R1
        elif curr_mate_strand.mate == su.ReadMate.R2:
            sam_flag += su.SAM_FLAG_R2

        if curr_mate_strand.strand == su.Strand.PLUS:
            sam_flag += su.SAM_FLAG_PLUS
        elif curr_mate_strand.strand == su.Strand.MINUS:
            sam_flag += su.SAM_FLAG_MINUS

        new_qname = str(read_umi_network.split("_")[0])
        res = new_qname, sam_flag
        return res

    @staticmethod
    def _construct_cigar(consensus_quals):
        """Generates a cigartuples list for read creation.

        :param list consensus_quals: BQs for the consensus, with None for del positions
        :return list: list for assignment to pysam.AlignedSegment.cigartuples attribute
        """

        # This is a somewhat messy way to generate the CIGAR string, and without creating insertions
        # There is probably a better implementation but this will suffice for now

        ctuples = []
        match_counter = 0
        del_counter = 0
        del_seen = False
        for cbq in consensus_quals:

            # Extend deletions
            if del_seen:
                if cbq is None:
                    del_counter += 1
                    continue
                else:
                    # If we no longer have deletion consensus bases, enumerate and move on
                    # Note this assumes we will always have a match after a deletion;
                    # a pretty darn safe assumption
                    ctuples.append((su.PYSAM_CIGARTUPLES_DEL, del_counter))
                    del_counter = 0
                    del_seen = False

            # Identify and initiate deletion operations
            if cbq is None and not del_seen:
                del_seen = True
                # Annotate the match as soon as we see a del
                ctuples.append((su.PYSAM_CIGARTUPLES_MATCH, match_counter))
                match_counter = 0
                del_counter += 1
                continue

            # All other match bases
            match_counter += 1

        # Once we have finished the loop, enumerate the last match segment
        ctuples.append((su.PYSAM_CIGARTUPLES_MATCH, match_counter))

        return ctuples

    def _construct_align_seg(self, read_umi_network, curr_mate_strand, start_pos,
                             consensus_seq, consensus_quals, n_duplicates):
        """Constructs a new read object for the consensus read.

        :param str read_umi_network: UMI network ID
        :param collections.namedtuple curr_mate_strand: read mate, strand, reference information
        :param int start_pos: aligned start coordinate
        :param list consensus_seq: base calls for the consensus
        :param list consensus_quals: BQs for the consensus
        :param int n_duplicates: number of duplicates contributing to the consensus
        :return pysam.AlignedSegment: consensus read object
        """

        # In future, instead of realigning:
        # Use samtools fixmate to fix SAM flags and the tlen field
        # Use samtools calmd to fix MD and NM tags, which are REQUIRED for variant calling

        qname, sam_flag = self._get_consensus_read_attrs(read_umi_network, curr_mate_strand)

        new_align_seg = pysam.AlignedSegment()
        new_align_seg.query_name = qname
        new_align_seg.flag = sam_flag

        # Unfortunately we need to construct a valid CIGAR
        # cigartuples = self._construct_cigar(consensus_quals)
        # new_align_seg.cigartuples = cigartuples

        # Convert missing bases in contig to N
        consensus_seq_update, consensus_quals_update = self._set_missing_bases(consensus_seq, consensus_quals)

        # At this point, for del operations we expect an empty string in consensus_seq and None in consensus_quals
        # Regions in R2 contigs between merged R2s that exceed the del threshold are converted to Ns with max BQ
        new_align_seg.query_sequence = "".join(consensus_seq_update)
        new_align_seg.query_qualities = [e for e in consensus_quals_update if e is not None]

        # TODO: ensure this works with all references
        new_align_seg.reference_id = 0
        # new_align_seg.reference_name = curr_mate_strand.ref
        new_align_seg.reference_start = start_pos
        new_align_seg.mapping_quality = self.DEFAULT_MAPQ
        new_align_seg.set_tag(self.N_DUPLICATES_TAG, n_duplicates)

        return new_align_seg

    @staticmethod
    def _add_consensus_data(cbase, cbq, consensus_seq, consensus_quals):
        """Adds the consensus information unless at an insertion position.

        :param str cbase: consensus base call
        :param int cbq: associated BQ
        :param list consensus_seq: consensus sequence to append to
        :param list consensus_quals: consensus quals to append to
        """

        # Deletions will have "" elements, and BQs None, for deletion ops
        # This will be used for CIGAR string re-generation
        consensus_seq.append(cbase)
        consensus_quals.append(cbq)

    def _write_consensus(self, out_af, consensus_dict, pos_list, read_umi_network):
        """Generates a consensus for each [mate x strand x UMI x position] combo.

        :param pysam.AlignmentFile out_af: output file to write consensus reads on the fly
        :param collections.OrderedDict consensus_dict: ordered dict keeping aligned bases for each read-strand-position
        :param pos_list pos_list: read start positions
        :param str read_umi_network: UMI network ID
        """

        last_mate_strand = list(consensus_dict.keys())[0]
        last_mate_strand = MATE_STRAND_POS_TUPLE(mate=last_mate_strand.mate, strand=last_mate_strand.strand,
                                                 pos=None, ref=last_mate_strand.ref)

        consensus_seq = []
        consensus_quals = []
        for k, v in consensus_dict.items():
            curr_mate_strand = MATE_STRAND_POS_TUPLE(mate=k.mate, strand=k.strand, pos=None, ref=k.ref)

            # Use this trick to determine when we have visited a new mate
            # This logic assumes no input mates align to same strand (discordant pair)
            if curr_mate_strand == last_mate_strand:

                # Generate the consensus base for the position
                consensus_base, consensus_bq = self._get_consensus(k, v)

                # Add the data and deal with deletions
                self._add_consensus_data(cbase=consensus_base, cbq=consensus_bq,
                                         consensus_seq=consensus_seq, consensus_quals=consensus_quals)

            else:
                # If we have started on a new mate/strand, write the results of the last mate-strand
                # and reset the consensus seq; this is needed for R2s following R1s with shared network/group ID
                # Our aligned start position of the consensus is the min start of all duplicates
                mate_strand_pos = min(pos_list)
                n_duplicates = len(pos_list)

                # Generate a new read object and write it
                new_align_seg = self._construct_align_seg(
                    read_umi_network, last_mate_strand, mate_strand_pos, consensus_seq, consensus_quals, n_duplicates)

                out_af.write(new_align_seg)

                # Reset the lists for the "next" UMI network (we are currently at it)
                reset_base, reset_qual = self._get_consensus(k, v)
                consensus_seq = [reset_base]
                consensus_quals = [reset_qual]

            last_mate_strand = curr_mate_strand

        # Write the second mate/strand in the dict once we've completed the loop
        mate_strand_pos = min(pos_list)
        n_duplicates = len(pos_list)

        new_align_seg = self._construct_align_seg(
            read_umi_network, last_mate_strand, mate_strand_pos, consensus_seq, consensus_quals, n_duplicates)

        out_af.write(new_align_seg)

    def _generate_consensus_reads(self):
        """Generates consensus reads from a UMI group-tag-sorted BAM input.

        :return str: name of the consensus BAM
        """

        with tempfile.NamedTemporaryFile("wb", suffix=".cdedup.bam", delete=False) as dedup_bam, \
                pysam.AlignmentFile(self.in_bam, "rb") as in_af:

            # We are adding a tag so update the header
            af_header = in_af.header.to_dict()
            af_header[su.SAM_CO_HEADER] = "%s: Number PCR duplicates." % vu.VCF_ND_ID
            new_header = pysam.AlignmentHeader.from_dict(af_header)

            with pysam.AlignmentFile(dedup_bam, "wb", header=new_header) as out_af:

                last_umi_network = "No_UMI"

                # TODO: convert this to an array
                # TODO: implement heuristic to avoid updating for non-duplicated UMIs
                consensus_dict = collections.OrderedDict()
                pos_list = []

                for i, align_seg in enumerate(in_af.fetch(until_eof=True)):

                    read_umi_network = self._extract_umi_network(align_seg)

                    if i == 0 or read_umi_network == last_umi_network:
                        # Store the per-base information for each read in a dict
                        self._update_consensus_dict(align_seg, consensus_dict, pos_list)
                    else:
                        # Generate the consensus for the last UMI network and write, then re-init dicts
                        self._write_consensus(out_af, consensus_dict, pos_list, last_umi_network)

                        # Regenerate the consensus dict and update with the current read
                        consensus_dict = collections.OrderedDict()
                        pos_list = []
                        self._update_consensus_dict(align_seg, consensus_dict, pos_list)

                    last_umi_network = read_umi_network

                # Write the last consensus read
                self._write_consensus(out_af, consensus_dict, pos_list, read_umi_network)

                return dedup_bam.name

    def _realign_consensus_reads(self, in_bam):
        """Realigns the consensus reads to re-generate CIGAR, MD tags, mate information.

        :param str in_bam: input BAM of deduplicated consensus reads
        :return str: output BAM name
        """

        # Write the reads to FASTQ
        r1_fastq, r2_fastq = su.bam_to_fastq(in_bam, None, True, self.nthreads, s=os.devnull)

        # Realign the reads
        bowtie2_nthreads = 1 if self.nthreads == 0 else self.nthreads
        baw(f1=r1_fastq, f2=r2_fastq, ref=self.ref, outbam=self.out_bam, nthreads=self.nthreads)

        fu.safe_remove((r1_fastq, r2_fastq,))

    def workflow(self):
        """Runs the ConsensusDeduplicator workflow."""

        logger.info("Started consensus read generation workflow for %s" % self.in_bam)
        consensus_bam = self._generate_consensus_reads()

        self._realign_consensus_reads(consensus_bam)
        logger.info("Completed consensus read generation workflow.")

        fu.safe_remove((consensus_bam,))


class ReadMasker(object):
    """Class for masking the synthetic primer region of reads."""

    MASKED_SUFFIX = "masked.bam"
    MATE_EXT_DELIM = "/"
    GROUPBY_QNAME_FIELD = 0
    GROUPBY_PRIMER_CONTIG_FIELD = 1
    GROUPBY_PRIMER_START_FIELD = 2
    GROUPBY_PRIMER_STOP_FIELD = 3
    GROUPBY_PRIMER_STRAND_FIELD = 4
    GROUPBY_DELIM = ","
    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    DEFAULT_RACE_LIKE = False
    DEFAULT_AMPLICON_START_BUFFER = 15
    DEFAULT_RACE_LIKE_START_BUFFER = 3

    def __init__(self, in_bam, feature_file, race_like=DEFAULT_RACE_LIKE, sort_cmd=QNAME_SORTS[ILLUMINA_FORMAT_INDEX],
                 outdir=DEFAULT_OUTDIR, nthreads=DEFAULT_NTHREADS):
        """Constructor for ReadMasker.

        :param str in_bam: BAM file to mask
        :param str feature_file: BED or GTF/GFF file of primer locations
        :param bool race_like: is the data produced by RACE-like (e.g. AMP) data? Default False.
        :param tuple sort_cmd: sort command and key dependent on the qname format
        :param str outdir: Optional output directory. Default current working directory.
        :param int nthreads: number threads to use for SAM/BAM file manipulations. Default 0 (autodetect).
        """

        self.in_bam = in_bam
        self.feature_file = feature_file
        self.race_like = race_like
        self.sort_cmd = sort_cmd
        self.nthreads = nthreads
        self.outdir = outdir
        self.out_bam = os.path.join(outdir, fu.replace_extension(os.path.basename(in_bam), self.MASKED_SUFFIX))

        self.feature_file_type = ffu.BED_FILETYPE if self.feature_file.endswith(ffu.BED_FILETYPE) else ffu.GFF_FILETYPE
        self.feature_file_nfields = len(pybedtools.BedTool(self.feature_file)[0].fields)

        self.primer_start_offset = ffu.BED_INTERSECT_WB_B_BED_START_OFFSET if self.feature_file_type == ffu.BED_FILETYPE \
            else ffu.BED_INTERSECT_WB_B_GFF_START_OFFSET

        self.primer_stop_offset = ffu.BED_INTERSECT_WB_B_BED_STOP_OFFSET if self.feature_file_type == ffu.BED_FILETYPE \
            else ffu.BED_INTERSECT_WB_B_GFF_STOP_OFFSET

        self.primer_strand_offset = ffu.BED_INTERSECT_WB_B_BED_STRAND_OFFSET if self.feature_file_type == ffu.BED_FILETYPE \
            else ffu.BED_INTERSECT_WB_B_GFF_STRAND_OFFSET

        # Store primer coordinates; use a more stringent buffer for RACE-like data as fragments may be amplified
        # by more than one R primer; use a loose buffer for amplicon data- it would be unusual for reads to start
        # within 15 nts of each other as primers are rarely tiled back to back on the same strand; thus, any reads
        # exhibiting this behavior are likely reads that have not been fully adapter-trimmed, or alignment artifacts
        start_buff = self.DEFAULT_RACE_LIKE_START_BUFFER if race_like else self.DEFAULT_AMPLICON_START_BUFFER
        self.primer_info = ffu.store_coords(feature_file=feature_file, use_name=False, start_buffer=start_buff)

    def _get_read_primer_intersection(self):
        """Applies a custom bedtools sort pipeline to get primers associated with each read.

        :return str: path of the output tempfile
        """

        intersect_cmd = ("bedtools", "intersect", "-bed", "-wa", "-wb", "-abam", self.in_bam, "-b", self.feature_file)

        # Determine which fields of the intersected BED to summarize based on the feature file type and expected number
        # of fields from the read portion (-wa) of the results
        primer_coord_fields = [ffu.BED_INTERSECT_WB_READ_NFIELDS,
                               ffu.BED_INTERSECT_WB_READ_NFIELDS + self.primer_start_offset,
                               ffu.BED_INTERSECT_WB_READ_NFIELDS + self.primer_stop_offset,
                               ffu.BED_INTERSECT_WB_READ_NFIELDS + self.primer_strand_offset]

        # Need 1-based field numbers for bedtools groupby
        primer_coord_fields = [str(e + 1) for e in primer_coord_fields]

        # Use groupby to collapse the intersecting primer coordinates for each read
        groupby_cmd = ("bedtools", "groupby", "-g", str(ffu.BED_INTERSECT_WB_B_BED_NAME_OFFSET + 1),
                       "-c", ",".join(primer_coord_fields), "-o", "collapse")

        # WARNING: determine if the number of fields in the read name can vary, as this could affect the sort
        # Often if we encounter improperly paired R1s/R2s during iteration in the VariantCaller, it is because
        # this sort is not working properly.

        # Finally run the custom intersect, groupby, and sort pipeline in a subprocess
        with tempfile.NamedTemporaryFile(mode="w", suffix=".read.primer.intersect.sort.bed", delete=False) as sort_out, \
                open(os.path.join(self.outdir, "%s.stderr.log" % self.__class__.__name__), "a") as masker_stderr:

            intersect_p = subprocess.Popen(intersect_cmd, stdout=subprocess.PIPE, stderr=masker_stderr)
            groupby_p = subprocess.Popen(groupby_cmd, stdin=intersect_p.stdout, stdout=subprocess.PIPE, stderr=masker_stderr)
            sort_p = subprocess.Popen(self.sort_cmd, stdin=groupby_p.stdout, stdout=sort_out, stderr=masker_stderr)
            sort_p.wait()
            intersect_p.stdout.close()
            groupby_p.stdout.close()
            return sort_out.name

    def _get_read_primer_associations(self, align_seg, groupby_res):
        """Determines originating primer(s) for each read, excluding primers where the alignment reads-through the primer.

        :param pysam.AlignedSegment: read object
        :param str groupby_res: bedtools groupby result line
        :return tuple: (str, set)- (qname, associated primer coord strings)
        """

        primer_assocs = set()

        groupby_res_split = groupby_res.rstrip(fu.FILE_NEWLINE).split(fu.FILE_DELIM)
        intersecting_qname = groupby_res_split[self.GROUPBY_QNAME_FIELD].split(self.MATE_EXT_DELIM)[0]

        # Split the contig again since we use groupby -o collapse instead of -o distinct
        intersecting_primer_contig = [groupby_res_split[self.GROUPBY_PRIMER_CONTIG_FIELD].split(self.GROUPBY_DELIM)[0]]
        intersecting_primer_starts = groupby_res_split[self.GROUPBY_PRIMER_START_FIELD].split(ffu.BED_GROUPBY_DELIM)
        intersecting_primer_stops = groupby_res_split[self.GROUPBY_PRIMER_STOP_FIELD].split(ffu.BED_GROUPBY_DELIM)
        intersecting_primer_strands = groupby_res_split[self.GROUPBY_PRIMER_STRAND_FIELD].split(ffu.BED_GROUPBY_DELIM)

        intersecting_primers = zip(intersecting_primer_contig * len(intersecting_primer_starts),
                                   intersecting_primer_starts, intersecting_primer_stops, intersecting_primer_strands)

        for intersecting_primer in intersecting_primers:

            primer_coord_str = su.COORD_FORMAT_STRAND.format(*intersecting_primer)
            primer_tuple = self.primer_info[primer_coord_str]

            # We should get no more than two matching primers. Primers in which the read has "read-through" should
            # not be associated, as we only want to mask the synthetic sequences underlying each amplicon.
            # AlignedSegment.reference_start is 0-based, whereas allowable_coords is a set of 1-based positions
            # Make sure we compare the right read coordinate with the primer coordinates
            if align_seg.reference_start + 1 in primer_tuple.allowable_coords or \
                    align_seg.reference_end in primer_tuple.allowable_coords:
                primer_assocs.add(primer_coord_str)

        return intersecting_qname, primer_assocs

    def _get_mask_base_indices(self, align_seg, associated_primers):
        """Finds the read indices to mask.

        :param pysam.AlignedSegment: read object
        :param set associated_primers: associated primers
        :return set: set of read indices to mask
        """

        base_indices_to_mask = []

        read_mate = su.ReadMate(align_seg.is_read1)
        read_strand = su.Strand(align_seg.is_reverse)

        # Convert BAM reference coordinates to 1-based for matching to all other 1-based coordinates
        # Need to keep the None values from any softclips for proper indexing into the read
        reference_positions = [p + 1 if isinstance(p, int) else p for p in
                               align_seg.get_reference_positions(full_length=True)]

        reference_positions_noclips = [p for p in reference_positions if p is not None]
        reference_positions_set = set(reference_positions_noclips)
        ref_pos_start = reference_positions_noclips[0]
        ref_pos_end = reference_positions_noclips[-1]

        for ap in associated_primers:

            primer_tuple = self.primer_info[ap]

            # 5' and 3' coordinates should be 1-based
            primer_fiveprime_coord = primer_tuple.start + 1
            primer_threeprime_coord = primer_tuple.stop
            if primer_tuple.strand == su.Strand.MINUS:
                primer_fiveprime_coord, primer_threeprime_coord = primer_threeprime_coord, primer_fiveprime_coord

            # One might think if read does not overlap the 3' end of the primer in question, there is no need to mask
            # However consider a short read with 3' BQ trimming that starts at a primer:
            # ---->  forward read
            # ------> forward primer
            # Handle this edge case first
            if primer_threeprime_coord not in reference_positions_set and \
                    ((primer_tuple.strand == su.Strand.PLUS and primer_fiveprime_coord == ref_pos_start)
                     or (primer_tuple.strand == su.Strand.MINUS and primer_fiveprime_coord == ref_pos_end)):

                # In this case the whole read should be masked only if it is Tile-seq or if it is a RACE-like R2
                if not self.race_like or (self.race_like and read_strand == su.ReadMate.R2):
                    base_indices_to_mask = set(range(0, align_seg.query_length))
                    return base_indices_to_mask

            # For any other primers whose 3' ends are not in the reference, there is no need to check
            if primer_threeprime_coord not in reference_positions_set:
                continue

            # index() will throw a ValueError if the index is not found, which is why we have the line above
            threeprime_read_index = reference_positions.index(primer_threeprime_coord)

            # For Tile-seq libraries, both ends of R1 and R2 should be masked
            if not self.race_like and ((read_strand == primer_tuple.strand and read_strand == su.Strand.PLUS) or
                                       (read_strand != primer_tuple.strand and read_strand == su.Strand.MINUS)):

                # Mask from start of read (0 index) to index of 3' end of primer (include stop in range)
                base_indices_to_mask += list(range(0, threeprime_read_index + 1))

            if not self.race_like and ((read_strand == primer_tuple.strand and read_strand == su.Strand.MINUS) or
                                       (read_strand != primer_tuple.strand and read_strand == su.Strand.PLUS)):

                # Mask from index of 3' end of primer to end of read (read length)
                base_indices_to_mask += list(range(threeprime_read_index, align_seg.query_length))

            # For RACE-like (e.g. AMP) chemistry make sure the 5' end of a R1 is never masked
            if self.race_like and read_mate == su.ReadMate.R1 and read_strand != primer_tuple.strand:

                if read_strand == su.Strand.PLUS and primer_fiveprime_coord == ref_pos_end:
                    base_indices_to_mask += list(range(threeprime_read_index, align_seg.query_length))

                if read_strand == su.Strand.MINUS and primer_fiveprime_coord == ref_pos_start:
                    base_indices_to_mask += list(range(0, threeprime_read_index + 1))

            # For RACE-like (e.g. AMP) chemistry make sure the 3' end of a R2 is never masked
            if self.race_like and read_mate == su.ReadMate.R2 and read_strand == primer_tuple.strand:

                if read_strand == su.Strand.PLUS and primer_fiveprime_coord == ref_pos_start:
                    base_indices_to_mask += list(range(0, threeprime_read_index + 1))

                if read_strand == su.Strand.MINUS and primer_fiveprime_coord == ref_pos_end:
                    base_indices_to_mask += list(range(threeprime_read_index, align_seg.query_length))

        return set(base_indices_to_mask)

    def _mask_reads(self, qname_sorted_bam, intersected_bed, masked_bam):
        """Iterate over the reads of the intersection and mask relevant bases.

        :param str qname_sorted_bam: input qname-sorted BAM
        :param str intersected_bed: input intersected with primers and groupedby
        :param file masked_bam: temporary BAM to write masked alignments to
        :raises RuntimeError: if the qname-sort between intersection results and alignments is not the same
        """

        # Iterate over the reads and mask the synthetic primer regions
        with pysam.AlignmentFile(qname_sorted_bam, check_sq=False) as in_af, \
                pysam.AlignmentFile(masked_bam, "wb", header=in_af.header) as out_af, \
                open(intersected_bed, "r") as intersected_bed_file:

            # Create a generator iterator to iterate over the qname-sorted BAM reads and corresponding bedtools groupby
            # results for the read (also sorted by qname)
            read_gen_iter = zip(in_af.fetch(until_eof=True), intersected_bed_file)

            for align_seg, intersected_groupby_res in read_gen_iter:

                if align_seg.is_unmapped:
                    out_af.write(align_seg)
                    continue

                # Get the associated (i.e. originating) primers for each read to inform masking regions
                qname, associated_primers = self._get_read_primer_associations(align_seg, intersected_groupby_res)

                if align_seg.query_name != qname:
                    raise RuntimeError(
                        "Inconsistent match between input BAM and primer intersection results. qname in input BAM was "
                        "%s. qname in intersected groupby result was %s." % (align_seg.query_name, qname))

                # Do the masking if primers were found associated with the read
                if len(associated_primers) > 0:
                    indices_to_mask = self._get_mask_base_indices(align_seg, associated_primers)
                    masked_quals = [
                        bq if i not in indices_to_mask else su.MASKED_BQ for i, bq in enumerate(align_seg.query_qualities)]
                    align_seg.query_qualities = masked_quals

                # Always write the read whether or not it needs masking
                out_af.write(align_seg)

    def workflow(self):
        """Runs the ReadMasker workflow."""

        logger.info("Started primer base quality masking for %s" % self.in_bam)

        logger.info("Intersecting alignments with primer features.")
        intersected_bed = self._get_read_primer_intersection()

        with tempfile.NamedTemporaryFile(suffix=".primer.intersected.bam") as primer_intersected:
            # We must make sure that all reads in the qname-sorted BAM to be masked will be found in the bedtools
            # groupby output such that the nreads is equivalent and we can iterate over both results simultaneously
            # Unfortunately this requires us to intersect with the primers again but output a BAM instead of BED
            ffu.intersect_features(ff1=self.in_bam, ff2=self.feature_file, outfile=primer_intersected.name)
            fu.flush_files((primer_intersected,))

            logger.info("Sorting BAM by qname to enable masking on read stream.")
            # samtools sort -n will sort read names numerically instead of lexicographically, and sorts R1s
            # prior to R2s
            qname_sorted_bam = su.sort_bam(primer_intersected.name, by_qname=True, nthreads=self.nthreads)

        logger.info("Masking synthetic primer regions in reads.")
        with tempfile.NamedTemporaryFile(suffix=".masked.bam") as masked_bam, \
                tempfile.NamedTemporaryFile(suffix=".primer.nonintersected.bam") as primer_nonintersected:

            self._mask_reads(qname_sorted_bam, intersected_bed, masked_bam)
            fu.flush_files((masked_bam,))

            # We need to add back in the alignments that did not intersect any primers and were not masked;
            # this is especially needed for the reads at the end of a cloned CDS, where an amplifying primer may
            # sit in the vector, and thus reads emanating from it may not have intersected CDS-specific primers.
            ffu.intersect_features(ff1=self.in_bam, ff2=self.feature_file, outfile=primer_nonintersected.name, v=True)
            fu.flush_files((primer_nonintersected,))

            # Add non-primer-intersecting reads back to the masked BAM here
            su.cat_bams(bams=(masked_bam.name, primer_nonintersected.name,), output_bam=self.out_bam)

        # Clean up the remaining tempfiles
        pybedtools.cleanup()
        fu.safe_remove(paths=(qname_sorted_bam, intersected_bed,))

        logger.info("Completed primer base quality masking.")


class VariantCallerPreprocessor(object):
    """Class for preprocessing alignments prior to variant calling."""

    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    QNAME_SUFFIX = "qname.sort.bam"
    R1_SUFFIX = "R1.call.bam"
    R2_SUFFIX = "R2.call.bam"

    def __init__(self, am, ref, output_dir=DEFAULT_OUTDIR, nthreads=DEFAULT_NTHREADS):
        r"""Constructor for VariantCallerPreprocessor.

        :param str am: SAM/BAM file to enumerate variants in
        :param str ref: path to reference FASTA used in alignment. Must be samtools faidx indexed.
        :param str targets: BED, GFF, or GTF file containing targeted regions to enumerate variants for
        :param str output_dir: Optional output directory. Default current working directory.
        :param int nthreads: number threads to use for SAM/BAM file manipulations. Default 0 (autodetect).
        """

        self.am = am
        self.ref = ref
        self.output_dir = output_dir
        self.nthreads = nthreads

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

        self.in_bam = self.am
        if self.am.endswith(su.SAM_SUFFIX) or (not os.path.exists(fu.add_extension(self.in_bam, su.BAM_INDEX_SUFFIX))):

            self.in_bam = os.path.join(self.output_dir, fu.replace_extension(
                os.path.basename(self.am), "in.bam"))

            logger.info("Converting SAM to BAM.")
            su.sort_and_index(am=self.am, output_am=self.in_bam, nthreads=self.nthreads)

        self.qname_sorted = os.path.join(self.output_dir, fu.replace_extension(
            os.path.basename(self.am), self.QNAME_SUFFIX))

        self.r1_calling_bam = os.path.join(self.output_dir, fu.replace_extension(
            os.path.basename(self.am), self.R1_SUFFIX))

        self.r2_calling_bam = os.path.join(self.output_dir, fu.replace_extension(
            os.path.basename(self.am), self.R2_SUFFIX))

        self.workflow()

    def workflow(self):
        """Preprocesses the alignments: unmapped pair filtering, qname-sorting, and splitting into R1, R2 BAMs."""

        logger.info("Started variant call preprocessing workflow.")

        logger.info("Sorting and splitting input BAM into R1 and R2.")
        su.sort_bam(bam=self.in_bam, output_am=self.qname_sorted, by_qname=True, nthreads=self.nthreads)

        su.sam_view(am=self.qname_sorted, output_am=self.r1_calling_bam, nthreads=self.nthreads,
                    f=su.SAM_FLAG_R1, F=su.SAM_FLAG_UNMAP + su.SAM_FLAG_MUNMAP)

        su.sam_view(am=self.qname_sorted, output_am=self.r2_calling_bam, nthreads=self.nthreads,
                    f=su.SAM_FLAG_R2, F=su.SAM_FLAG_UNMAP + su.SAM_FLAG_MUNMAP)

        logger.info("Completed variant call preprocessing workflow.")
