#!/usr/bin/env python3
"""Read editing objects."""

import abc
import collections
import logging
import os
import pysam
import random
import tempfile

from analysis.read_preprocessor import QnameVerification, ReadMasker
from analysis import seq_utils as su
from core_utils import file_utils as fu
from core_utils import vcf_utils as vu
from satmut_utils.definitions import QNAME_SORTS, DEFAULT_TEMPDIR
from scripts.run_bowtie2_aligner import workflow as align_workflow

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"


VARIANT_CONFIG_TUPLE = collections.namedtuple("VARIANT_CONFIG_TUPLE", "type, contig, pos, ref, alt, af")
EDIT_KEY_TUPLE = collections.namedtuple("EDIT_KEY_TUPLE", "qname, mate")
EDIT_CONFIG_TUPLE = collections.namedtuple("EDIT_CONFIG_TUPLE", "contig, pos, ref, alt, read_pos")

tempfile.tempdir = DEFAULT_TEMPDIR
logger = logging.getLogger(__name__)


class ReadEditorPreprocessor(object):
    """Class for preparing alignment files for editing of variants."""

    DEFAULT_RACE_LIKE = False
    DEFAULT_PRIMERS = None
    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    EDIT_INPUT_SUFFIX = "edit.input.bam"
    QNAME_INPUT_SUFFIX = "qname.sort.bam"

    def __init__(self, bam, ref, race_like=DEFAULT_RACE_LIKE, primers=DEFAULT_PRIMERS, outdir=DEFAULT_OUTDIR,
                 nthreads=DEFAULT_NTHREADS):
        r"""Constructor for ReadEditorPreprocessor.

        :param str bam: alignments to edit into
        :param str ref: samtools faidx indexed reference FASTA
        :param bool race_like: is the data produced by RACE-like (e.g. AMP) data? Default False.
        :param str | None primers: BED, GFF, or GTF file containing primers; for masking synthetic sequences for \
        accurate AF of variants under primer regions. This feature file should contain the strand of the primer. \
        Set to None for no masking.
        :param str outdir: Optional output directory. Default current directory.
        :param int nthreads: Number of threads to use for SAM/BAM operations. Default 0 (autodetect).
        """

        self.input_bam = bam
        self.ref = ref
        self.race_like = race_like
        self.primers = primers
        self.outdir = outdir
        self.nthreads = nthreads
        self.tempdir = tempfile.mkdtemp(suffix=".editor.preprocessor.tmp")

        # We consider masking synthetic primer regions to enable facile detection of variants "under" primers. In these
        # cases we want to ensure we don't edit into a read such that the variant appears to be a synthesis error
        input_masked_bam = self.input_bam
        if self.primers is not None:

            # Determine if the read names are compatible with masking
            qv = QnameVerification(bam=self.input_bam)
            sort_cmd = QNAME_SORTS[qv.format_index]

            rm = ReadMasker(in_bam=self.input_bam, feature_file=self.primers, race_like=race_like, sort_cmd=sort_cmd,
                            outdir=self.tempdir, nthreads=self.nthreads)
            rm.workflow()
            input_masked_bam = rm.out_bam

        # Note: do not intersect input alignments with the VCF as this could lead to unequal R1-R2 pairs that result
        # in unequal-length output FASTQs, which by definition would be invalid.
        am_basename = os.path.basename(self.input_bam)

        # Filter supplementary and secondary alignments
        # This must occur because if we choose the primary alignment for editing and writing, we should not try to edit
        # or write the supplement or secondary alignment
        preprocessed_input_bam = su.sam_view(
            input_masked_bam, None, "BAM", self.nthreads, F=su.SAM_FLAG_SUPPL + su.SAM_FLAG_SECONDARY)

        # We will need a coordinate- and qname- sorted BAM for the pileup and then the editing
        self.edit_background = os.path.join(self.tempdir, fu.replace_extension(am_basename, self.EDIT_INPUT_SUFFIX))
        su.sort_and_index(am=preprocessed_input_bam, output_am=self.edit_background, nthreads=nthreads)

        self.qname_bam = os.path.join(self.tempdir, fu.replace_extension(am_basename, self.QNAME_INPUT_SUFFIX))
        self.qname_sorted_bam = su.sort_bam(
            bam=self.edit_background, output_am=self.qname_bam, by_qname=True, nthreads=nthreads)

        temp_files = [preprocessed_input_bam]
        if self.primers is not None:
            temp_files.append(input_masked_bam)

        fu.safe_remove(tuple(temp_files))


class NonconfiguredVariant(NotImplementedError):
    """Exception to handle cases of missing configurations."""
    pass


class InvalidVariantConfig(NotImplementedError):
    """Exception to handle cases of invalid variant configurations."""
    pass


class ReadEditor(object):
    """Class for editing concordant variants into paired-end sequencing alignments."""

    DEFAULT_REFERENCE_DIR = "./references"
    DEFAULT_ENSEMBL_ID = None
    DEFAULT_REF = None
    DEFAULT_PRIMERS = None
    DEFAULT_OUTDIR = "."
    DEFAULT_PREFIX = None
    DEFAULT_FORCE = False
    DEFAULT_SEED = 9
    DEFAULT_NTHREADS = 0
    MAX_DP = 100000000
    MIN_BQ = 1  # omit primer-masked bases where BQ = 0
    VAR_TAG_DELIM = "_"
    DEFAULT_MAX_ERROR_RATE = 0.03
    DEFAULT_BUFFER = 6
    DEFAULT_MAX_NM = 10
    DEFAULT_MIN_BQ = 30

    TRUTH_VCF_SUFFIX = "truth.vcf"
    NORM_VCF_SUFFIX = "norm.sort.vcf"
    EDIT_BAM_SUFFIX = "edit.realign.bam"

    def __init__(self, bam, variants, ref, race_like=ReadEditorPreprocessor.DEFAULT_RACE_LIKE,
                 primers=DEFAULT_PRIMERS, output_dir=DEFAULT_OUTDIR, output_prefix=DEFAULT_PREFIX,
                 buffer=DEFAULT_BUFFER, max_nm=DEFAULT_MAX_NM, min_bq=DEFAULT_MIN_BQ, random_seed=DEFAULT_SEED,
                 force_edit=DEFAULT_FORCE, nthreads=DEFAULT_NTHREADS):
        r"""Constructor for ReadEditor.

        :param str bam: alignments to edit into.
        :param str variants: VCF/BCF specifying variants to edit; use the AF tag to specify AF, e.g. AF=0.1.
        :param str ref: reference FASTA. Default APPRIS primary annotation transcriptome.
        :param bool race_like: is the data produced by RACE-like (e.g. AMP) data? Default False.
        :param str | None primers: BED, GFF, or GTF file containing primers; for masking synthetic sequences for \
        accurate AF of variants under primer regions. This feature file should contain the strand of the primer. \
        Set to None for no masking.
        :param str output_dir: Optional output directory to store generated FASTQs and BAM. Default current directory.
        :param str | None output_prefix: Optional output prefix for the FASTQ(s) and BAM; If None, use bam basename.
        :param int buffer: buffer about the edit span (position + REF len) to ensure lack of error before editing. Default 3.
        :param int max_nm: max edit distance to consider a read/pair for simulation. Default 10.
        :param int min_bq: min base quality to consider a read pair for editing. Default 30.
        :param int random_seed: seed for random qname sampling. Default 9.
        :param bool force_edit: flag to attempt editing of variants despite a NonconfiguredVariant exception.
        :param int nthreads: Number of threads to use for SAM/BAM operations and alignment. Default 0 (autodetect) \
        for samtools operations. If 0, will pass 1 to bowtie2 --threads.
        """

        self.bam = bam
        self.variants = variants
        self.ref = ref
        self.race_like = race_like
        self.primers = primers
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.buffer = buffer
        self.max_nm = max_nm
        self.min_bq = min_bq
        self.random_seed = random_seed
        self.force_edit = force_edit
        self.nthreads = nthreads

        logger.info("Validating variant configurations.")
        self._verify_variant_freqs()

        logger.info("Left-normalizing variants, splitting any multi-allelic records, and sorting.")
        self.norm_sort_vcf = os.path.join(output_dir, fu.replace_extension(
            os.path.basename(variants), self.NORM_VCF_SUFFIX))

        norm_vcf = vu.VcfNormalizer(ref=ref, split_multiallelics=True).run_norm(in_vcf=variants)
        vu.VcfSorter(in_vcf=norm_vcf, out_vcf=self.norm_sort_vcf)
        self.variant_tbi = vu.tabix_index(self.norm_sort_vcf)
        fu.safe_remove((norm_vcf,))

        logger.info("Pre-processing input files for editing.")
        self.editor_preprocessor = ReadEditorPreprocessor(
            bam=self.bam, primers=primers, ref=ref, race_like=race_like, outdir=output_dir, nthreads=nthreads)

        logger.info("Getting variant configs.")
        self.variant_configs = self._get_variant_configs()

        self.variant_config_ids = [
            su.COORD_FORMAT.format(variant_config.contig, variant_config.pos, variant_config.pos)
            for variant_config in self.variant_configs]

        self.variant_config_id_set = set(self.variant_config_ids)

        self.qname_lookup = collections.defaultdict(int)

        random.seed(self.random_seed)

        if self.output_prefix is None:
            self.output_prefix = fu.remove_extension(os.path.basename(self.bam))

        self.out_path = os.path.join(self.output_dir, self.output_prefix)
        self.temp_edit_bam = tempfile.NamedTemporaryFile(mode="wb", suffix=".temp.edit.bam", delete=False).name
        self.output_bam = fu.add_extension(self.out_path, self.EDIT_BAM_SUFFIX)
        self.truth_vcf = fu.add_extension(self.out_path, self.TRUTH_VCF_SUFFIX)

    def _verify_variant_freqs(self):
        """Tests if the variant frequency sum exceeds 1 and if so, raises an exception.

        :raises analysis.read_editor.NonconfiguredVariant: if a variant does not have an AF tag
        :raises analysis.read_editor.InvalidVariantConfig: if the frequency sum across all positions exceeds 1.
        """

        with pysam.VariantFile(self.variants) as vf:

            af_sum = 0.0
            for var in vf.fetch():
                # Get the configurations for the variant to be edited
                if vu.VCF_AF_ID not in var.info:
                    var_id = "".join((var.contig, var.pos, var.ref, var.alts[0],))
                    raise NonconfiguredVariant("Please provide an %s tag for variant %s." % (vu.VCF_AF_ID, var_id))

                af_sum += var.info[vu.VCF_AF_ID]

            if af_sum > (1.0 - self.DEFAULT_MAX_ERROR_RATE):
                logger.warning(
                    "Sum of all variant frequencies exceeds (1 - the default max error rate of %f). Some variants "
                    "may not be edited due to preservation of errors." % self.DEFAULT_MAX_ERROR_RATE)

            if af_sum > 1.0 and not self.force_edit:
                raise InvalidVariantConfig(
                    "Sum of all variant frequencies exceeds 1. Due to sim design constraints, "
                    "all variants cannot be edited.")

            vf.reset()

    def _get_qname_alias(self, qname):
        """Maps qnames to integers for reduced memory footprint.

        :param str qname: read ID to lookup
        :return int: integer qname ID
        """

        if qname not in self.qname_lookup:
            new_id = len(self.qname_lookup) + 1
            self.qname_lookup[qname] = new_id
            return new_id

        qname_alias = self.qname_lookup[qname]
        return qname_alias

    def _get_variant_configs(self):
        """Gets the variant editing configurations from the VCF.

        :return list variant_tuples: variant tuples containing configurations for each variant
        """

        variant_tuples = list()

        # Use pysam's auto file mode detections
        with pysam.VariantFile(self.norm_sort_vcf, index_filename=self.variant_tbi) as vf:

            # Iterate over the variants in the file and store relevant coordinates
            for var in vf.fetch():

                # Get the configurations for the variant to be edited
                af_val = var.info[vu.VCF_AF_ID]

                vct = VARIANT_CONFIG_TUPLE(
                    type=vu.get_variant_type(var.ref, var.alts[0]), contig=var.contig, pos=var.pos,
                    ref=var.ref, alt=var.alts[0], af=af_val)

                variant_tuples.append(vct)

            vf.reset()

        return variant_tuples

    @staticmethod
    def _get_edit_qnames(amenable_qnames, variant_config, total_qnames):
        """Gets the set of qnames to edit based on depth and allele frequency.

        :param set amenable_qnames: qnames that can be edited
        :param collections.namedtuple variant_config: variant data
        :param int total_qnames: total qnames at the coordinate (less those arising from a primer at the column)
        :return set: qname IDs to edit
        """

        amenable_qname_len = len(amenable_qnames)

        # TODO: simulate sampling effects instead of a hard calculation for AF
        num_qnames_to_edit = int(float(total_qnames) * variant_config.af)

        # Make sure to protect ourselves in cases where the VcfSplitter has not been used prior and we are editing
        # multiple variants at the same position with somewhat high AFs
        if num_qnames_to_edit > amenable_qname_len:
            num_qnames_to_edit = amenable_qname_len

        # Ensure at least 1 qname can be edited if int(DP * AF) == 0
        num_qnames_to_edit = 1 if num_qnames_to_edit == 0 else num_qnames_to_edit

        # Could have no amenable qnames if the sum of AFs for variants at a position exceeds total_amenable_qnames
        qnames_to_edit = set()
        if amenable_qname_len > 0 and num_qnames_to_edit > 0:
            qnames_to_edit = set(random.sample(amenable_qnames, num_qnames_to_edit))

        return qnames_to_edit

    def _iterate_over_pileup_reads(self, pileup_column, variant_config, edit_configs, amenable_qnames,
                                   total_qnames, qname_blacklist):
        """Iterates over reads at a single column to find amenable qnames for editing.

        :param pysam.PileupColumn pileup_column: iterator over PileupRead objects
        :param read_editor.VARIANT_CONFIG_TUPLE variant_config: config for the variant
        :param dict edit_configs: dict containing editing data
        :param set amenable_qnames: amenable qnames at the coordinate
        :param int total_qnames: number total qnames at the coordinate
        :param set qname_blacklist: blacklist of qnames to avoid, as they have already been configured for editing
        :return float: expected allele frequency of the variant specified by variant_config
        """

        qnames_to_edit = self._get_edit_qnames(amenable_qnames, variant_config, total_qnames)

        # Add qnames that are configured to be edited to the blacklist, so we never edit them again
        # this is mostly required for InDel generation which upon edit, changes the edit position of the edit_config
        qname_blacklist |= qnames_to_edit

        # Once we edit a read, do not allow editing on it again at the same column so we don't glob variants
        amenable_qnames -= qnames_to_edit

        # Now we have our set of qnames to edit, so store these as edit_configs
        for pileup_read in pileup_column.pileups:

            qname_alias = self._get_qname_alias(pileup_read.alignment.query_name)

            if qname_alias in qnames_to_edit:

                edit_key = EDIT_KEY_TUPLE(qname=qname_alias, mate=su.ReadMate(pileup_read.alignment.is_read1))

                edit_config = EDIT_CONFIG_TUPLE(
                    contig=variant_config.contig, pos=variant_config.pos,
                    ref=variant_config.ref, alt=variant_config.alt, read_pos=pileup_read.query_position)

                edit_configs[edit_key] = edit_config

        # Determine the expected frequency for writing the truth VCF
        expected_cao = len(qnames_to_edit)
        expected_caf = float(expected_cao / total_qnames)
        return expected_cao, expected_caf

    @staticmethod
    def _get_concordant_qnames(all_qnames, amenable_qnames):
        """Determines the concordant qnames (overlapping mate coverage) for all and amenable qname subsets.

        :param set all_qnames: all non-masked qnames at the column
        :param list amenable_qnames: all non-masked, error-free qnames at the column
        :return tuple: (int, set) total_qnames, amenable_qnames
        """

        # Determine the total concordant depth at the column; this will be multiplied by individual variant AFs
        # to compute the integer number of qnames to edit
        total_qnames = len(all_qnames)

        # Ensure both mates overlap the POS; satmut_utils call requires both mates for variant calls
        amenable_qname_counter = collections.Counter(amenable_qnames)
        amenable_qnames = {qname for (qname, qname_count) in amenable_qname_counter.items() if qname_count == 2}

        return total_qnames, amenable_qnames

    def _get_window_indices(self, pos, ref_len):
        """Generates slice indices for a window about the position.

        :param int pos: coordinate position of the variant, local to the read or relative to the reference
        :param int ref_len: span of the variant REF
        :return tuple: (int, int) min and max indices for the window
        """

        idx_min = pos - self.buffer
        # idx_max is 1 past the end index to be checked (to accomodate python range behavior)
        idx_max = pos + ref_len + self.buffer
        indices = (idx_min, idx_max)
        return indices

    def _ref_matches_window(self, contig, query_seq, query_pos, ref_len, ref_pos):
        """Determines if the read sequence matches the reference within a window about the edit position.

        :param str contig: contig reference name
        :param str query_seq: read sequence
        :param int query_pos: position of the edit relative to the read
        :param int ref_len: length of the edited variant REF field
        :param int ref_pos: 1-based coordinate position of the edited variant
        :return bool: whether or not the read sequence matches the reference within the window
        """

        query_seq_len = len(query_seq)
        read_idx_min, read_idx_max = self._get_window_indices(query_pos, ref_len)

        # Handle cases where the read position is near the termini of the read
        read_idx_min = 0 if read_idx_min < 0 else read_idx_min
        read_idx_max = query_seq_len if read_idx_max > query_seq_len else read_idx_max
        read_seq_window = query_seq[read_idx_min:read_idx_max]

        ref_idx_min, ref_idx_max = self._get_window_indices(ref_pos, ref_len)
        ref_seq_window = su.extract_seq(contig=contig, start=ref_idx_min, stop=ref_idx_max - 1, ref=self.ref)

        if read_seq_window == ref_seq_window:
            return True

        return False

    @staticmethod
    def _reads_overlap(align_seg):
        """Determines if a read overlaps its mate.

        :param pysam.AlignedSegment align_seg: read object
        :return bool: whether the read overlaps its mate or not
        """

        if align_seg.reference_start + align_seg.query_length < align_seg.next_reference_start:
            return False

        return True

    def _get_edit_configs(self):
        """Gets a dictionary of read editing configurations and writes variants to the truth VCF.

        :return collections.defaultdict: dict of lists of reads to edit
        """

        edit_configs = dict()

        with pysam.AlignmentFile(self.editor_preprocessor.edit_background, "rb") as edited_background_af, \
                pysam.VariantFile(self.variants, "r") as in_vcf:

            # Add pertinent INFO tag headers for the truth VCF
            info_ids = (
                (vu.VCF_CAO_ID, 1, "Integer", "Concordant alternate observations, i.e. observations in both mates."),
                (vu.VCF_CAF_ID, 1, "Float", "Concordant allele frequency in range (0,1). Calculated as CAO/DP.")
            )

            vu.update_header(in_vcf.header, info_ids=info_ids)

            with pysam.VariantFile(self.truth_vcf, "w", header=in_vcf.header) as truth_vcf:

                # Note this only works for single contig alignments! This is used for target arg of the pileup call
                contig = list(self.variant_configs)[0].contig
                contig_positions = {variant_config.pos for variant_config in self.variant_configs}
                contig_pos_min = min(contig_positions)
                contig_pos_max = max(contig_positions)
                target = su.COORD_FORMAT.format(contig, contig_pos_min, contig_pos_max)

                qname_blacklist = set()

                # We must use PileupColumns in iteration only
                for pc in edited_background_af.pileup(
                        region=target, truncate=True, max_depth=self.MAX_DP, stepper="all",
                        ignore_overlaps=False, ignore_orphans=False,
                        min_base_quality=self.MIN_BQ, min_mapping_quality=su.DEFAULT_MAPQ):

                    # Skip positions that don't match any variant to be edited
                    pc_coord = su.COORD_FORMAT.format(pc.reference_name, pc.reference_pos + 1, pc.reference_pos + 1)
                    if pc_coord not in self.variant_config_id_set:
                        continue

                    # Get the variant configs associated with the current coordinate
                    var_indices = {i for i, e in enumerate(self.variant_config_ids) if e == pc_coord}
                    var_configs = [e for i, e in enumerate(self.variant_configs) if i in var_indices]

                    # First calculate the total number of amenable (non-aberrant) qnames at the coordinate;
                    # this is the depth for multiplying by variant AFs to determine the integer number of qnames to edit
                    # Here apply several filters to determine amenable qnames

                    # Find the longest REF span among the variant configs at the current position
                    # Use this REF sequence for matching against the read
                    var_configs_ref_lens = [len(vc.ref) for vc in var_configs]
                    ref_len_max = max(var_configs_ref_lens)
                    ref_len_max_idx = [i for i, e in enumerate(var_configs_ref_lens) if e == ref_len_max][0]
                    vc_max_ref_len = var_configs_ref_lens[ref_len_max_idx]

                    # Keep track of the set of all fragment qnames at the position (this is the DP denominator to AF)
                    all_qnames = set()

                    # The amenable list will additionally enforce concordance after the filters below
                    amenable_qnames = []

                    for pileup_read in pc.pileups:

                        # If the current edited base intersects any synthetic primer regions, do not edit or count for DP
                        if not (pileup_read.is_del or pileup_read.is_refskip) and su.MASKED_BQ in \
                                {pileup_read.alignment.query_qualities[pileup_read.query_position]}:
                            continue

                        # Do not consider reads/pairs that exceed the edit distance threshold
                        # To mirror the variant caller, this should check both mates simultaneously and continue for the pair
                        # However, this would require significant re-factoring, so just check each read individually
                        # Note that if one mate exceeds the thresh but the other does not, fragment will contribute to DP
                        if su.get_edit_distance(pileup_read.alignment) > self.max_nm:
                            continue

                        # Do not consider reads/pairs that do not overlap
                        if not self._reads_overlap(pileup_read.alignment):
                            continue

                        # Now at this point, consider any concordant pair with/without error at the column for the
                        # total column depth,
                        qname_alias = self._get_qname_alias(pileup_read.alignment.query_name)
                        all_qnames.add(qname_alias)

                        # If the current read position has a InDel, skip for editing
                        # This does not check nearby positions; which is the purpose of _ref_matches_window
                        if pileup_read.is_del or pileup_read.is_refskip:
                            continue

                        # If the base qualities across the REF span are not all above the min BQ, do not edit
                        ref_span_quals = list(pileup_read.alignment.query_qualities[
                                              pileup_read.query_position:pileup_read.query_position + vc_max_ref_len])

                        ref_span_quals_pass = [bq for bq in ref_span_quals if bq >= self.min_bq]

                        if len(ref_span_quals_pass) < vc_max_ref_len:
                            continue

                        # If the editable read position does not match the reference sequence in a window about the variant
                        # position, do not edit; this preserves existing errors and prohibits variant conversion by phasing
                        # with nearby errors
                        if not self._ref_matches_window(
                                contig=pileup_read.alignment.reference_name, query_seq=pileup_read.alignment.query_sequence,
                                query_pos=pileup_read.query_position, ref_len=vc_max_ref_len,
                                ref_pos=pileup_read.alignment.reference_start + pileup_read.query_position + 1):

                            continue

                        amenable_qnames.append(qname_alias)

                    # Determine count and the set of amenable qnames for editing
                    total_qnames, amenable_qnames = self._get_concordant_qnames(all_qnames, amenable_qnames)

                    # At this point amenable_qnames are all the "clean" reads at the current column
                    # Make sure to not edit into reads that have been configured to contain a variant at lower coord
                    amenable_qnames -= qname_blacklist

                    if len(amenable_qnames) == 0:
                        logger.warning("No amenable qnames for editing at %s." % pc_coord)
                        continue

                    # Get the edit configs for all variants at the column; this will update the qname_blacklist and
                    # amenable_qnames internally
                    for var_config in var_configs:

                        # Update the edit_configs dict and determine the truth frequency
                        expected_cao, expected_caf = self._iterate_over_pileup_reads(
                            pc, var_config, edit_configs, amenable_qnames, total_qnames, qname_blacklist)

                        # Write the variant and its expected frequency to the truth VCF
                        self._write_variant_to_truth_vcf(truth_vcf, var_config, expected_cao, expected_caf)

            return edit_configs

    @staticmethod
    def _unmask_quals(quals):
        """Reassigns masked primer regions to non-zero BQs so the reads do not undergo 3' quality trimming.

        :param list quals: list of BQs to unmask
        :return list: list of BQs, all non-zero
        """

        unmasked_quals = [su.DEFAULT_MAX_BQ if bq == su.MASKED_BQ else bq for bq in quals]
        return unmasked_quals

    def _edit(self, align_seg, variant):
        """Edits a single variant in the read object.

        :param pysam.AlignedSegment align_seg: read object
        :param read_editor.EDIT_CONFIG_TUPLE variant: editing config for a single variant
        """

        editor_obj = Editor.read_factory(
            seq=list(align_seg.query_sequence), quals=list(align_seg.query_qualities),
            ref=variant.ref, alt=variant.alt, read_pos=variant.read_pos)

        new_seq, new_quals = editor_obj.edit_read()

        # Sequence always must be reassigned first as per pysam docs
        align_seg.query_sequence = "".join(new_seq)
        align_seg.query_qualities = new_quals

        # Set an alignment tag so we know which reads were edited and what variant was generated
        var_tag = self.VAR_TAG_DELIM.join([variant.contig, str(variant.pos), variant.ref, variant.alt])
        align_seg.set_tag(su.SAM_EDITED_TAG, var_tag)

    def _iterate_over_reads(self, edit_configs):
        """Iterate over the reads and edit.

        :param dict edit_configs: dict with list of variants to edit at a specific read and position."""

        with pysam.AlignmentFile(self.editor_preprocessor.qname_sorted_bam, "rb") as in_af, \
                pysam.AlignmentFile(self.temp_edit_bam, mode="wb", header=in_af.header) as out_af:

            # Iterate over the qname-sorted reads, that way we write the BAM in that order for FASTQ conversion
            for align_seg in in_af.fetch(until_eof=True):

                if align_seg.is_unmapped:
                    out_af.write(align_seg)
                    continue

                qname_alias = self._get_qname_alias(align_seg.query_name)
                candidate_key = EDIT_KEY_TUPLE(qname=qname_alias, mate=su.ReadMate(align_seg.is_read1))

                # Make sure to unmask BQs if they were masked to begin with
                align_seg.query_qualities = self._unmask_quals(align_seg.query_qualities)

                if candidate_key in edit_configs:
                    self._edit(align_seg, edit_configs[candidate_key])

                # Unset the CIGAR string prior to writing as this interferes with proper object construction
                align_seg.cigarstring = None
                out_af.write(align_seg)

    def _write_fastqs(self):
        r"""Writes and gzips FASTQs for re-alignment.

        :return tuple: (str, str) paths of the R1 and R2 FASTQ files
        """

        r1_fastq, r2_fastq = su.bam_to_fastq(self.temp_edit_bam, self.out_path, True, self.nthreads)
        zipped_r1_fastq = fu.gzip_file(r1_fastq, force=True)
        zipped_r2_fastq = fu.gzip_file(r2_fastq, force=True)
        return zipped_r1_fastq, zipped_r2_fastq

    @staticmethod
    def _write_variant_to_truth_vcf(truth_vcf_fh, vc, expected_cao, expected_caf):
        """Writes a variant to the truth VCF containing the expected frequencies for input variants.

        :param pysam.VariantFile truth_vcf_fh: truth VCF VariantFile object
        :param analysis.read_editor.VARIANT_CONFIG_TUPLE vc: namedtuple containing variant info
        :param int expected_cao: expected counts for the variants
        :param float expected_caf: expected concordant allele frequency of the variant
        """

        info_dict = {vu.VCF_CAO_ID: expected_cao, vu.VCF_CAF_ID: expected_caf}

        # start should be 0-based whereas the variant configs have 1-based coordinates
        new_variant_record = truth_vcf_fh.new_record(
            contig=vc.contig, start=vc.pos - 1, alleles=(vc.ref, vc.alt), info=info_dict)

        truth_vcf_fh.write(new_variant_record)

    def workflow(self):
        """Runs the ReadEditor workflow.

        :return tuple: (str, str, str) paths of the edited and realigned BAM, R1 FASTQ, R2 FASTQ
        """

        logger.info("Getting edit configs. This could take time if the number of target positions and/or the depth "
                    "across target positions is high.")
        edit_configs = self._get_edit_configs()

        logger.info("Editing variants.")
        self._iterate_over_reads(edit_configs)

        logger.info("Writing and gzipping FASTQs.")
        zipped_r1_fastq, zipped_r2_fastq = self._write_fastqs()

        # We need to realign to re-generate CIGAR and MD tags and for proper visualization of alignments in browsers
        logger.info("Globally re-aligning edited reads.")
        nthreads = self.nthreads if self.nthreads != 0 else 1
        align_workflow(f1=zipped_r1_fastq, f2=zipped_r2_fastq, ref=self.ref, outdir=self.output_dir,
                       outbam=self.output_bam, local=False, nthreads=nthreads)

        # Remove temp files
        fu.safe_remove((self.editor_preprocessor.tempdir, self.temp_edit_bam,), force_remove=True)

        return self.output_bam, zipped_r1_fastq, zipped_r2_fastq


class Editor(object):
    """Base class for editing variant types in sequencing data."""

    __metaclass__ = abc.ABCMeta

    def __init__(self, seq, quals, ref, alt, read_pos):
        """Constructor for Editor.

        :param list seq: sequence
        :param list quals: qualities
        :param str ref: ref field
        :param str alt: alt field
        :param int read_pos: read position to edit into
        """

        self.seq = seq
        self.quals = quals
        self.orig_len = len(self.seq)
        self.ref = ref
        self.alt = alt
        self.read_pos = read_pos
        self.ref_stop_read_index = self.read_pos + len(ref)

    def paste_components(self, full_list, insert_list):
        """Puts together read pieces for an insert.

        :param list full_list: seq or qual list to edit
        :param list insert_list: insert elements
        :return list: new list including left and right arms separated by the insert
        """

        left_arm = full_list[:self.read_pos]
        right_arm = full_list[self.ref_stop_read_index:]
        new_list = left_arm + insert_list + right_arm
        end_list = new_list[:self.orig_len]
        return end_list

    @abc.abstractmethod
    def edit_seq(self):
        """Method for editing sequence. Must be implemented by child."""

        self.seq = self.paste_components(self.seq, list(self.alt))

    @abc.abstractmethod
    def edit_quals(self):
        """Method for editing quals. Must be implemented by child."""

        new_quals = [random.randint(su.DEFAULT_MIN_BQ, su.DEFAULT_MAX_BQ) for _ in range(len(self.alt))]
        self.quals = self.paste_components(self.quals, new_quals)

    def edit_read(self):
        """Creates a new seq and qual field of the AlignedSegment object.

        :return (list, list): new sequence and qualities
        """

        self.edit_seq()
        # Quals must be edited after the read
        self.edit_quals()
        return self.seq, self.quals

    @staticmethod
    def read_factory(seq, quals, ref, alt, read_pos):
        """Generate a read object based on the type of variant to edit.

        :param list seq: sequence
        :param list quals: qualities
        :param str ref: ref field
        :param str alt: alt field
        :param int read_pos: read position to edit into
        :return analysis.read_editor.ReadEditor: object for editing a variant
        """

        len_ref = len(ref)
        len_alt = len(alt)

        if len_ref == len_alt:
            if len_ref == 1:
                # Type SNP
                return SnpEditor(seq, quals, ref, alt, read_pos)
            else:
                # Type MNP; could be used for inversion too
                return MnpEditor(seq, quals, ref, alt, read_pos)
        elif len_ref < len_alt:
            # Type ins
            return InsEditor(seq, quals, ref, alt, read_pos)
        elif len_ref > len_alt:
            # Type del
            return DelEditor(seq, quals, ref, alt, read_pos)


class SnpEditor(Editor):
    """Class for creating SNPs."""

    __doc__ += Editor.__doc__

    def __init__(self, seq, quals, ref, alt, read_pos):
        """Constructor for SnpEditor."""

        super(SnpEditor, self).__init__(seq, quals, ref, alt, read_pos)

    def edit_seq(self):

        self.seq[self.read_pos] = self.alt

    def edit_quals(self):
        pass


class InsEditor(Editor):
    """Class for creating insertions."""

    __doc__ += Editor.__doc__

    def __init__(self, seq, quals, ref, alt, read_pos):
        """Constructor for InsEditor."""

        super(InsEditor, self).__init__(seq, quals, ref, alt, read_pos)

    def edit_seq(self):
        super(InsEditor, self).edit_seq()

    def edit_quals(self):
        super(InsEditor, self).edit_quals()


class DelEditor(Editor):
    """Class for creating deletions."""

    __doc__ += Editor.__doc__

    def __init__(self, seq, quals, ref, alt, read_pos):
        """Constructor for DelEditor."""

        super(DelEditor, self).__init__(seq, quals, ref, alt, read_pos)

    def edit_seq(self):
        super(DelEditor, self).edit_seq()

    def edit_quals(self):
        super(DelEditor, self).edit_quals()


class MnpEditor(Editor):
    """Class for creating MNPs."""
    
    __doc__ += Editor.__doc__

    def __init__(self, seq, quals, ref, alt, read_pos):
        """Constructor for MnpEditor."""

        super(MnpEditor, self).__init__(seq, quals, ref, alt, read_pos)

    def construct_mnp_insert(self):
        """Generates the alternate insert without a copy-paste mechanism (preserves native error).

        :return list: alternate insert
        """

        # Collect the indices and ALT bases at positions of differences between REF and ALT
        mm_data = []
        for i, (e1, e2) in enumerate(zip(self.ref, self.alt)):
            if e1 != e2:
                mm_data.append((i, e2))

        # Then reconstruct the ALT, preserving native error at non-changing positions
        mnp_insert = self.seq[self.read_pos:self.ref_stop_read_index]
        mnp_insert_len = len(mnp_insert)
        for i, e in mm_data:

            # We need to check if the index is out of bounds in cases where edited MNP spans right terminus of the read
            if i < mnp_insert_len:
                mnp_insert[i] = e

        return mnp_insert

    def edit_seq(self):

        mnp_insert = self.construct_mnp_insert()
        self.seq = self.paste_components(self.seq, mnp_insert)

    def edit_quals(self):
        pass
