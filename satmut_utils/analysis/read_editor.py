#!/usr/bin/env/python
"""Read editing objects."""

import abc
import collections
import logging
import os
import pysam
import random
import tempfile

from analysis.read_preprocessor import ReadMasker
from analysis import seq_utils as su
from core_utils import file_utils as fu
from core_utils import vcf_utils as vu
from scripts.run_bowtie2_aligner import workflow as align_workflow

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"


VARIANT_CONFIG_TUPLE = collections.namedtuple("VARIANT_CONFIG_TUPLE", "type, contig, pos, ref, alt, af")
EDIT_KEY_TUPLE = collections.namedtuple("EDIT_KEY_TUPLE", "qname, mate")
EDIT_CONFIG_TUPLE = collections.namedtuple("EDIT_CONFIG_TUPLE", "contig, pos, ref, alt, read_pos")

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")

_logger = logging.getLogger(__name__)


class ReadEditorPreprocessor(object):
    """Class for preparing alignment files for editing of variants."""

    DEFAULT_PRIMERS = None
    DEFAULT_OUTDIR = "."
    DEFAULT_NTHREADS = 0
    EDIT_INPUT_SUFFIX = "edit.input.bam"
    QNAME_INPUT_SUFFIX = "qname.sort.bam"

    def __init__(self, am, ref, primers=DEFAULT_PRIMERS, outdir=DEFAULT_OUTDIR, nthreads=DEFAULT_NTHREADS):
        r"""Constructor for ReadEditorPreprocessor.

        :param str am: alignments to edit into
        :param str ref: samtools faidx indexed reference FASTA
        :param str | None primers: BED, GFF, or GTF file containing primers; for masking synthetic sequences for \
        accurate AF of variants under primer regions. This feature file should contain the strand of the primer. \
        Set to None for no masking.
        :param str outdir: Optional output directory. Default current directory.
        :param int nthreads: Number of threads to use for SAM/BAM operations. Default 0 (autodetect).
        """

        self.am = am
        self.primers = primers
        self.ref = ref
        self.outdir = outdir
        self.nthreads = nthreads

        input_bam = self.am
        if self.am.endswith(su.SAM_SUFFIX):
            _logger.info("Converting SAM to BAM.")
            input_bam = su.sam_view(self.am)

        # We consider masking synthetic primer regions to enable facile detection of variants "under" primers. In these
        # cases we want to ensure we don't edit into a read such that the variant appears to be a synthesis error
        self.preprocessed_input_bam = input_bam
        if self.primers is not None:

            self.preprocessed_input_bam = ReadMasker(
                in_bam=input_bam, feature_file=self.primers, outdir=outdir).out_bam

        # Note: do not intersect input alignments with the VCF as this could lead to unequal R1-R2 pairs that result
        # in unequal-length output FASTQs, which by definition would be invalid.
        am_basename = os.path.basename(am)

        # We will need a coordinate- and qname- sorted BAM for the pileup and then the editing
        self.edit_background = os.path.join(outdir, fu.replace_extension(am_basename, self.EDIT_INPUT_SUFFIX))
        su.sort_and_index(am=self.preprocessed_input_bam, output_am=self.edit_background, nthreads=nthreads)

        qname_bam = os.path.join(outdir, fu.replace_extension(am_basename, self.QNAME_INPUT_SUFFIX))
        self.qname_sorted_bam = su.sort_bam(self.edit_background, output_am=qname_bam, by_qname=True, nthreads=nthreads)


class ReadEditor(object):
    """Class for editing concordant variants into paired-end sequencing alignments."""

    DEFAULT_REFERENCE_DIR = "./references"
    DEFAULT_ENSEMBL_ID = None
    DEFAULT_REF = None
    DEFAULT_PRIMERS = None
    DEFAULT_OUTDIR = "."
    DEFAULT_PREFIX = None
    DEFAULT_REALIGN = False
    DEFAULT_FILTER = True
    DEFAULT_SEED = 9
    DEFAULT_NTHREADS = 0

    MAX_DP = 100000000
    MIN_BQ = 1  # omit primer-masked bases where BQ = 0
    MATCH_REQ_LIM = 3
    VAR_TAG_DELIM = "_"

    TRUTH_VCF_SUFFIX = "truth.vcf"
    NORM_VCF_SUFFIX = "norm.vcf"
    EDIT_BAM_SUFFIX = "edit.bam"

    def __init__(self, am, variants, ref, primers=DEFAULT_PRIMERS, output_dir=DEFAULT_OUTDIR,
                 output_prefix=DEFAULT_PREFIX, random_seed=DEFAULT_SEED, nthreads=DEFAULT_NTHREADS):
        r"""Constructor for ReadEditor.

        :param str am: alignments to edit into.
        :param str variants: VCF/BCF specifying variants to edit; use the AF tag to specify AF, e.g. AF=0.1.
        :param str ref: reference FASTA. Default APPRIS primary annotation transcriptome.
        :param str | None primers: BED, GFF, or GTF file containing primers; for masking synthetic sequences for \
        accurate AF of variants under primer regions. This feature file should contain the strand of the primer. \
        Set to None for no masking.
        :param str output_dir: Optional output directory to store generated FASTQs and BAM. Default current directory.
        :param str | None output_prefix: Optional output prefix for the FASTQ(s) and BAM
        :param int random_seed: seed for random qname sampling
        :param int nthreads: Number of threads to use for SAM/BAM operations. Default 0 (autodetect).

        Note: Reads that have been edited will be tagged with an IN alignment tag.
        """

        self.am = am
        self.variants = variants
        self.primers = primers
        self.ref = ref
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.random_seed = random_seed
        self.nthreads = nthreads

        _logger.info("Left-normalizing variants, splitting any multi-allelic records, and sorting.")
        self.norm_sort_vcf = os.path.join(output_dir, fu.replace_extension(
            os.path.basename(variants), self.NORM_VCF_SUFFIX))

        norm_vcf = vu.VcfNormalizer(ref=ref, split_multiallelics=True).run_norm(in_vcf=variants)
        vu.VcfSorter(in_vcf=norm_vcf, out_vcf=self.norm_sort_vcf)
        self.variant_tbi = vu.tabix_index(self.norm_sort_vcf)

        if self.output_prefix is None:
            self.output_prefix = fu.remove_extension(os.path.basename(am))

        self.out_path = os.path.join(self.output_dir, self.output_prefix)
        self.output_bam = fu.add_extension(self.out_path, self.EDIT_BAM_SUFFIX)
        self.truth_vcf = fu.add_extension(self.out_path, self.TRUTH_VCF_SUFFIX)

        _logger.info("Pre-processing input files for editing.")
        self.editor_preprocessor = ReadEditorPreprocessor(am=am, primers=primers, ref=ref, outdir=output_dir)

        _logger.info("Getting variant configs.")
        self.variant_configs = self._get_variant_configs()

        self.variant_config_ids = [
            su.COORD_FORMAT.format(variant_config.contig, variant_config.pos, variant_config.pos)
            for variant_config in self.variant_configs]

        self.variant_config_id_set = set(self.variant_config_ids)

        self.qname_lookup = {}

        random.seed(self.random_seed)

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
                af_val = 1.0
                if vu.VCF_AF_ID in var.info:
                    af_val = var.info[vu.VCF_AF_ID]

                vct = VARIANT_CONFIG_TUPLE(
                    type=vu.get_variant_type(var.ref, var.alts[0]), contig=var.contig, pos=var.pos,
                    ref=var.ref, alt=var.alts[0], af=af_val)

                variant_tuples.append(vct)

            vf.reset()

        return variant_tuples

    @staticmethod
    def _get_edit_qnames(amenable_qnames, variant_config, total_amenable_qnames):
        """Gets the set of qnames to edit based on depth and allele frequency.

        :param set amenable_qnames: qnames that can be edited
        :param collections.namedtuple variant_config: variant data
        :param int total_amenable_qnames: total amenable qnames at the coordinate
        :return set: qname IDs to edit
        """

        amenable_qname_len = len(amenable_qnames)

        # TODO: simulate sampling effects instead of a hard calculation for AF
        num_qnames_to_edit = int(float(total_amenable_qnames) * variant_config.af)

        # Make sure to protect ourselves in cases where the VcfSplitter has not been used prior and we are editing
        # multiple variants at the same position with somewhat high AFs
        if num_qnames_to_edit > amenable_qname_len:
            num_qnames_to_edit = int(float(amenable_qname_len) * variant_config.af)

        # Ensure at least 1 qname can be edited if int(DP * AF) == 0
        num_qnames_to_edit = 1 if num_qnames_to_edit == 0 else num_qnames_to_edit

        # Could have no amenable qnames if the sum of AFs for variants at a position exceeds total_amenable_qnames
        if amenable_qname_len != 0:
            qnames_to_edit = set(random.sample(amenable_qnames, num_qnames_to_edit))
        else:
            qnames_to_edit = set()

        return qnames_to_edit

    def _iterate_over_pileup_reads(self, pileup_column, variant_config, edit_configs, amenable_qnames,
                                   total_amenable_qnames):
        """Iterates over reads at a single column to find amenable qnames for editing.

        :param pysam.PileupColumn pileup_column: iterator over PileupRead objects
        :param read_editor.VARIANT_CONFIG_TUPLE variant_config: config for the variant
        :param dict edit_configs: dict containing editing data
        :param set amenable_qnames: amenable qnames as the coordinate
        :param int total_amenable_qnames: number amenable qnames as the coordinate
        :return float: expected allele frequency of the variant specified by variant_config
        """

        qnames_to_edit = self._get_edit_qnames(amenable_qnames, variant_config, total_amenable_qnames)

        # Once we edit a read, do not allow editing on it again so we don't glob variants;
        # also required for InDel generation which upon edit, changes the edit position of the edit_config
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
        expected_af = len(qnames_to_edit) / total_amenable_qnames
        return expected_af

    def _get_edit_configs(self):
        """Gets a dictionary of read variant editing configurations. Editing should be deterministic.

        :return collections.defaultdict: dict of lists of reads to edit
        """

        edit_configs = dict()
        random.seed(self.random_seed)

        with pysam.AlignmentFile(self.editor_preprocessor.edit_background, "rb") as edited_background_af, \
                pysam.VariantFile(self.variants, "r") as in_vcf, \
                pysam.VariantFile(self.truth_vcf, "w", header=in_vcf.header) as truth_vcf:

            # Note this only works for single contig alignments!
            # This is used for target arg of the pileup call
            contig = list(self.variant_configs)[0].contig
            contig_positions = {variant_config.pos for variant_config in self.variant_configs}
            contig_pos_min = min(contig_positions)
            contig_pos_max = max(contig_positions)
            target = su.COORD_FORMAT.format(contig, contig_pos_min, contig_pos_max)

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
                pc_ref_base = var_configs[0].ref
                len_ref = len(pc_ref_base)
                amenable_qnames = []

                for pileup_read in pc.pileups:

                    # If the current read position has a InDel, skip for editing
                    # WARNING: this does not check nearby positions, and may cause FNs if certain reads are chosen
                    # In this case, try changing the seed for low frequency dels
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue

                    # If the editable read position does not match the reference sequence, do not edit
                    # This should handle cases where editing InDels at the termini of read does not raise an IndexError
                    end_pos = pileup_read.query_position + len_ref
                    if end_pos > pileup_read.alignment.query_alignment_length:

                        # If we match up to the end position of the REF (e.g. in a del at terminus), we may edit
                        # Otherwise, the slice will create a string that does not match the pc_ref_base
                        end_pos = pileup_read.alignment.query_alignment_length

                    obs_bps = pileup_read.alignment.query_sequence[pileup_read.query_position:end_pos]

                    if pileup_read.alignment.query_sequence[
                       pileup_read.query_position:pileup_read.query_position + len_ref] != obs_bps:
                        continue

                    # If the edited bases intersect any synthetic primer regions, do not edit
                    if su.MASKED_BQ in pileup_read.alignment.query_qualities[
                                       pileup_read.query_position:pileup_read.query_position + len_ref]:
                        continue

                    qname_alias = self._get_qname_alias(pileup_read.alignment.query_name)
                    amenable_qnames.append(qname_alias)

                # Finally ensure both mates overlap the POS; satmut_utils call requires both mates for variant calls
                amenable_qname_counter = collections.Counter(amenable_qnames)
                amenable_qnames = {qname for (qname, qname_count) in amenable_qname_counter.items() if qname_count == 2}
                total_amenable_qnames = len(amenable_qnames)

                if total_amenable_qnames == 0:
                    _logger.warning("No amenable qnames for editing at %s." % pc_coord)
                    continue

                # Get the edit configs for all variants at the column
                for var_config in var_configs:

                    # Update the edit_configs dict and determine the truth frequency
                    expected_af = self._iterate_over_pileup_reads(
                        pc, var_config, edit_configs, amenable_qnames, total_amenable_qnames)

                    # Write the variant and its expected frequency to the truth VCF
                    self._write_variant_to_truth_vcf(truth_vcf, var_config, expected_af)

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
        """Induces a single variant in the read object.

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

        already_written = set()

        with pysam.AlignmentFile(self.editor_preprocessor.qname_sorted_bam, "rb") as in_af, \
                pysam.AlignmentFile(self.output_bam, "wb", header=in_af.header) as out_af:

            # Iterate over the qname-sorted reads, that way we write the BAM in that order for FASTQ conversion
            for align_seg in in_af.fetch(until_eof=True):

                if align_seg.is_unmapped:
                    out_af.write(align_seg)
                    continue

                qname_alias = self._get_qname_alias(align_seg.query_name)
                candidate_key = EDIT_KEY_TUPLE(qname=qname_alias, mate=su.ReadMate(align_seg.is_read1))

                # We do this for RNA reads where we have supplementary alignments.
                # The first encountered read will be edited, and we must not write a later supplement (so no dups arise
                # that would cause the resultant FASTQs to be differing lengths).
                # Note this could potentially cause issues if the first supplementary read encountered hard clips
                # and the later supplement does not, as IndexErrors could arise when trying to edit a variant
                # in the hard-clipped read (whereas the index would have been present in the later supplement).
                # However, this logic greatly simplifies logic; thus filter hard-clipped reads prior to editing to
                # avoid such edge cases
                if candidate_key in already_written:
                    continue

                # Make sure to unmask BQs if they were masked to begin with
                align_seg.query_qualities = self._unmask_quals(align_seg.query_qualities)

                if candidate_key in edit_configs:
                    self._edit(align_seg, edit_configs[candidate_key])

                # Unset the CIGAR string prior to writing as this interferes with proper object construction
                align_seg.cigarstring = None
                out_af.write(align_seg)
                already_written.add(candidate_key)

    def _write_fastqs(self):
        r"""Writes FASTQs for re-alignment.

        :return tuple: (str, str) paths of the R1 and R2 FASTQ files
        """

        r1_fastq, r2_fastq = su.bam_to_fastq(bam=self.output_bam, out_prefix=self.out_path)
        return r1_fastq, r2_fastq

    @staticmethod
    def _write_variant_to_truth_vcf(truth_vcf_fh, vc, expected_af):
        """Writes a variant to the truth VCF containing the expected frequencies for input variants.

        :param pysam.VariantFile truth_vf: truth VCF VariantFile object
        :param analysis.read_editor.VARIANT_CONFIG_TUPLE vc: namedtuple containing variant info
        :param float expected_af: expected allele frequency of the variant
        """

        info_dict = {vu.VCF_VAR_ID: vu.VCF_VAR_ID_DELIM.join(vc.contig, vc.pos, vc.ref, vc.alt),
                     vu.VCF_AF_ID: str(expected_af)}

        new_variant_record = truth_vcf_fh.new_record(
            contig=vc.contig, start=vc.pos, alleles=(vc.ref, vc.alt), info=info_dict)

        truth_vcf_fh.write(new_variant_record)

    def workflow(self):
        r"""Runs the ReadEditor workflow.

        :return tuple: (str, str, str) paths of the edited BAM, R1 FASTQ, R2 FASTQ

        Note: the edited BAM should not be used for variant calling as its CIGAR and MD tags have not been updated. \
        I prefer re-alignment rather than CIGAR and MD tag update for simplicity. However alignment incurs a cost so \
        in the future implement BAM tag updates with samtools.
        """

        _logger.info("Getting edit configs. This could take time if aligned depth across target positions is high.")
        edit_configs = self._get_edit_configs()

        _logger.info("Editing variants.")
        self._iterate_over_reads(edit_configs)

        _logger.info("Writing and gzipping FASTQs.")
        r1_fastq, r2_fastq = self._write_fastqs()
        zipped_r1_fastq = fu.gzip_file(r1_fastq)
        zipped_r2_fastq = fu.gzip_file(r2_fastq)

        # We need to realign to re-generate CIGAR and MD tags and for proper visualization of alignments in browsers
        _logger.info("Locally realigning all reads.")
        align_workflow(f1=r1_fastq, f2=r2_fastq, ref=self.ref, outdir=self.output_dir, outbam=self.output_bam, local=True)

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
