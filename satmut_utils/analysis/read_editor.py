#!/usr/bin/env/python
"""Read editing objects."""

import abc
import collections
import logging
import os
import pysam
import random
import tempfile
import pickle

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


VARIANT_CONFIG_TUPLE = collections.namedtuple("VARIANT_CONFIG_TUPLE", "type, contig, pos, ref, alt, af, ie, ir")
EDIT_KEY_TUPLE = collections.namedtuple("EDIT_KEY_TUPLE", "qname, mate")
EDIT_CONFIG_TUPLE = collections.namedtuple("EDIT_CONFIG_TUPLE", "contig, pos, ref, alt, read_pos, ie, ir")

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
    """Class for editing variants into sequencing data."""

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
    MIN_BQ = 1  # omit primer-masked bases of BQ = 0
    MATCH_REQ_LIM = 3
    VAR_TAG_DELIM = "_"

    TRUTH_VCF_SUFFIX = "truth.vcf"
    NORM_VCF_SUFFIX = "norm.vcf"
    EDIT_BAM_SUFFIX = "edit.bam"
    REALIGN_BAM_SUFFIX = "realigned.bam"
    EDIT_FILT_BAM_SUFFIX = "edit.filt.bam"
    UNEDIT_VCF_SUFFIX = "unedit.vcf"

    def __init__(self, am, variants, ref, primers=DEFAULT_PRIMERS, output_dir=DEFAULT_OUTDIR,
                 output_prefix=DEFAULT_PREFIX, filter_edited=DEFAULT_FILTER, realign=DEFAULT_REALIGN,
                 random_seed=DEFAULT_SEED, nthreads=DEFAULT_NTHREADS):
        r"""Constructor for ReadEditor.

        :param str am: alignments to edit into.
        :param str variants: VCF/BCF specifying variants to edit; use the AF tag to specify AF, e.g. AF=0.1.
        :param str ref: reference FASTA. Default APPRIS primary annotation transcriptome.
        :param str | None primers: BED, GFF, or GTF file containing primers; for masking synthetic sequences for \
        accurate AF of variants under primer regions. This feature file should contain the strand of the primer. \
        Set to None for no masking. Induce variants on select strands using the IE tag, e.g. IE="+"|"-". \
        Induce variants on select reads by using the RI tag, e.g. IR="R1"|"R2".
        :param str output_dir: Optional output directory to store generated FASTQs and BAM. Default current directory.
        :param str | None output_prefix: Optional output prefix for the FASTQ(s) and BAM
        :param bool filter_edited: filter the edited BAM for those read pairs that were edited? Default True.
        :param bool realign: should full FASTQs be realigned? Needed for visualization of reads and variant calling. Default False.
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
        self.filter_edited = filter_edited
        self.realign = realign
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

                ie_val = None
                ir_val = None
                if vu.VCF_IE_ID in var.info:
                    ie_val = su.Strand(var.info[vu.VCF_IE_ID][0])

                if vu.VCF_IR_ID in var.info:
                    ir_val = su.ReadMate(var.info[vu.VCF_IR_ID][0])

                vct = VARIANT_CONFIG_TUPLE(type=vu.get_variant_type(var.ref, var.alts[0]),
                                           contig=var.contig, pos=var.pos, ref=var.ref, alt=var.alts[0],
                                           af=af_val, ie=ie_val, ir=ir_val)

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
        :param collections.defaultdict edit_configs: dict containing editing data
        :param set amenable_qnames: amenable qnames as the coordinate
        :param int total_amenable_qnames: number amenable qnames as the coordinate
        """

        qnames_to_edit = self._get_edit_qnames(amenable_qnames, variant_config, total_amenable_qnames)
        amenable_qnames -= qnames_to_edit

        # Now we have our set of qnames to edit, so store these as edit_configs
        for pileup_read in pileup_column.pileups:

            # Though we have already checked these criteria in the initial determination of variants to edit,
            # if running with IE or IR tags, we must check again as one mate may map while the other mate may have
            # a del at the variant coordinate. The same applies to the match criterion: one mate may match the
            # reference while the other has an error which we would like to preserve.
            if not pileup_read.is_del and not pileup_read.is_refskip:

                qpos = pileup_read.query_position
                qpos_end = pileup_read.query_position + len(variant_config.ref)

                # Consider a qname amenable only if it matches the reference sequence at the pileup
                obs_bp = pileup_read.alignment.query_sequence[qpos:qpos_end]

                if obs_bp != variant_config.ref or su.MASKED_BQ in pileup_read.alignment.query_qualities[qpos:qpos_end]:
                    continue

                qname_alias = self._get_qname_alias(pileup_read.alignment.query_name)

                if qname_alias in qnames_to_edit:

                    edit_key = EDIT_KEY_TUPLE(qname=qname_alias, mate=su.ReadMate(pileup_read.alignment.is_read1))

                    edit_config = EDIT_CONFIG_TUPLE(
                        contig=variant_config.contig, pos=variant_config.pos,
                        ref=variant_config.ref, alt=variant_config.alt,
                        read_pos=pileup_read.query_position,
                        ie=variant_config.ie, ir=variant_config.ir)

                    edit_configs[edit_key].append(edit_config)

    def _get_edit_configs(self):
        """Gets a dictionary of read variant editing configurations. Editing should be deterministic.

        :return collections.defaultdict: dict of lists of reads to edit
        """

        edit_configs = collections.defaultdict(list)
        random.seed(self.random_seed)

        _logger.info("Getting read pileups.")
        with pysam.AlignmentFile(self.editor_preprocessor.edit_background, "rb") as edited_background_af:

            contig = list(self.variant_configs)[0].contig
            contig_positions = {variant_config.pos for variant_config in self.variant_configs}
            contig_pos_min = min(contig_positions)
            contig_pos_max = max(contig_positions)
            target = su.COORD_FORMAT.format(contig, contig_pos_min, contig_pos_max)

            # We must use PileupColumns in iteration only
            # See https://github.com/pysam-developers/pysam/issues/746
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

                # First calculate the total number of amenable qnames at the coordinate; this is used for multiplying by
                # per-variant AFs to determine number of qnames to edit; we do this as opposed to calculating later
                # as amenable_qnames for each variant are dynamically updated based on qnames already selected for
                # editing of prior variants.
                pc_ref_base = var_configs[0].ref

                # The expr at the end says any base involved in the variant should not overlap a primer
                amenable_qnames = [
                    self._get_qname_alias(pileup_read.alignment.query_name) for pileup_read in pc.pileups
                    if not pileup_read.is_del and not pileup_read.is_refskip and
                       pileup_read.alignment.query_sequence[
                       pileup_read.query_position:pileup_read.query_position+len(pc_ref_base)] == pc_ref_base and
                       su.MASKED_BQ not in pileup_read.alignment.query_qualities[
                                           pileup_read.query_position:pileup_read.query_position + len(pc_ref_base)]
                ]

                # Ensure that both mates overlap the position as in calling satmut_utils requires both reads support
                amenable_qname_counter = collections.Counter(amenable_qnames)
                amenable_qnames = {qname for (qname, qname_count) in amenable_qname_counter.items() if qname_count == 2}
                total_amenable_qnames = len(amenable_qnames)

                if total_amenable_qnames == 0:
                    _logger.warning("No amenable qnames for editing at %s." % pc_coord)
                    continue

                # Finally get the edit configs
                for var_config in var_configs:
                    self._iterate_over_pileup_reads(pc, var_config, edit_configs, amenable_qnames, total_amenable_qnames)

        return edit_configs

    @staticmethod
    def _unmask_quals(quals):
        """Reassigns masked primer regions to non-zero BQs so the reads do not undergo 3' quality trimming.

        :param list quals: list of BQs to unmask
        :return list: list of BQs, all non-zero
        """

        unmasked_quals = [su.DEFAULT_MAX_BQ if bq == su.MASKED_BQ else bq for bq in quals]
        return unmasked_quals

    @staticmethod
    def _sort_key(edit_config_tuple):
        """Custom sort key for variants.

        :param analysis.read_editor.EDIT_CONFIG_TUPLE edit_config_tuple: config for a variant
        :return int: sort order for elements
        """

        len_ref = len(edit_config_tuple.ref)
        len_alt = len(edit_config_tuple.alt)

        # We should always attempt to edit SNPs and MNPs before InDels
        if len_ref == len_alt:
            return 0
        elif len_ref < len_alt:
            return 1
        elif len_alt > len_ref:
            return 2

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

        # Set an alignment tag so we know which reads were editd into and what variant was editd in the read
        var_tag = self.VAR_TAG_DELIM.join([variant.contig, str(variant.pos), variant.ref, variant.alt])
        align_seg.set_tag(su.SAM_EDITED_TAG, var_tag)

    def _iterate_over_reads(self, edit_configs):
        """Iterate over the reads and edit.

        :param collections.defaultdict edit_configs: dict with list of variants to edit at a specific read and position."""

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

                # We do this for RNA reads where we have supplementary alignments in the context of genomic alignment.
                # The first encountered read will be edited, and we must not write a later supplement (so no dups arise
                # that would cause the resultant FASTQs to be differing lengths).
                # Note this could potentially cause issues if the first supplementary read encountered hard clips
                # and the later supplement does not, as IndexErrors could arise when trying to edit a variant
                # in the hard-clipped read (whereas the index would have been present in the later supplement).
                # However, doing it this way greatly simplifies logic; thus filter hard-clipped reads prior to
                # editing to avoid such edge cases. The benefit in keeping the key as just the qname and mate allows
                # us to edit variants on different exons when RNA reads were mapped to genomic space.
                if candidate_key in already_written:
                    continue

                # Make sure to unmask BQs if they were masked; because we could use a masked BAM as input and
                # provide primers=None, always apply the unmasking to simplify logic
                align_seg.query_qualities = self._unmask_quals(align_seg.query_qualities)

                if candidate_key in edit_configs:

                    vars_to_edit = sorted(edit_configs[candidate_key], key=self._sort_key)
                    for variant in vars_to_edit:

                        # Here we do strand and read checking on the fly
                        if variant.ie is None and variant.ir is None:
                            self._edit(align_seg, variant)

                        elif variant.ie is not None and variant.ir is not None and \
                                variant.ie == su.Strand(align_seg.is_reverse) and \
                                variant.ir == su.ReadMate(align_seg.is_read1):
                            self._edit(align_seg, variant)

                        elif variant.ie is not None and variant.ir is None and \
                                variant.ie == su.Strand(align_seg.is_reverse):
                            self._edit(align_seg, variant)

                        elif variant.ir is not None and variant.ie is None and \
                                variant.ir == su.ReadMate(align_seg.is_read1):
                            self._edit(align_seg, variant)

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

    def _get_edited_read_pairs(self, edit_configs, filtered_bam=None):
        """Filters the edited BAM for read pairs that were edited.

        :param collections.defaultdict edit_configs: dict with list of variants to edit at a specific read and position
        :param str | None filtered_bam: output BAM filename to write to; if None, will write to configured output path
        :return str: BAM file containing edited reads.
        """

        filt_bam = filtered_bam
        if filtered_bam is None:
            filt_bam = fu.add_extension(self.out_path, self.EDIT_FILT_BAM_SUFFIX)

        edited_qnames = {k.qname for k in edit_configs.keys()}

        # Now stream through the edited BAM and filter out the edited read pairs, which are qname-sorted.
        with tempfile.NamedTemporaryFile("wb", suffix=".temp.filtered.bam", delete=False) as temp_bam, \
                pysam.AlignmentFile(self.output_bam, "rb") as in_af, \
                pysam.AlignmentFile(temp_bam, "wb", header=in_af.header) as out_af:

            for align_seq in in_af.fetch(until_eof=True):

                qname_alias = self._get_qname_alias(align_seq.query_name)

                if qname_alias in edited_qnames:
                    out_af.write(align_seq)

            temp_bam_name = temp_bam.name

            # Finally coordinate sort and index the BAM for visualization in a browser
            su.sort_and_index(temp_bam_name, filt_bam)
            fu.safe_remove((temp_bam_name,))

    def _get_unedited_vars(self, edit_configs):
        """Determines if variants not called from the truth set were due to no editing for various reasons.

        :param collections.defaultdict edit_configs: dict with list of variants to edit at a specific read and position
        :return str | None: output VCF file containing unedited variants if any exist, else None
        """

        # Get the set of intended variants
        input_var_ids = {
            self.VAR_TAG_DELIM.join([vct.contig, str(vct.pos), vct.ref, vct.alt]) for vct in self.variant_configs}

        # Get the set of variants with edit configurations; to exist in the edit config dict, there must
        # be amenable qnames for editing
        config_var_ids = {
            self.VAR_TAG_DELIM.join([edit_config_tuple.contig, str(edit_config_tuple.pos),
                                     edit_config_tuple.ref, edit_config_tuple.alt])
            for e in edit_configs.values() for edit_config_tuple in e}

        # Take the difference between input variants and variants with edit configurations
        unedited_var_ids = input_var_ids - config_var_ids

        if len(unedited_var_ids) > 0:
            unedited_vcf = fu.add_extension(self.out_path, self.UNEDIT_VCF_SUFFIX)

            with pysam.VariantFile(self.variants, "r") as in_vcf, \
                    pysam.VariantFile(unedited_vcf, "w", header=in_vcf.header) as out_vcf:

                for var_record in in_vcf.fetch():
                    if self.VAR_TAG_DELIM.join(
                        [str(var_record.chrom), str(var_record.pos),
                         str(var_record.ref), str(var_record.alts[0])]) in unedited_var_ids:

                        out_vcf.write(var_record)

            return unedited_vcf

        return None

    def _write_truth_vcf(self):
        """Writes a truth VCF containing the expected frequencies for input variants.

        :return str: filename of the truth VCF
        """

        pass

    def workflow(self):
        r"""Runs the ReadEditor workflow.

        :return tuple: (str, str, str) paths of the edited BAM, R1 FASTQ, R2 FASTQ

        Note: the edited BAM should not be used for variant calling as its CIGAR and MD tags have not been updated.
        I prefer re-alignment rather than CIGAR and MD tag update for simplicity. However alignment incurs a cost so
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

        if self.filter_edited:
            _logger.info("Filtering edited BAM for edited read pairs.")
            self._get_edited_read_pairs(edit_configs)

        # Determine if there were any variants that were not edited and write them out
        self._get_unedited_vars(edit_configs)

        # The first output BAM is useful for post-editing investigation, but we really need a realigned BAM
        # for proper visualization of alignments in some genome browsers like IGV
        # Realign reads to re-generate CIGAR and MD tags
        if self.realign:
            _logger.info("Locally realigning all reads.")
            realigned_bam_name = fu.add_extension(self.out_path, self.REALIGN_BAM_SUFFIX)
            align_workflow(
                f1=r1_fastq, f2=r2_fastq, ref=self.ref, outdir=self.output_dir, outbam=realigned_bam_name, local=True)

            return realigned_bam_name, zipped_r1_fastq, zipped_r2_fastq

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
