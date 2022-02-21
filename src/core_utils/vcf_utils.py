#!/usr/bin/env python3
"""Collection of VCF/BCF related utilities."""

import aenum
import collections
import math
import logging
import numpy
import pysam
import random
import shutil
import subprocess
import tempfile
import warnings

import analysis.seq_utils as su
import core_utils.feature_file_utils as ffu
import core_utils.file_utils as fu
from satmut_utils.definitions import *

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")
logger = logging.getLogger(__name__)

VCF_FILETYPE = "vcf"
VCF_HEADER_CHAR = "#"
VCF_METADATA_CHAR = "##"
VCF_HEADER_ID_KEY = "ID"
VCF_HEADER_INFO_KEY = "INFO"
VCF_HEADER_CONTIG_KEY = "contig"
VCF_HEADER_FILTER_KEY = "FILTER"
VCF_HEADER_DESCRIPTION = "Description"
VCF_ID_HEADER_FORMAT = VCF_METADATA_CHAR + VCF_HEADER_ID_KEY + "=<Description=\"{}\">"
VCF_INFO_HEADER_FORMAT = VCF_METADATA_CHAR + VCF_HEADER_INFO_KEY + "=<ID={},Number={},Type={},Description=\"{}\">"
VCF_CONTIG_HEADER_FORMAT = VCF_METADATA_CHAR + VCF_HEADER_CONTIG_KEY + "=<ID={}>"
VAR_CONTIG = "#CHROM"
VAR_POS = "POS"
VAR_REF = "REF"
VAR_ALT = "ALT"
VARIANT_FORMAT = "{}:{}:{}:{}"
VCF_SUMMARY_EXT = "vcf.summary.txt"
VCF_HEADER_FIELDS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")

# VCF field indices
VCF_CONTIG_INDEX = 0
VCF_POS_INDEX = 1
VCF_ID_INDEX = 2
VCF_REF_INDEX = 3
VCF_ALT_INDEX = 4
VCF_QUAL_INDEX = 5
VCF_FILTER_INDEX = 6
VCF_INFO_INDEX = 7

# Variant INFO field IDs

VCF_ND_ID = "ND"  # number PCR duplicates

# VariantCaller constants
VCF_END_ID = "END"
VCF_DP_ID = "DP"
VCF_AO_ID = "AO"  # alternate observations
VCF_R1_AO_ID = "R1_AO"
VCF_R2_AO_ID = "R2_AO"
VCF_PLUS_AO_ID = "PLUS_AO"
VCF_MINUS_AO_ID = "MINUS_AO"
VCF_R1_PLUS_AO_ID = "R1_PLUS_AO"
VCF_R2_PLUS_AO_ID = "R2_PLUS_AO"
VCF_R1_MINUS_AO_ID = "R1_MINUS_AO"
VCF_R2_MINUS_AO_ID = "R2_MINUS_AO"
VCF_CAO_ID = "CAO"  # concordant alternate observations
VCF_NORM_CAO_ID = "NORM_CAO"
VCF_AF_ID = "AF"
VCF_CAF_ID = "CAF"  # concordant allele frequency
VCF_R1_MED_POS_ID = "R1_MED_RP"  # R1 median read position
VCF_R2_MED_POS_ID = "R2_MED_RP"  # R2 median read position
VCF_R1_MOD_POS_ID = "R1_MOD_RP"  # R1 mode read position
VCF_R2_MOD_POS_ID = "R2_MOD_RP"  # R2 mode read position
VCF_R1_MED_BQ_ID = "R1_MED_BQ"  # R1 median BQ
VCF_R2_MED_BQ_ID = "R2_MED_BQ"  # R2 median BQ
VCF_R1_MOD_BQ_ID = "R1_MOD_BQ"  # R1 mode BQ
VCF_R2_MOD_BQ_ID = "R2_MOD_BQ"  # R2 mode BQ
VCF_R1_MED_NM_ID = "R1_MED_NM"  # R1 median NM
VCF_R2_MED_NM_ID = "R2_MED_NM"  # R2 median NM
VCF_R1_MOD_NM_ID = "R1_MOD_NM"  # R1 mode NM
VCF_R2_MOD_NM_ID = "R2_MOD_NM"  # R2 mode NM
VCF_VARTYPE_ID = "VARTYPE"

VCF_R1_PLUS_MED_POS_ID = "R1_PLUS_MED_RP"
VCF_R1_MINUS_MED_POS_ID = "R1_MINUS_MED_RP"
VCF_R2_PLUS_MED_POS_ID = "R2_PLUS_MED_RP"
VCF_R2_MINUS_MED_POS_ID = "R2_MINUS_MED_RP"

VCF_R1_PLUS_MED_BQ_ID = "R1_PLUS_MED_BQ"
VCF_R1_MINUS_MED_BQ_ID = "R1_MINUS_MED_BQ"
VCF_R2_PLUS_MED_BQ_ID = "R2_PLUS_MED_BQ"
VCF_R2_MINUS_MED_BQ_ID = "R2_MINUS_MED_BQ"

VCF_R1_PLUS_MED_NM_ID = "R1_PLUS_MED_NM"
VCF_R1_MINUS_MED_NM_ID = "R1_MINUS_MED_NM"
VCF_R2_PLUS_MED_NM_ID = "R2_PLUS_MED_NM"
VCF_R2_MINUS_MED_NM_ID = "R2_MINUS_MED_NM"

VCF_POS_NT_ID = "POS_NT"
VCF_REF_NT_ID = "REF_NT"
VCF_ALT_NT_ID = "ALT_NT"
VCF_UP_REF_NT_ID = "UP_REF_NT"
VCF_DOWN_REF_NT_ID = "DOWN_REF_NT"

# VariantAnnotator constants
VCF_VAR_ID = "VAR_ID"
VCF_MM_ID = "MISMATCH_ID"
VCF_VAR_ID_DELIM = ":"
VCF_ANNOT_GENE_NAME = "GENE_NAME"
VCF_ANNOT_GENE_ID = "GENE_ID"
VCF_ANNOT_TRANSCRIPT_NAME = "TRX_NAME"
VCF_ANNOT_TRANSCRIPT_ID = "TRX_ID"
VCF_ANNOT_EXON_ID = "EXON_ID"
VCF_MUT_SIG_MATCH = "MATCHES_MUT_SIG"

# AAM=AminoAcidMapper constants
VCF_AAM_LOCATION_ID = "LOCATION"
VCF_AAM_CODON_REF_ID = "REF_CODON"
VCF_AAM_CODON_ALT_ID = "ALT_CODON"
VCF_AAM_AA_REF_ID = "REF_AA"
VCF_AAM_AA_ALT_ID = "ALT_AA"
VCF_AAM_AA_CHANGE_ID = "AA_CHANGE"
VCF_AAM_AA_POS_ID = "AA_POS"

# VG=VariantGenerator constants
VCF_VG_VAR_ALTS_ID = "CODON_ALTS"
VCF_VG_AA_ALTS_ID = "AA_ALTS"
VCF_VG_REF_AA_ID = "REF_AA"
VCF_VG_ALT_AA_ID = "ALT_AA"
VCF_VG_ALT_AAS_ID = "ALT_AAS"


def get_variant_type(ref, alt, split_mnps=False):
    """Gets the variant type based on ref and alt.

    :param str ref: reference base
    :param str alt: alternate
    :param bool split_mnps: should MNPs be split into di- and tri- MNP types? Default False.
    :return vcf_utils.VariantType: one of the variant enums
    """

    len_ref = len(ref)
    len_alt = len(alt)

    if len_ref == len_alt:

        if len_ref == 1:
            return VariantType.SNP

        elif split_mnps:

            if len_ref == 2:
                return VariantType.DI_NT_MNP

            n_mms = len([e1 for e1, e2 in zip(ref, alt) if e1 != e2])

            if len_ref == 3:
                if n_mms == 2:
                    return VariantType.DI_NT_MNP
                else:
                    return VariantType.TRI_NT_MNP
            else:
                # The number of these options comes from combining two variants that may individually be
                # either SNPs, di_nt_MNPs, or tri_nt_MNPs
                if n_mms == 2:
                    return VariantType.HAPLO_TWO
                if n_mms == 3:
                    return VariantType.HAPLO_THREE
                if n_mms == 4:
                    return VariantType.HAPLO_FOUR
                if n_mms == 5:
                    return VariantType.HAPLO_FIVE
                if n_mms == 6:
                    return VariantType.HAPLO_SIX
        else:
            return VariantType.MNP

    elif len_ref > len_alt:
        return VariantType.DEL

    elif len_alt > len_ref:
        return VariantType.INS


def tabix_index(feature_file, force_overwrite=True, as_csi=False):
    """Tabix index a tabular file.

    :param str feature_file: BED, GFF, or GTF
    :param bool force_overwrite: force overwrite of an existing index file? Default True
    :param bool as_csi: create a csi index? Default create tabix (tbi) index
    :return str: path of the tabix-indexed file
    :raises NotImplementedError: if the feature file is not valid for tabix indexing
    """

    input_filepath = os.path.abspath(feature_file)
    filetype = fu.get_extension(input_filepath).lower()

    # According to the pysam.tabix_index docs, a gtf filetype is not valid (though gff is); since they are the same
    # just make sure to use the valid extension so we don't run into issues in the call
    if filetype == ffu.GTF_FILETYPE:
        filetype = ffu.GFF_FILETYPE

    index_compatible_filetypes = {VCF_FILETYPE, ffu.BED_FILETYPE, su.SAM_SUFFIX, ffu.GFF_FILETYPE, "psltbl", "pileup"}

    if filetype not in index_compatible_filetypes:
        raise NotImplementedError(
            "Invalid filetype for tabix index. Valid filetypes are %s" % ",".join(index_compatible_filetypes))

    # First sort the feature file, as this is required for proper tabix indexing
    sorted_filepath = ".".join([fu.remove_extension(feature_file), "sorted", filetype])
    with open(sorted_filepath, "w") as sorted_fh:
        sorted_ff = sorted_fh.name
        ffu.sort_feature_file(feature_file, output=sorted_ff, header=True)

    zerobased = False
    if filetype == ffu.BED_FILETYPE:
        zerobased = True

    tbi_ff = pysam.tabix_index(
        filename=sorted_ff, force=force_overwrite, preset=filetype, keep_original=True,
        zerobased=zerobased, csi=as_csi, meta_char=None)

    return tbi_ff


def remove_end_info_tag(in_vcf, out_vcf=None):
    """Removes the INFO END tag that is added by pysam.VariantFile.new_record.

    :param str in_vcf: input VCF with END=0 INFO tag-value
    :param str | None out_vcf: optional output VCF; if None, will create a tempfile
    :return str: temp VCF with INFO END tag removed

    This is a patch which removes END INFO tags that interfere in visualization of non-symbolic variants in IGV.
    See # https://github.com/pysam-developers/pysam/issues/718.
    """

    # del variant.info[vu.VCF_END_ID] does not work as pysam protects this member of a VariantRecord object
    output_vcf = out_vcf
    if out_vcf is None:
        output_vcf = tempfile.NamedTemporaryFile(mode="w", suffix=".noend.vcf", delete=False).name

    with open(in_vcf, "r") as in_fh, \
            open(output_vcf, "w") as out_fh:
        subprocess.call(["sed", "s/END=0;//g"], stdin=in_fh, stdout=out_fh)

    return output_vcf


def table_from_vcf(vcf, output_filename=None):
    """Creates a table from VCF records.

    :return str: output summary filename
    """

    outname = output_filename
    if output_filename is None:
        outname = fu.replace_extension(vcf, VCF_SUMMARY_EXT)

    with pysam.VariantFile(vcf, "r") as in_vf, \
            open(outname, "w") as out_fh:

        for i, variant in enumerate(in_vf):

            if i == 0:
                table_header = list(VCF_HEADER_FIELDS)
                table_header += variant.info.keys()
                out_fh.write(fu.FILE_DELIM.join(table_header) + fu.FILE_NEWLINE)

            var_id = variant.id if variant.id is not None else su.R_COMPAT_NA
            var_qual = variant.qual if variant.qual is not None else su.R_COMPAT_NA
            var_filter = su.R_COMPAT_NA

            output_fields = [
                variant.contig, variant.pos, var_id, variant.ref, ",".join(variant.alts), var_qual, var_filter]

            output_fields += [",".join(v) if isinstance(v, tuple) else v for v in variant.info.values()]
            output_fields = list(map(str, output_fields))
            output_res = fu.FILE_DELIM.join(output_fields) + fu.FILE_NEWLINE
            out_fh.write(output_res)

    return outname


def update_header(header, info_ids=None, contigs=None):
    """Updates a pysam VCF header object with new IDs.

    :param pysam.VariantHeader header: header object to update
    :param tuple | None info_ids: list of tuples to pass to the VCF_INFO_HEADER_FORMAT string
    :param tuple | None contigs: Contigs to include in the header
    """

    if info_ids is not None:
        for new_info_id in info_ids:
            header.add_line(VCF_INFO_HEADER_FORMAT.format(*new_info_id))

    if contigs is not None:
        # Add the contig IDs
        for contig in contigs:
            header.add_line(VCF_CONTIG_HEADER_FORMAT.format(contig))


class VariantType(aenum.MultiValueEnum):
    """Enum for variant representation."""

    SNP = "snp", "SNP", "snv", "SNV"
    MNP = "mnp", "MNP", "mnv", "MNV"
    INS = "ins", "INS", "insertion", "INSERTION"
    DEL = "del", "DEL", "deletion", "DELETION"
    COMPLEX = "complex", "COMPLEX"
    INV = "inv", "INV", "inversion", "INVERSION"
    DI_NT_MNP = "di_nt_MNP"
    TRI_NT_MNP = "tri_nt_MNP"
    MULTI_NT_MNP = "multi_nt_MNP"
    HAPLO_TWO = "2_nt_HAPLO"
    HAPLO_THREE = "3_nt_HAPLO"
    HAPLO_FOUR = "4_nt_HAPLO"
    HAPLO_FIVE = "5_nt_HAPLO"
    HAPLO_SIX = "6_nt_HAPLO"


class VcfSorter(object):
    """Variant-type-aware VCF sorter."""

    VARTYPE_LIST = (VariantType.SNP, VariantType.MNP, VariantType.INV, VariantType.INS, VariantType.DEL, VariantType.COMPLEX)
    VARTYPE_DICT = dict(zip(VARTYPE_LIST, range(0, len(VARTYPE_LIST))))

    def __init__(self, in_vcf, out_vcf=None):
        """Constructor for VcfSorter.

        :param str in_vcf: input VCF
        :param str | None out_vcf: optional output VCF. If None, will return a tempfile VCF
        """

        self.in_vcf = in_vcf
        self.out_vcf = out_vcf

        if self.out_vcf is None:
            self.out_vcf = tempfile.NamedTemporaryFile(suffix=".sorted.vcf", delete=False).name

        logger.info("Sorting %s and writing to %s" % (self.in_vcf, self.out_vcf))
        self.sort_vcf()

    def _sort_key(self, variant_record):
        """Creates a custom key for sorting.

        :param pysam.VariantRecord variant_record: variant object
        :return tuple: sort order for the variant type
        """

        var_type = get_variant_type(ref=variant_record.ref, alt=variant_record.alts[0])
        type_key = self.VARTYPE_DICT[var_type]
        sort_tuple = variant_record.contig, variant_record.pos, type_key
        return sort_tuple

    def sort_vcf(self):
        """Sorts a VCF based on a custom key."""

        with pysam.VariantFile(self.in_vcf) as in_vf:

            with pysam.VariantFile(self.out_vcf, "w", header=in_vf.header) as out_vf:

                var_list = [v for v in in_vf.fetch()]
                var_list.sort(key=self._sort_key)
                [out_vf.write(v) for v in var_list]
                in_vf.reset()


class VcfNormalizer(object):
    """Left-normalizes variants and splits or joins multiallelic records, via bcftools norm."""

    def __init__(self, ref, split_multiallelics=True, join_multiallelics=False):
        r"""Constructor for VcfNormalizer.

        :param str ref: reference FASTA corresponding to the VCF
        :param bool split_multiallelics: should multiallelic variant records be split into multiple records? Default True.
        :param bool join_multiallelics: should multiallelic variant records be joined into one? Default False.
        :raises NotImplementedError: if the ref is not a FASTA file or if both split_multiallelics and \
        join_multiallelics are True.
        """

        self.ref = ref
        self.split_multiallelics = split_multiallelics
        self.join_multiallelics = join_multiallelics

        if not fu.get_extension(self.ref) in su.FASTA_FILETYPES:
            raise NotImplementedError("ref must specify a reference FASTA file.")

        if self.split_multiallelics is True and self.join_multiallelics is True:
            raise NotImplementedError(
                "Either split_multiallelics=True or join_multiallelics=True, but not both. Ensure only one is True.")

    def run_norm(self, in_vcf, out_vcf=None):
        """Runs bcftools norm on the input VCF.

        :param str in_vcf: input VCF
        :param str | None out_vcf: optional output VCF. If None, will return a tempfile VCF
        :return str: path of the output VCF

        Warning: does not handle duplicate records!
        """

        if out_vcf is None:
            out_vcf = tempfile.NamedTemporaryFile(suffix=".norm.vcf", delete=False).name

        logger.info("Running bcftools norm on %s and writing to %s" % (in_vcf, out_vcf))

        # -c checks REF base, e will exit if a REF base does not match the coordinate
        norm_call = ["bcftools", "norm", "-c", "e", "-f", self.ref, "-o", out_vcf]

        if self.split_multiallelics:
            norm_call.extend(["-m", "-"])
        elif self.join_multiallelics:
            norm_call.extend(["-m", "+"])

        norm_call.extend([in_vcf])
        subprocess.call(norm_call)

        return out_vcf


class VcfPreprocessor(object):
    """Preprocesses as a VCF by left-normalizing variants, sorting, filtering, and annotating with INFO tags."""

    DEFAULT_EXT = "norm.annot.vcf"
    DEFAULT_AF = 0.001
    DEFAULT_OMIT_TYPES = (VariantType.INS, VariantType.DEL, VariantType.COMPLEX, VariantType.INV)
    INFO_TAG_CONFIG_TUPLE = collections.namedtuple("info_tag_config", "tag_value, tag_header_fields")
    VAR_COORDS_FORMAT = "{}:{}"

    def __init__(self, in_vcf, caf_estimates, ref, outdir=".", overwrite_existing=False, omit_types=DEFAULT_OMIT_TYPES):
        r"""Constructor for VcfPreprocessor.

        :param str in_vcf: input VCF, likely with multiple variants at a number of coordinates
        :param dict caf_estimates: dict keyed by VARTYPE and valued by two-tuple of (median log10 CAF, SD log10 CAF)
        :param str ref: reference FASTA corresponding to the VCF
        :param str outdir: optional output directory to write the split VCFs. Default current directory
        :param bool overwrite_existing: Should existing INFO tags be overwritten by values. Default False.
        :param tuple omit_types: variant types to omit from the output VCF
        """

        self.in_vcf = in_vcf
        self.caf_estimates = caf_estimates
        self.ref = ref
        self.outdir = outdir
        self.overwrite_existing = overwrite_existing
        self.omit_types = set(omit_types)

        self.info_tag_dict = {
            VCF_AF_ID: self.INFO_TAG_CONFIG_TUPLE(
                self.DEFAULT_AF, (VCF_AF_ID, 1, "Float", "Allele frequency in range (0,1)."))
        }

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        self.out_vcf = os.path.join(self.outdir, os.path.basename(fu.replace_extension(self.in_vcf, self.DEFAULT_EXT)))

        self.workflow()

    def add_info_tags(self, norm_vcf):
        """Add the tags in the info_tag_dict to each variant record (universally).

        :param str norm_vcf: sorted and normalized VCF
        """

        # Update the header with the new INFO tag IDs
        new_info_ids = [itct.tag_header_fields for itct in self.info_tag_dict.values()]
        # this class will always add the variant type
        new_info_ids.append((VCF_VARTYPE_ID, 1, "String", "Variant type"))

        # Use pysam's auto file mode detections
        with pysam.VariantFile(norm_vcf) as in_vf:

            class_module_path = ".".join(
                [PROJECT_LAB, PROJECT_AUTHOR, os.path.basename(PROJECT_ROOT), __name__, self.__class__.__name__])

            in_vf.header.add_line(VCF_METADATA_CHAR + "source=" + class_module_path)

            update_header(header=in_vf.header, info_ids=tuple(new_info_ids))

            with pysam.VariantFile(self.out_vcf, "w", header=in_vf.header) as out_vf:

                for var in in_vf.fetch():

                    var_type = get_variant_type(ref=var.ref, alt=var.alts[0], split_mnps=True)

                    if var_type.value in self.omit_types:
                        continue

                    var_info = var.info
                    if VCF_VARTYPE_ID not in var_info or (VCF_VARTYPE_ID in var_info and self.overwrite_existing):
                        var_info.update({VCF_VARTYPE_ID: var_type.value})

                    for info_tag, itct in self.info_tag_dict.items():

                        if (info_tag not in var_info or (info_tag in var_info and self.overwrite_existing)) \
                                and itct.tag_value is not None:

                            # For the AF tag, sample from a normal distribution for the value
                            tag_value = itct.tag_value
                            if info_tag == VCF_AF_ID:

                                # Need a little hack to ensure AFs are never more than 1 (100%).
                                tag_value = 2
                                while tag_value > 1:
                                    # Sometimes we may have a variant type that we do not have
                                    # estimates for in the mutant fraction. In this case set to the SNP parameters.
                                    if str(var_type) not in self.caf_estimates:
                                        caf_median, caf_sd = self.caf_estimates[str(VariantType.SNP)]
                                    else:
                                        caf_median, caf_sd = self.caf_estimates[str(var_type)]
                                    log10_tag_value = numpy.random.normal(caf_median, caf_sd)
                                    tag_value = float(math.pow(10, log10_tag_value))

                            var_info.update({info_tag: tag_value})

                    new_variant_record = out_vf.new_record(
                        contig=var.contig, start=var.pos - 1, alleles=(var.ref, var.alts[0]), info=var_info)

                    out_vf.write(new_variant_record)

                in_vf.reset()

        # Patch for removing END INFO tag additions by pysam
        patch_vcf = remove_end_info_tag(self.out_vcf)
        shutil.copy2(patch_vcf, self.out_vcf)

    def workflow(self):
        """Runs the VcfPreprocessor workflow."""

        logger.info("Starting %s workflow" % __class__.__name__)

        # First sort the VCF so that records are specially sorted based on variant type
        sorted_vcf = VcfSorter(self.in_vcf).out_vcf

        # Second do a basic check and normalization on the input VCF (requires sorted input)
        vn = VcfNormalizer(ref=self.ref)
        normed_vcf = vn.run_norm(in_vcf=sorted_vcf)

        # Third add VCF INFO tags (default read editor configurations)
        logger.info("Adding INFO tags to %s and writing to %s" % (normed_vcf, self.out_vcf))
        self.add_info_tags(normed_vcf)

        logger.info("Completed %s workflow" % __class__.__name__)


class VcfSubsampler(object):
    """Class for subsampling variants from a VCFs."""

    DEFAULT_OUTDIR = "."
    DEFAULT_OUTFILE = None
    DEFAULT_EXT = "subsamp.vcf"
    DEFAULT_SEED = 9
    DEFAULT_NBASES = None
    DEFAULT_SNP_PROP = 0.25

    def __init__(self, cf, outdir=DEFAULT_OUTDIR, random_seed=DEFAULT_SEED, snp_prop=DEFAULT_SNP_PROP):
        """Constructor for VariantSubsampler.

        :param str cf: path to a BCF or VCF (possibly gzipped) file
        :param str outdir: output directory to write subsampled VCF to
        :param int random_seed: seed for variant sampling
        :param float snp_prop: SNP proportion; MNPs will be uniformly drawn from the complement. Default 0.25.
        """

        self.cf = cf
        self.outdir = outdir
        self.random_seed = random_seed
        self.snp_prop = snp_prop
        random.seed(self.random_seed)

    def _get_nvars(self):
        """Gets the number of variants in the VCF.

        :return int: number variants
        """

        tbi_fn = tabix_index(self.cf)
        nvars_call = ("bcftools", "index", "--nrecords", tbi_fn)

        with tempfile.NamedTemporaryFile(suffix=".nrecords.res") as res_file:

            subprocess.run(nvars_call, stdout=res_file)
            fu.flush_files((res_file,))
            nrecords = res_file.readlines()[0]
            return int(nrecords)

    def _get_mut_sig_match_nvars(self):
        """Gets the number of variants in the VCF that match the mutagenesis signature.

        :return int: number variants matching the signature
        """

        with pysam.VariantFile(self.cf) as in_vcf:

            mut_sig_match_counts = len(
                [var_record for var_record in in_vcf.fetch() if var_record.info[VCF_MUT_SIG_MATCH] == "True"])

            in_vcf.reset()
            return mut_sig_match_counts

    def _mix_distr(self, snp_vars, mnp_vars, nvars):
        """Mixes SNPs and MNPs according to their proportions.

        :param [pysam.VariantRecord] snp_vars: list of SNP VariantRecords
        :param [pysam.VariantRecord] mnp_vars: list of MNP VariantRecords
        :param int nvars: number of variants requested
        :return list: list of SNPs and MNPs mixed according to the snp_prop
        """

        # Determine the number of SNPs and MNPs to generate
        n_snps = int(self.snp_prop * nvars)
        n_mnps = (int((1.0 - self.snp_prop) * nvars))
        snp_samp = random.sample(snp_vars, n_snps)
        mnp_samp = random.sample(mnp_vars, n_mnps)
        all_vars = snp_samp + mnp_samp
        random.shuffle(all_vars)
        return all_vars

    def _get_snp_mnp_variants(self, nvars):
        """Gets SNP and MNP variants at the configured proportion.

        :param int nvars: number of variants requested
        :return list: SNP and MNP variants mixed according to proportion.
        """

        with pysam.VariantFile(self.cf, mode="r") as in_vcf:

            snp_variants = []
            mnp_variants = []

            for var_record in in_vcf.fetch():

                if var_record.info[VCF_MUT_SIG_MATCH] == "False":
                    continue

                var_type = get_variant_type(var_record.ref, var_record.alts[0], split_mnps=False)

                if var_type == VariantType.SNP:
                    snp_variants.append(var_record)
                elif var_type == VariantType.MNP:
                    mnp_variants.append(var_record)

            all_variants = self._mix_distr(snp_variants, mnp_variants, nvars)
            in_vcf.reset()
            return all_variants

    def subsample_variants(self, nvars, outfile=DEFAULT_OUTFILE):
        """Subsamples a specified number of variants.

        :param int nvars: number of variants requested
        :param str | None outfile: output filename
        :return str: output filename
        """

        if outfile is None:
            out_prefix = os.path.basename(self.cf)
            out_suffix = fu.add_extension("nvars" + str(nvars), self.DEFAULT_EXT)
            out_filename = os.path.join(self.outdir, fu.replace_extension(out_prefix, out_suffix))
        else:
            out_filename = os.path.join(self.outdir, outfile)

        # Get the number of variants in the input file
        emp_nvars = self._get_nvars()

        if nvars > emp_nvars:
            warnings.warn("Number of variants requested (%i) is greater than the number of variants in the VCF: %i."
                          % (nvars, emp_nvars))

        weighted_variants = self._get_snp_mnp_variants(nvars)
        vars_to_sample = set(random.sample(range(len(weighted_variants)), k=nvars))

        with pysam.VariantFile(self.cf) as in_vcf, \
                pysam.VariantFile(out_filename, "w", header=in_vcf.header) as out_vcf:

            for i, var_record in enumerate(weighted_variants):

                if vars_to_sample is not None and i not in vars_to_sample:
                    continue

                out_vcf.write(var_record)

        return out_filename

    def subsample_bases(self, nvars, nbases, outfile=DEFAULT_OUTFILE):
        """Subsamples an approximate number of mismatched bases.

        :param int nvars: number of variants to use for SNP/MNP mixing
        :param int nbases: number of mismatched bases to return
        :param str | None outfile: output filename
        :return str: output filename

        WARNING: This method returns an approximate nbases.
        """

        if outfile is None:
            out_prefix = os.path.basename(self.cf)
            out_suffix = fu.add_extension("nbases" + str(nbases), self.DEFAULT_EXT)
            out_filename = os.path.join(self.outdir, fu.replace_extension(out_prefix, out_suffix))
        else:
            out_filename = os.path.join(self.outdir, outfile)

        # Get the number of variants in the input file
        emp_nvars = self._get_nvars()

        if nvars > emp_nvars:
            warnings.warn("Number of variants requested (%i) is greater than the number of variants in the VCF: %i."
                          % (nvars, emp_nvars))

        # Make a mixture distr of SNPs and MNPs
        weighted_variants = self._get_snp_mnp_variants(nvars)

        with pysam.VariantFile(self.cf, mode="r") as in_vcf, \
                pysam.VariantFile(out_filename, mode="w", header=in_vcf.header) as out_vcf:

            base_counter = 0
            for var_record in weighted_variants:

                if base_counter >= nbases:
                    break

                vartype = get_variant_type(ref=str(var_record.ref), alt=str(var_record.alts[0]), split_mnps=True)

                if vartype == VariantType.SNP:
                    base_counter += 1
                elif vartype == VariantType.DI_NT_MNP:
                    base_counter += 2
                elif vartype == VariantType.TRI_NT_MNP:
                    base_counter += 3
                else:
                    base_counter += len(
                        [e1 for (e1, e2) in zip(str(var_record.ref), str(var_record.alts[0])) if e1 != e2])

                out_vcf.write(var_record)

        return out_filename
