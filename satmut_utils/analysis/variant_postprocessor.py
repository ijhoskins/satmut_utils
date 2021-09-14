#!/usr/bin/env/python
"""Objects for post-processing variants."""

import logging
import pysam
import re
import tempfile

import analysis.coordinate_mapper as cm
import core_utils.feature_file_utils as ffu
import core_utils.file_utils as fu
import core_utils.vcf_utils as vu
from definitions import *

__author__ = "Ian_Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

tempfile.tempdir = os.getenv("SCRATCH", "/tmp")
_logger = logging.getLogger(__name__)


class NnkFilterer(object):
    """Class for filtering out variants that do not match an expected POPcode mutagenesis NNK signature."""

    VALID_TRX_IDS = re.compile("|".join(["ENST", "NM"]))

    def __init__(self, vcf, gff, outfile=None, outdir="."):
        """Ctor for NnkFilterer.

        :param str vcf: VCF of variants WITH RESPECT TO THE TRANSCRIPTOME
        :param str gff: GFF/GTF to create a mapper for; must have "transcript_id" and "CDS", "stop_codon" features
        :param str | None outfile: output filename
        :param str outdir: Output directory to write results
        """

        vcf_file_ext = fu.get_extension(vcf)
        if not re.match(vu.VCF_FILETYPE, vcf_file_ext):
            raise NotImplementedError("VCF file is not of correct type.")

        gff_file_ext = fu.get_extension(gff)
        if not (re.match(ffu.GFF_FILETYPE, gff_file_ext) or re.match(ffu.GTF_FILETYPE, gff_file_ext)):
            raise NotImplementedError("GFF file is not of correct type.")

        self.vcf = vcf
        self.gff = gff
        self.outfile = outfile
        self.outdir = outdir

        if self.outfile is None:
            self.outfile = fu.replace_extension(os.path.basename(self.vcf), "NNK.filtered.vcf")

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        self.output_filepath = os.path.join(self.outdir, self.outfile)

        _logger.info("Co-opting and loading mapper object for NNK variant filtering")
        self.aam = cm.AminoAcidMapper(gff)

        _logger.info(
            "Running NNK variant filtering for %s and outputting results to %s" % (self.vcf, self.output_filepath))

        with tempfile.NamedTemporaryFile(suffix=".for.patch.vcf") as patch_vcf:
            self.filter_variants(patch_vcf)
            fu.flush_files((patch_vcf,))

            # Always have to remove the annoying INFO END tag when using pysam for creating VCFs
            vu.remove_end_info_tag(patch_vcf.name, self.output_filepath)
            vu.table_from_vcf(self.output_filepath)
            vu.tabix_index(self.output_filepath)

    def filter_variants(self, out_vcf):
        """Filters the variants that do not match the expected true-positive signature.

        :param file out_vcf: tempfile VCF for writing filtered results, prior to INFO tag patch
        """

        with pysam.VariantFile(self.vcf, "r") as in_vf:

            class_module_path = ".".join(
                [PROJECT_LAB, PROJECT_AUTHOR, os.path.basename(PROJECT_ROOT), __name__, self.__class__.__name__])

            in_vf.header.add_line(vu.VCF_METADATA_CHAR + "source=" + class_module_path)

            with pysam.VariantFile(out_vcf, "w", header=in_vf.header) as out_vf:

                for variant in in_vf.fetch():

                    # Need the transcript ID from the (potentially) multi-field contig name
                    contig_split = variant.contig.split("|")
                    if len(contig_split) == 1:
                        trx_id = contig_split[0]
                    else:
                        for contig_field in contig_split:
                            if self.VALID_TRX_IDS.match(contig_field):
                                trx_id = contig_field
                                break
                        else:
                            raise NotImplementedError("No valid transcript ID found for variant {}".format(
                                vu.VARIANT_FORMAT(variant.contig, variant.pos, variant.ref, variant.alts[0])))

                    mut_info_tuple = self.aam.get_codon_and_aa_changes(
                        trx_id=trx_id, pos=variant.pos, ref=variant.ref, alt=variant.alts[0])

                    # aa_changes will be "" if we did not have an expected NNK signature; additionally, we don't expect
                    # multi-AA changes due to variants spanning multiple codons
                    if not mut_info_tuple.matches_nnk or \
                            len(mut_info_tuple.aa_changes.split(self.aam.MUT_INFO_DELIM)) > 1:
                        continue

                    new_variant_record = out_vf.new_record(
                        contig=variant.contig, start=variant.start, alleles=(variant.ref, variant.alts[0]), info=variant.info)

                    out_vf.write(new_variant_record)
