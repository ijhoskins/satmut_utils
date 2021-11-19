#!/usr/bin/env/python
"""Objects for generating simulated training data for error correction models."""

import logging
import numpy as np
import pandas as pd
import tempfile

from analysis.read_editor import ReadEditor
import core_utils.file_utils as fu
import core_utils.vcf_utils as vu
import prototype.variant_generator as vg

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


class ErrorCorrectionDataGenerator(object):
    """Generates data for training error correction models. Uses satmut_utils call output summary files."""

    DEFAULT_TARGETS = None
    DEFAULT_RACE_LIKE = False
    DEFAULT_PRIMERS = None
    DEFAULT_OUTDIR = "."
    DEFAULT_PREFIX = None
    DEFAULT_NVARS = None

    def __init__(self, negative_summary, mutant_summary, negative_bam, ref, gff, race_like=DEFAULT_RACE_LIKE,
                 primers=DEFAULT_PRIMERS, outdir=DEFAULT_OUTDIR, output_prefix=DEFAULT_PREFIX, nvars=DEFAULT_NVARS,
                 haplotypes=vg.VariantGenerator.DEFAULT_HAPLO, haplotype_len=vg.VariantGenerator.DEFAULT_HAPLO_LEN,
                 random_seed=vu.VcfSubsampler.DEFAULT_SEED):
        """Constructor for ErrorCorrectionDataGenerator.

        :param str negative_summary: vcf.summary.txt file for the negative control library
        :param str mutant_summary: vcf.summary.txt file for a mutant library
        :param str negative_bam: endogenous BAM file to induce into
        :param str ref: reference FASTA used in alignment/variant calling
        :param str gff: reference GTF of the transcript that was mutagenized
        :param bool race_like: is the data produced by RACE-like (e.g. AMP) data? Default False.
        :param str | None primers: primer bed file for BQ masking
        :param str outdir: Optional output directory to write generated VCFs and edited FASTQs, BAM
        :param str | None output_prefix: Optional output prefix for the FASTQ(s) and BAM; if None, use same prefix as VCF
        :param int | None nvars: number of variants requested; set to None for equivalent nbases as in the mutant_summary
        :param bool haplotypes: should haplotypes be created with uniform number to codon variants? Default False.
        :param int haplotype_len: max length to create haplotypes. Must be no longer than read length. Default 12.
        :param int random_seed: seed for variant sampling. Default 9.
        """

        self.negative_summary = negative_summary
        self.mutant_summary = mutant_summary
        self.negative_bam = negative_bam
        self.ref = ref
        self.gff = gff
        self.race_like = race_like
        self.primers = primers
        self.outdir = outdir
        self.output_prefix = output_prefix
        self.haplotypes = haplotypes
        self.haplotype_len = haplotype_len
        self.random_seed = random_seed
        self.nvars = nvars
        if nvars is None:
            self.nvars = self._get_fp_count()

        self.vg = vg.VariantGenerator(
            gff=gff, ref=ref, haplotypes=haplotypes, haplotype_len=haplotype_len, outdir=outdir, random_seed=random_seed)

    def _get_fp_count(self):
        """Gets the number of false positive mismatches in the negative control.

        :return int: number of false positive mismatches in the background
        """

        nc_df = pd.read_csv(self.negative_summary, sep=fu.FILE_DELIM)

        nc_df[vu.VCF_MM_ID] = nc_df[vu.VAR_POS].astype(str) + vu.VCF_VAR_ID_DELIM + \
                              nc_df[vu.VCF_REF_NT_ID] + vu.VCF_VAR_ID_DELIM + \
                              nc_df[vu.VCF_ALT_NT_ID]

        num_mismatches = int(nc_df[vu.VCF_MM_ID].drop_duplicates().count())
        return num_mismatches

    def _get_mutant_caf_estimates(self):
        """Gets the median and SD of log10 CAFs of each variant type in a mutant fraction.

        :return dict: {str: tuple} median and SD of log10 CAFs keyed by the variant type
        """

        caf_dict = {}
        mut_df = pd.read_csv(self.mutant_summary, sep=fu.FILE_DELIM)

        # Iterate over variants for each type and store the CAF estimates
        mut_df[vu.VCF_VARTYPE_ID] = mut_df.apply(
            lambda x: str(vu.get_variant_type(ref=x[vu.VAR_REF], alt=x[vu.VAR_ALT], split_mnps=True)), axis=1)

        # Filter out variants that don't match the mutagenesis signature
        filter_df = mut_df[mut_df["MATCHES_MUT_SIG"] == "True"]

        groupby_df = filter_df[[vu.VCF_VARTYPE_ID, vu.VCF_CAF_ID]].groupby(vu.VCF_VARTYPE_ID)

        for vartype, cafs in groupby_df:
                log10_cafs = np.log10(cafs[vu.VCF_CAF_ID])
                caf_dict[vartype] = (float(log10_cafs.median()), float(log10_cafs.std()))

        return caf_dict

    def workflow(self, trx_id, targets=vg.VariantGenerator.DEFAULT_TARGETS, out_vcf=vg.VariantGenerator.DEFAULT_OUTFILE,
                 var_type=vg.VariantGenerator.DEFAULT_VAR_TYPE, mnp_bases=vg.VariantGenerator.DEFAULT_MNP_BASES):
        """Runs the ErrorCorrectionDataGenerator workflow.

        :param str trx_id: transcript ID to generate variants for; only one version may be available in the input GFF
        :param str | None targets: optional target feature file. Only variants intersecting the target will be generated.
        :param str | None out_vcf: optional output VCF name for all codon permutation variants
        :param str var_type: one of {"snp", "mnp", "total"}
        :param int mnp_bases: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
        :return tuple: (str, str, str, str) paths of the truth VCF, edited BAM, R1 FASTQ, R2 FASTQ
        """

        _logger.info("Starting error correction data generation workflow.")

        _logger.info("Estimating truth set parameters from positive and negative control variant calls.")
        fp_nbases = self._get_fp_count()
        caf_estimates = self._get_mutant_caf_estimates()

        _logger.info("%i false positive mismatched bases counted in %s" % (fp_nbases, self.negative_summary))

        _logger.info("Estimates of the median and standard deviation of log10 concordant frequencies are %s"
                     % str(caf_estimates))

        _logger.info("Generating codon permutation variants.")
        all_codon_permuts = self.vg.workflow(
            trx_id=trx_id, targets=targets, outfile=out_vcf, var_type=var_type, mnp_bases=mnp_bases)

        _logger.info("Subsampling variants/bases for true positives.")
        vs = vu.VcfSubsampler(all_codon_permuts, outdir=self.outdir, random_seed=self.random_seed)
        subsamp_vcf = vs.subsample_bases(nbases=fp_nbases)

        _logger.info("Annotating variants using estimated frequency parameters from the mutant summary file.")
        vp = vu.VcfPreprocessor(in_vcf=subsamp_vcf, caf_estimates=caf_estimates, ref=self.ref, outdir=self.outdir)

        _logger.info("Editing variants into the negative control alignments.")
        outbam, r1_fastq, r2_fastq = ReadEditor(
            bam=self.negative_bam, variants=vp.out_vcf, ref=self.ref, race_like=self.race_like, primers=self.primers,
            output_dir=self.outdir, output_prefix=self.output_prefix).workflow()

        _logger.info("Completed error correction data generation workflow.")

        return subsamp_vcf, outbam, r1_fastq, r2_fastq
