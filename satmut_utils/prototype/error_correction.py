#!/usr/bin/env python3
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
logger = logging.getLogger(__name__)


class ErrorCorrectionDataGenerator(object):
    """Generates data for training error correction models. Uses satmut_utils call output summary files."""

    DEFAULT_TARGETS = None
    DEFAULT_RACE_LIKE = False
    DEFAULT_PRIMERS = None
    DEFAULT_OUTDIR = "."
    DEFAULT_PREFIX = None
    DEFAULT_NTHREADS = 0

    def __init__(self, negative_summary, mutant_summary, negative_bam, ref, gff, race_like=DEFAULT_RACE_LIKE,
                 primers=DEFAULT_PRIMERS, outdir=DEFAULT_OUTDIR, nthreads=DEFAULT_NTHREADS):
        """Constructor for ErrorCorrectionDataGenerator.

        :param str negative_summary: vcf.summary.txt file for the negative control library
        :param str mutant_summary: vcf.summary.txt file for a mutant library
        :param str negative_bam: endogenous BAM file to induce into
        :param str ref: reference FASTA used in alignment/variant calling
        :param str gff: reference GTF of the transcript that was mutagenized
        :param bool race_like: is the data produced by RACE-like (e.g. AMP) data? Default False.
        :param str | None primers: primer bed file for BQ masking
        :param str outdir: Optional output directory to write generated VCFs and edited FASTQs, BAM
        :param int nthreads: Number of threads to use for BAM operations. Default 0 (autodetect).
        """

        self.negative_summary = negative_summary
        self.mutant_summary = mutant_summary
        self.negative_bam = negative_bam
        self.ref = ref
        self.gff = gff
        self.race_like = race_like
        self.primers = primers
        self.outdir = outdir
        self.nthreads = nthreads

    def _get_fp_nvars(self):
        """Gets the number of false positive variants in the negative control.

        :return int: number of false positive variants
        """

        nc_df = pd.read_csv(self.negative_summary, sep=fu.FILE_DELIM)

        nc_df[vu.VCF_VAR_ID] = nc_df[vu.VAR_POS].astype(str) + vu.VCF_VAR_ID_DELIM + \
                               nc_df[vu.VAR_REF] + vu.VCF_VAR_ID_DELIM + \
                               nc_df[vu.VAR_ALT]

        num_vars = int(nc_df[vu.VCF_VAR_ID].drop_duplicates().count())
        return num_vars

    def _get_fp_nbases(self):
        """Gets the number of false positive mismatches in the negative control.

        :return int: number of false positive mismatches
        """

        nc_df = pd.read_csv(self.negative_summary, sep=fu.FILE_DELIM)
        n_mismatches = len(nc_df)
        return n_mismatches

    def _get_mutant_caf_estimates(self):
        """Gets the mean and SD of log10 CAFs of each variant type in the positive control.

        :return dict: {str: tuple} mean and SD of log10 CAFs keyed by the variant type
        """

        caf_dict = {}
        mut_df = pd.read_csv(self.mutant_summary, sep=fu.FILE_DELIM)

        # Iterate over variants for each type and store the CAF estimates
        mut_df[vu.VCF_VARTYPE_ID] = mut_df.apply(
            lambda x: str(vu.get_variant_type(ref=x[vu.VAR_REF], alt=x[vu.VAR_ALT], split_mnps=False)), axis=1)

        # Filter out variants that don't match the mutagenesis signature
        filter_df = mut_df[mut_df["MATCHES_MUT_SIG"] == "True"]

        groupby_df = filter_df[[vu.VCF_VARTYPE_ID, vu.VCF_CAF_ID]].groupby(vu.VCF_VARTYPE_ID)

        for vartype, cafs in groupby_df:
                log10_cafs = np.log10(cafs[vu.VCF_CAF_ID])
                caf_dict[vartype] = (float(log10_cafs.mean()), float(log10_cafs.std()))

        return caf_dict

    def workflow(self, trx_id, targets=vg.VariantGenerator.DEFAULT_TARGETS, mut_sig=DEFAULT_MUT_SIG,
                 out_vcf=vg.VariantGenerator.DEFAULT_OUTFILE, var_type=vg.VariantGenerator.DEFAULT_VAR_TYPE,
                 mnp_bases=vg.VariantGenerator.DEFAULT_MNP_BASES, output_prefix=DEFAULT_PREFIX,
                 haplotypes=vg.VariantGenerator.DEFAULT_HAPLO, haplotype_len=vg.VariantGenerator.DEFAULT_HAPLO_LEN,
                 random_seed=vu.VcfSubsampler.DEFAULT_SEED, buffer=ReadEditor.DEFAULT_BUFFER,
                 force_edit=ReadEditor.DEFAULT_FORCE):
        """Runs the ErrorCorrectionDataGenerator workflow.

        :param str trx_id: transcript ID to generate variants for; only one version may be available in the input GFF
        :param str | None targets: optional target feature file. Only variants intersecting the target will be generated.
        :param str mut_sig: mutagenesis signature- one of {NNN, NNK, NNS}. Default NNN.
        :param str | None out_vcf: optional output VCF name for all codon permutation variants
        :param str var_type: one of {"snp", "mnp", "total"}
        :param int mnp_bases: report for di- or tri-nt MNP? Must be either 2 or 3. Default 3.
        :param str | None output_prefix: Optional output prefix for the FASTQ(s) and BAM; if None, use same prefix as VCF
        :param bool haplotypes: should haplotypes be created with uniform number to codon variants? Default False.
        :param int haplotype_len: max length to create haplotypes. Must be no longer than read length. Default 12.
        :param int random_seed: seed for variant sampling. Default 9.
        :param int buffer: buffer about the edit span (position + REF len) to ensure lack of error before editing. Default 6.
        :param bool force_edit: flag to attempt editing of variants despite a NonconfiguredVariant exception.
        :return tuple: (str, str, str, str) paths of the truth VCF, edited BAM, R1 FASTQ, R2 FASTQ
        """

        logger.info("Starting error correction data generation workflow.")

        logger.info("Estimating truth set parameters from positive and negative control variant calls.")
        fp_nvars = self._get_fp_nvars()
        fp_nbases = self._get_fp_nbases()
        caf_estimates = self._get_mutant_caf_estimates()

        logger.info("%i false positive variants counted in %s." % (fp_nvars, self.negative_summary))
        logger.info("%i false positive mismatched bases counted in %s." % (fp_nbases, self.negative_summary))
        logger.info("Estimates of the mean and standard deviation of log10 concordant frequencies are %s"
                     % str(caf_estimates))

        logger.info("Generating all codon-permutated variants.")
        variant_generator = vg.VariantGenerator(
            gff=self.gff, ref=self.ref, mut_sig=mut_sig, haplotypes=haplotypes, haplotype_len=haplotype_len,
            outdir=self.outdir, random_seed=random_seed)

        all_codon_permuts = variant_generator.workflow(
            trx_id=trx_id, targets=targets, outfile=out_vcf, var_type=var_type, mnp_bases=mnp_bases)

        logger.info("Subsampling variants/bases for true positives.")
        vs = vu.VcfSubsampler(cf=all_codon_permuts, outdir=self.outdir, random_seed=random_seed)
        subsamp_vcf = vs.subsample_bases(nvars=fp_nvars, nbases=fp_nbases)

        logger.info("Annotating variants using estimated frequency parameters from the mutant summary file.")
        vp = vu.VcfPreprocessor(in_vcf=subsamp_vcf, caf_estimates=caf_estimates, ref=self.ref, outdir=self.outdir)

        logger.info("Editing variants into the negative control alignments.")
        outbam, r1_fastq, r2_fastq = ReadEditor(
            bam=self.negative_bam, variants=vp.out_vcf, ref=self.ref, race_like=self.race_like, primers=self.primers,
            output_dir=self.outdir, output_prefix=output_prefix, random_seed=random_seed,
            buffer=buffer, force_edit=force_edit, nthreads=self.nthreads).workflow()

        logger.info("Completed error correction data generation workflow.")

        return subsamp_vcf, outbam, r1_fastq, r2_fastq
