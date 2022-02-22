#!/usr/bin/env python3
"""Aligner interfaces."""

import os
import datetime
import logging
import re
import subprocess
import tempfile

import analysis.seq_utils as su
import core_utils.file_utils as fu
from satmut_utils.definitions import DEFAULT_QUALITY_OFFSET, PRE_V1p8_QUALITY_OFFSET, DEFAULT_TEMPDIR

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPLv3"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"


tempfile.tempdir = DEFAULT_TEMPDIR
logger = logging.getLogger(__name__)


class BowtieConfig(object):
    """Class for configuring bowtie2 call."""

    DEFAULT_LOCAL = True
    DEFAULT_NTHREADS = 1
    INDEX_EXTENSIONS_RE = re.compile(r".[0-9].bt2")
    DEFAULT_FLAGS = ["--maxins", "1000", "--no-discordant", "--fr"]
    DEFAULT_SCORES = ["--mp", "4", "--rdg", "6,4", "--rfg", "6,4"]

    def __init__(self, ref, local=DEFAULT_LOCAL, nthreads=DEFAULT_NTHREADS, quality_encoding=DEFAULT_QUALITY_OFFSET,
                 *args, **kwargs):
        """Constructor for BowtieConfig.

        :param str ref: path of indexed reference FASTA
        :param bool local: should a local alignment be done instead of global alignment (default True)
        :param int nthreads: number of threads to use in alignments. Default 1.
        :param int quality_encoding: base quality encoding for ASCII characters. Default 33.
        :param sequence args: single flags to pass, use no - prefix
        :param kwargs: two-field configuration parameters, must use long arg format as key (--arg)
        """

        self.ref = ref
        self.local = local
        self.nthreads = nthreads
        self.quality_encoding = quality_encoding
        self.args = args
        self.kwargs = kwargs

    def test_build(self):
        """Tests that valid index files exist for the reference.

        :raises RuntimeError: if no FM index files are found for the reference FASTA
        """

        fasta_dir = os.path.dirname(self.ref)
        fasta_basename = os.path.basename(self.ref)
        dir_files = os.listdir(fasta_dir)
        matches = [self.INDEX_EXTENSIONS_RE.search(f) for f in dir_files if re.match(fasta_basename, f)]
        if not any(matches):
            raise RuntimeError("No FM-index files found for the reference FASTA.")

    def build_fm_index(self):
        """Builds and FM index with bowtie2.

        :param str ref: reference FASTA
        """
        build_call = ("bowtie2-build", "--quiet", self.ref, self.ref)
        subprocess.call(build_call)


class Bowtie2(object):
    """Class for bowtie2 alignment."""

    DEFAULT_OUTDIR = "."
    DEFAULT_OUTBAM = None

    def __init__(self, config, f1, f2=None, output_dir=DEFAULT_OUTDIR, output_bam=DEFAULT_OUTBAM):
        """Constructor for Bowtie2.

        :param aligners.BowtieConfig config: config object
        :param str f1: path to FASTA or FASTQ 1
        :param str | None f2: optional path to FASTA or FASTQ 2
        :param str output_dir: optional output directory to write the output BAM to. Default current working directory.
        :param str | None output_bam: optional output BAM filename. Default None, use f1 basename.
        """

        self.config = config
        self.f1 = f1
        self.f2 = f2
        self.output_dir = output_dir

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        self.output_bam = output_bam
        if output_bam is None:
            if f2 is None:
                out_name = fu.replace_extension(self.f1, su.BAM_SUFFIX)
            else:
                out_name = fu.add_extension(os.path.basename(os.path.commonprefix([f1, f2])), su.BAM_SUFFIX)

            self.output_bam = os.path.join(output_dir, out_name)

        if f2 is None:
            rg_id = os.path.basename(fu.remove_extension(f1))
        else:
            rg_id = os.path.basename(os.path.commonprefix([fu.remove_extension(f1), fu.remove_extension(f2)]))

        self.alignment_kwargs = {"rg-id": rg_id}
        self.alignment_kwargs.update(config.kwargs)
        self.workflow()

    def _align(self):
        """Aligns reads.

        :return tuple: return codes of each process in the pipeline
        """

        with open(self.output_bam, "wb") as out_file, \
                open(os.path.join(os.path.dirname(os.path.abspath(self.output_bam)),
                                  "Bowtie2.stderr.log"), "a") as bowtie2_stderr:

            call = ["bowtie2", "-p", str(self.config.nthreads)]
            call.extend(self.config.DEFAULT_FLAGS + self.config.DEFAULT_SCORES)

            if self.config.quality_encoding == PRE_V1p8_QUALITY_OFFSET:
                call.append("--phred64")

            if self.config.local:
                call.append("--local")
            else:
                call.append("--end-to-end")

            # Add the configuration parameters
            for f in self.config.args:
                call.extend(["-" + str(f)])

            for k, i in self.alignment_kwargs.items():
                call.extend(["--{} {}".format(k, i)])

            call.extend(["-x", self.config.ref])
            if self.f2 is not None:
                call.extend(["-1", self.f1, "-2", self.f2])
            else:
                call.extend(["-U", self.f1])

            # Record a time stamp for each new alignment
            bowtie2_stderr.write(
                " ".join(["{}.{}".format(__name__, self.__class__.__name__),
                          datetime.datetime.now().strftime("%d%b%Y %I%M%p")]) + fu.FILE_NEWLINE)

            bowtie2_stderr.write(" ".join(call) + fu.FILE_NEWLINE)

            align_p = subprocess.Popen(call, stdout=subprocess.PIPE, stderr=bowtie2_stderr)

            tobam_p = subprocess.Popen(("samtools", "view", "-u", "-"),
                                       stdin=align_p.stdout, stdout=subprocess.PIPE, stderr=bowtie2_stderr)

            sort_p = subprocess.Popen(("samtools", "sort", "-"),
                                      stdin=tobam_p.stdout, stdout=out_file, stderr=bowtie2_stderr)
            sort_p.wait()
            align_p.stdout.close()
            tobam_p.stdout.close()

            return align_p.poll(), tobam_p.poll(), sort_p.poll()

    def workflow(self):
        """Runs the bowtie2 alignment workflow.

        :return str: name of the output BAM.
        """

        logger.info("Started bowtie2 aligner workflow.")

        logger.info("Writing output BAM %s" % self.output_bam)
        _ = self._align()

        su.index_bam(self.output_bam)

        logger.info("Completed bowtie2 aligner workflow.")