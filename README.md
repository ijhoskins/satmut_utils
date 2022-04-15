![](./docs/satmut_utils_logo.png)

satmut\_utils is a Python package for simulation and variant calling of saturation mutagenesis data. The two main subcommands are:
1. 'sim'
2. 'call'

satmut\_utils commands are designed to simulate and call variants in paired-end, targeted sequencing reads. Alignments to a mature mRNA reference (contiguous, spliced coding sequence with possible untranslated regions) are expected. Genome-wide and transcriptome-wide variant calling is not supported.

[Installation](#Installation)

[Reference files](#Reference-files)

[Code examples](#Code-examples)

[Run sim](#Run-sim)

[Run call](#Run-call)

[Tests](#Tests)

## Installation

Currently, only Unix/Linux and MacOSX operating systems are supported. To get started, follow these steps:

1. If conda is not installed, install miniconda for managing environments. See this [link](https://docs.conda.io/en/latest/miniconda.html) for installation on your particular architecture.

2. Clone the satmut\_utils repository:
```
git clone https://github.com/CenikLab/satmut_utils.git
SATMUT_ROOT="$PWD/satmut_utils"
```

3. Execute the provided shell script to generate the satmut\_utils environment, install the package, and optionally download curated reference files, which are required if using Ensembl identifiers ([see Reference files](#Reference-files)). Finally, activate the satmut\_utils environment.
```
REF_DIR="$HOME/satmut_utils_refs"
$SATMUT_ROOT/install_satmut_utils.sh -h
$SATMUT_ROOT/install_satmut_utils.sh -t -g -r "$REF_DIR" "$SATMUT_ROOT"
conda activate satmut_utils
```

You are now ready to call the command-line executable ```satmut_utils```

satmut\_utils is the primary command, with subcommands 'sim' and 'call'.

## Code examples

See the [satmut_utils manual](https://github.com/CenikLab/satmut_utils/blob/master/docs/satmut_utils_manual.md) for more detailed usage information.

Parameter help:
```
satmut_utils -h
satmut_utils sim -h
satmut_utils call -h
```

Common arguments to both 'sim' and 'call' subcommands should be provided first, then the subcommand, and then the subcommand-specific arguments.

It is recommended that a new output directory is created for each job. Default is to output to the current directory.
```OUTPUT_DIR="/tmp/satmut_utils_test"```

### Run 'sim'

Run 'sim' on *in silico* alignments to generate SNPs, MNPs, and InDels. Structural variants and gene fusions are not currently supported.
```
TEST_DIR="$SATMUT_ROOT/src/tests/test_data"
OUTPUT_DIR="/tmp/satmut_utils_test"
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed sim -f -a $TEST_DIR/CBS_sim.bam -v $TEST_DIR/CBS_sim.vcf
```

The 'sim' workflow outputs paired FASTQs, a realigned BAM file, and a truth VCF containing simulated variants and their expected counts and frequencies.

### Run 'call'

Currently, only SNP and MNP calling is supported. InDels and long-range haplotypes are not called.

Run 'call' on the simulated data by specifying an Ensembl transcript/gene ID and the directory containing curated reference files.
```
TEST_DIR="$SATMUT_ROOT/src/tests/test_data"
OUTPUT_DIR="/tmp/satmut_utils_test"
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed call -1 $TEST_DIR/CBS_sim.R1.fq.gz -2 $TEST_DIR/CBS_sim.R2.fq.gz -v -m 1
```

Here, we call variants on a simulated dataset with no adapter sequences, so we pass -v. However, in most cases the user should provide 5' and 3' adapters for trimming.

If the Ensembl ID is not in the curated set of primary transcripts, or if you want to align to a custom reference, several reference files most be provided. [Reference files](#Reference-files).
```
satmut_utils -r $TEST_DIR/CBS.fa -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed call -1 $TEST_DIR/CBS_sim.R1.fq.gz -2 $TEST_DIR/CBS_sim.R2.fq.gz -v -m 1 -g $TEST_DIR/CBS.gff -k $REF_DIR/GRCh38.fa.gz
```

The 'call' workflow produces a VCF of candidate variant calls as well as a bedgraph file reporting fragment coverage across the transcript reference. The output VCF and its corresponding tab-delimited summary.txt file contain records for each mismatched base in a MNP. See the [satmut_utils manual](https://github.com/CenikLab/satmut_utils/blob/master/docs/satmut_utils_manual.md) or the corresponding VCF header for column/field descriptions.

## Reference files

For convenience, a curated transcriptome of primary human transcripts from [APPRIS](https://apprisws.bioinfo.cnio.es/landing_page/) is provided, allowing the user to pass an Ensembl gene or transcript ID to satmut\_utils. However, if the requested Ensembl ID is not found in this set, custom reference files must be passed, which include:

A. Transcript reference (FASTA)

B. Transcript annotations (GFF)

C. GFF reference (FASTA)

Common transcript annotations in GFF format (file B) map coordinates in the genome. For this case, file A should specify a transcript FASTA and file C should specify the genome FASTA.

In typical saturation mutagenesis datasets, an intron-less coding sequence is expressed from a vector. In this case, set file A and C to the same composite (vector + coding sequence) reference FASTA, then make a custom GFF annotation (file B) with a coding sequence exon. See the user manual for more details on creating custom reference files.

## Tests

To run unit tests, execute the following from the satmut\_utils repository:

```nose2 -q```

