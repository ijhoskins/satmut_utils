![](./docs/satmut_utils_logo.png)

satmut_utils is a Python package for simulation and variant calling of saturation mutagenesis data. The two main subcommands are:
1. sim
2. call

satmut_utils commands are designed to simulate and call variants in paired-end, targeted RNA-sequencing reads. That is, alignments for a single transcript are expected. 

[Installation](#Installation)

[Reference files](#Reference-files)

[Code examples](#Code-examples)

[Run sim](#Run-sim)

[Run call](#Run-call)

[Tests](#Tests)

## Installation

Currently, only Linux and MacOSX operating systems are supported. To get started, follow these steps:

1. If conda is not installed, install miniconda for managing environments. See this [link](https://docs.conda.io/en/latest/miniconda.html) for installation on your particular architecture.

2. Clone the satmut\_utils repository:
```
git clone https://github.com/ijhoskins/satmut_utils.git
```

3. Execute the bash script to generate the satmut\_utils environment, install the package, and optionally download curated reference files, which are required if using Ensembl identifiers ([see Reference files](#Reference-files)). Finally, activate the satmut\_utils environment.
```
REF_DIR="~/satmut_utils_refs"
satmut_utils/install_satmut_utils.sh -t -g -r $REF_DIR
conda activate satmut_utils
```

You are now ready to call the command-line executable ```satmut_utils```

satmut\_utils is the primary command, with subcommands sim and call.
 

## Code examples

See docs/satmut\_utils\_manual.md for more detailed usage information.

Parameter help:
```
satmut_utils -h
satmut_utils sim -h
satmut_utils call -h
```

It is recommended that a new output directory is created for each job. Default is to output to the current directory.
```OUTPUT_DIR="/tmp/satmut_utils_test"```

Common arguments to both sim and call subcommands should be provided first, then the subcommand, and then the subcommand-specific arguments.

### Run sim

Run sim on *in silico* alignments to generate SNPs, MNPs, and InDels. Structural variants and gene fusions are not currently supported.
```
TEST_DIR="src/tests/test_data"
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed sim -f -a $TEST_DIR/CBS_sim.bam -v $TEST_DIR/CBS_sim.vcf
```

The sim workflow outputs paired FASTQs, a realigned BAM file, and a truth VCF containing simulated variants and their expected counts and frequencies.

### Run call

Run call on the simulated data by specifying an Ensembl transcript/gene ID and the directory containing curated reference files.
```
TEST_DIR="src/tests/test_data"
satmut_utils -i ENST00000398165.7 -x $REF_DIR -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed call -1 $TEST_DIR/CBS_sim.R1.fq.gz -2 $TEST_DIR/CBS_sim.R2.fq.gz -v -m 1
```

Here, we call variants on a simulated dataset with no adapter sequences, so we pass -v. However, in most cases the user should provide 5' and 3' adapters for trimming.

If the Ensembl ID is not in the curated set of primary transcripts, or if you want to align to a custom reference, several reference files most be provided. [Reference files](#Reference-files).
```
satmut_utils -r $TEST_DIR/CBS.fa -o $OUTPUT_DIR -p $TEST_DIR/CBS_sim_primers.bed call -1 $TEST_DIR/CBS_sim.R1.fq.gz -2 $TEST_DIR/CBS_sim.R2.fq.gz -v -m 1 -g $TEST_DIR/CBS.gff -k $REF_DIR/GRCh38.fa
```

The call workflow produces a VCF of candidate variant calls as well as a bedgraph file reporting fragment coverage across the transcript reference. The output VCF and its corresponding tab-delimited summary.txt file contain records for each mismatched base in an MNP. See the corresponding VCF header for column/field descriptions.

A number of useful R functions exist in src/prototype/summarization_utils.r for parsing and summarizing the output VCF summary files.

## Reference files

For convenience, a curated transcriptome of primary human transcripts from [APPRIS](https://apprisws.bioinfo.cnio.es/landing_page/) is provided, allowing the user to pass an Ensembl gene or transcript ID to satmut_utils. However, if the requested Ensembl ID is not found in this set, custom reference files must be passed, which include:

A. Transcript reference (FASTA)

B. Transcript annotations (GFF)

C. GFF reference (FASTA)

Common transcript annotations in GFF format (file B) map coordinates in the genome. For this typical case, file A should specify a transcript FASTA and file C should specify the genome FASTA.

In typical saturation mutagenesis datasets, an intron-less coding sequence, lacking endogenous untranslated regions, is expressed from a vector. In this case, set file A and C to the same composite (vector + coding sequence) reference FASTA, then make a custom GFF annotation (file B) with a single exon and CDS feature. See the user manual for more details on creating custom reference files.

## Tests

To run unit tests, execute the following from the satmut_utils repository:

```nose2 -q```

