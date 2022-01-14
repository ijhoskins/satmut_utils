#!/bin/bash

SATMUT_UTILS_REPO="https://github.com/ijhoskins/satmut_utils.git"

SATMUT_UTILS_REFS_REPO="https://github.com/ijhoskins/satmut_utils_refs.git"

# FTP link to Ensembl-formatted genome 
GENOME_URL="ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

REF_DIR=$(mktemp -d -t satmut_utils_refs_XXXXXXXXXX)

usage() {
cat << EOF  
Usage: ./install_satmut_utils [-htg] [-r <REF_DIR>]
Install satmut_utils and download reference files.

-h	Help
-t	Download curated human transcriptome files
-g	Download human genome FASTA
-r REF_DIR	Download files to REF_DIR

EOF
}

while getopts ":tgr:" o; do
	case "${o}" in
		t)
			GET_TRANSCRIPTOME="True"
			;;
		g)
			GET_GENOME="True"
			;;
		r)
			REF_DIR=${OPTARG}
			;;
		*)
			usage
			exit
			;;
	esac
done
shift $((OPTIND-1))

echo "Started $0"

# Need miniconda for managing environments and packages. 
if [[ ! -x $(which conda) ]]
then
	echo "conda required. See https://docs.conda.io/en/latest/miniconda.html for installation."
	exit
fi

# Create and activate the conda environment
echo "Creating satmut_utils environment"
conda env create -f satmut_utils/satmut_utils_env.yaml
conda init bash && source ~/.bash_profile
conda activate satmut_utils

if [ ! -z "$GET_TRANSCRIPTOME" ]
then
	echo "Getting curated human transcriptome files"
	mkdir -p $REF_DIR
	echo "Writing to ${!REF_DIR}"
	git clone $SATMUT_UTILS_REFS_REPO
	cp satmut_utils_refs/* $REF_DIR && gunzip $REF_DIR/*gz
fi

if [ ! -z "$GET_GENOME" ]
then
	echo "Downloading human genome FASTA"
	mkdir -p $REF_DIR
	echo "Writing to ${!REF_DIR}"
	curl -L -R -o $REF_DIR/GRCh38.fa.gz $GENOME_URL
	gunzip $REF_DIR/GRCh38.fa.gz && samtools faidx $REF_DIR/GRCh38.fa
fi

# Navigate to the satmut_utils repo and install the package
echo "Building and installing satmut_utils"
cd satmut_utils
python3 -m pip install --upgrade build && python3 -m build
python3 -m pip install -e .

echo "Activate the satmut_utils environment with \"conda activate satmut_utils\""
echo "Completed $0"

