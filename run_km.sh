#!/bin/bash

# help page
if [ "$#" != "4" ]; then
	echo "Usage: `basename $0` [FASTQ_R1] [FASTQ_R2] [THREADS] [TARGET_DIR]"
  	echo "Incorrect number of arguments"
	echo "- Required: FASTQ_R1 = FASTQ PAIRED-END READ 1"
	echo "- Required: FASTQ_R2 = FASTQ PAIRED-END READ 2"
	echo "- Required: THREADS = Number of threads for running jellyfish"
	echo "- Required: TARGET_DIR = Location of directory containing targets to test against"
	exit 0
fi

# Setup all directories needed to run code and for results
BASE=`pwd`
VIRTUAL_ENV=${BASE}/.virtualenvs/km

TARGET_BASE=$4
TARGET_TYPE=`echo $(basename $4)`

# Directory setup
OUTPUTS=${BASE}/output
SAMPLE=`basename $1 | sed 's/_[^_]*$//'`

# km and jellyfish installed in virtual environment.
# must source the environment
source ${VIRTUAL_ENV}/bin/activate

# If statment that identifies the presence of `basename $1 .fastq.gz`.jf in the output dir 
if [ -f ${OUTPUTS}/${SAMPLE}/${SAMPLE}_countTable31.jf ]; then
	echo "Jellyfish Sample file exists"
else 
	# Make result directory if it doesnt exist
	mkdir -p ${OUTPUTS}/${SAMPLE}
	
	# Run jellyfish
	echo -e "generating count table for ${SAMPLE}"
	jellyfish count -m 31 -o ${OUTPUTS}/${SAMPLE}/${SAMPLE}_countTable31.jf -s 100M -t $3 -C -L 2 <(zcat $1) <(zcat $2)
fi

# Run km find mutation
wc -l ${OUTPUTS}/${SAMPLE}/${SAMPLE}_countTable31.jf
km find_mutation ${TARGET_BASE}/ ${OUTPUTS}/${SAMPLE}/${SAMPLE}_countTable31.jf > ${OUTPUTS}/${SAMPLE}/${SAMPLE}_${TARGET_TYPE}.txt
