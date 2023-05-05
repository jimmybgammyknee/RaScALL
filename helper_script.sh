#!/bin/bash

### Run km on multiple samples against all target sets except IGH ###
  ### Also runs Rscript to collate & filter output ###

# help page
if (( $# < 1 )); then
	echo "Usage: `basename $0` [SAMPLE_SHEET] [THREADS]"
  	echo "Incorrect number of arguments"
	echo "- Required: SAMPLE_SHEET = Tab-delimited Text file containing FASTQ R1 and R2 locations"
	echo ""
	echo "<READ1 PATH> <--TAB--> <READ2 PATH>"
	echo ""
	echo "- Optional: THREADS = Number of threads for running jellyfish. Default is 12"
	exit 0
fi

# Set base directory
BASE=`pwd`
THREAD=${2:-12}

### For loop; extracts sample data from txt file
while read line; do

  R1=$(echo $line | awk '{ print $1 }')
  R2=$(echo $line | awk '{ print $2 }')
  
  # Run_km
  bash run_km.sh ${R1} ${R2} ${THREAD} ALL_targets/Fusion
  bash run_km.sh ${R1} ${R2} ${THREAD} ALL_targets/DUX4
  bash run_km.sh ${R1} ${R2} ${THREAD} ALL_targets/SNV
  bash run_km.sh ${R1} ${R2} ${THREAD} ALL_targets/focal_deletions
  bash run_km.sh ${R1} ${R2} ${THREAD} ALL_targets/IGH_fusion
  
  # Directory setup for running Rscript
  SAMPLE=`basename ${R1} | sed 's/_[^_]*$//'`
  
  # Run filter_km_output.R
  Rscript bin/filter_km_output.R output/${SAMPLE}
  
done < $1
