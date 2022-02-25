#!/bin/bash

# use samtools to extract sequences for given genomic regions and write to fasta file as SNV target sequences

# help page
if [ "$#" != "4" ]; then
	echo "Usage: `basename $0` [SNV_FILE] [GTF] [REF] [OUT_DIR]"
  echo "Incorrect number of arguments"
	echo ""
	echo "- Required: SNV_FILE = Tab separated text file containing gene names, transcriptID and amino acid positions for target"
	echo ""
	echo "<GENE_SYMBPL> <--TAB--> <TRANSCRIPT_ID> <--TAB--> <AA_POS_1> <--TAB--> <GAA_POS_2>"
	echo ""
	echo "- Required: GTF = Location of GTF annotation file"
	echo "- Required: REF = Location of Reference genome"
	echo "- Required: OUT_DIR = Location for writing custom targets."
	echo ""
	exit 0
fi

# Set base directory
BASE=`pwd`

# Set variables
GTF=$2
REF=$3
# OUTPUT=${4:-${BASE}/custom_targets/Fusion}
OUT_DIR=$4

# create output directory if does not already exist
if [ -d ${OUT_DIR} ]; then
  echo -e "\nTargets will be written to ${OUT_DIR}\n"
else
  # Make result directory if it doesnt exist
	mkdir -p ${OUT_DIR}
	echo -e "\nTargets will be written to ${OUT_DIR}\n"
fi

# Run make_variant_region_file.R to create text file containing genomic regions for targets
echo -e "Running R script to extract genomic regions for each target"
Rscript ${BASE}/bin/make_variant_region_file.R $1 $GTF

## For loop to run samtools on each of the created fusionRegionFiles
for regionFile in variantRegions*
  do

  ### For loop; extracts target region information from fusion_region txt file
  while read line; do

    gene_name=$(echo $line | awk '{ print $1 }' )
    transcript_id=$(echo $line | awk '{ print $2 }' )
    strand=$(echo $line | awk '{ print $3 }' )
    region1=$(echo $line | awk '{ print $4 }' )
    region2=$(echo $line | awk '{ print $5 }' )
    target_aa_pos=$(echo $line | awk '{ print $8 }' )

    if [ ${strand} == "+" ] && [ ${region2} != "Nil" ]; then
  
    # run samtools to extract sequence for each region
    echo -e "making target ${gene_name}_pos${target_aa_pos}.fa"
    samtools faidx ${REF} ${region1} ${region2} > ${OUT_DIR}/${gene_name}_pos${target_aa_pos}.fa
  
    elif [ ${strand} == "+" ] && [ ${region2} == "Nil" ]; then
  
    # run samtools to extract reverse complement sequence for each region
    echo -e "making target ${gene_name}_pos${target_aa_pos}.fa"
    samtools faidx ${REF} ${region1} > ${OUT_DIR}/${gene_name}_pos${target_aa_pos}.fa

    elif [ ${strand} == "-" ] && [ ${region2} != "Nil" ]; then
  
    # run samtools to extract reverse complement sequence for each region
    echo -e "making target ${gene_name}_pos${target_aa_pos}.fa"
    samtools faidx ${REF} -i ${region1} ${region2} > ${OUT_DIR}/${gene_name}_pos${target_aa_pos}.fa
  
    elif [ ${strand} == "-" ] && [ ${region2} == "Nil" ]; then
  
    # run samtools to extract reverse complement sequence for each region
    echo -e "making target ${gene_name}_pos${target_aa_pos}.fa"
    samtools faidx ${REF} -i ${region1} > ${OUT_DIR}/${gene_name}_pos${target_aa_pos}.fa
  
    fi
  
  
  done < ${regionFile}

  # remove FusionRegions file
  rm ${regionFile}

done

echo -e "\nTarget Generation Complete\n"
