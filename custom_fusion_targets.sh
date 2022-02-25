#!/bin/bash

# use samtools to extract sequences for given genomic regions and write to fasta file as gene fusion target sequences

# help page
if [ "$#" != "4" ]; then
	echo "Usage: `basename $0` [FUSION_FILE] [GTF] [REF] [OUT_DIR]"
  echo "Incorrect number of arguments"
	echo ""
	echo "- Required: FUSION_FILE = Tab separated text file containing gene names, transcript_id, exons and target length"
	echo ""
	echo "<GENE1_NAME> <--TAB--> <GENE1_TRANSCRIPT_ID> <--TAB--> <GENE1_EXONS> <--TAB--> <GENE2_NAME> <--TAB--> <GENE2_TRANSCRIPT_ID> <--TAB--> <GENE2_EXONS> <--TAB--> <TARGET_LENGTH>"
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
OUT_DIR=$4

# create output directory if does not already exist
if [ -d ${OUT_DIR} ]; then
  echo -e "\nTargets will be written to ${OUT_DIR}\n"
else
  # Make result directory if it doesnt exist
	mkdir -p ${OUT_DIR}
	echo -e "\nTargets will be written to ${OUT_DIR}\n"
fi


# Run make_fusion_region_file.R to create text file containing genomic regions for targets
echo -e "Running R script to extract genomic regions for each target"
Rscript ${BASE}/bin/make_fusion_region_file.R $1 $GTF


## For loop to run samtools on each of the created fusionRegionFiles
for regionFile in FusionRegions*
  do

  ### For loop; extracts target region information from fusion_region txt file
  while read line; do

    gene1_region=$(echo $line | awk '{ print $1 }' )
    gene1_name=$(echo $line | awk '{ print $2 }' )
    gene1_strand=$(echo $line | awk '{ print $3 }' )
    gene1_exon=$(echo $line | awk '{ print $4 }' )
    gene2_region=$(echo $line | awk '{ print $5 }' )
    gene2_name=$(echo $line | awk '{ print $6 }' )
    gene2_strand=$(echo $line | awk '{ print $7 }' )
    gene2_exon=$(echo $line | awk '{ print $8 }' )
  
    if [ $gene1_strand = "+" ] && [ $gene2_strand = "+" ]; then
  
    # run samtools to extract sequence for each region
    echo -e "making target ${gene1_name}_${gene1_exon}_${gene2_name}_${gene2_exon}.fa"
    samtools faidx ${REF} ${gene1_region} ${gene2_region} > ${OUT_DIR}/${gene1_name}_${gene1_exon}_${gene2_name}_${gene2_exon}.fa
  
    elif [ $gene1_strand = "-" ] && [ $gene2_strand = "-" ]; then
  
    # run samtools to extract reverse complement sequence for each region
    echo -e "making target ${gene1_name}_${gene1_exon}_${gene2_name}_${gene2_exon}.fa"
    samtools faidx ${REF} -i ${gene1_region} ${gene2_region} > ${OUT_DIR}/${gene1_name}_${gene1_exon}_${gene2_name}_${gene2_exon}.fa
  
    elif [ $gene1_strand = "+" ] && [ $gene2_strand = "-" ]; then
  
    out1=$(samtools faidx ${REF} ${gene1_region})
    out2=$(samtools faidx ${REF} -i ${gene2_region})
  
    echo -e "making target ${gene1_name}_${gene1_exon}_${gene2_name}_${gene2_exon}.fa"
    echo -e "${out1}\n${out2}" > ${OUT_DIR}/${gene1_name}_${gene1_exon}_${gene2_name}_${gene2_exon}.fa
  
    else
  
    out1=$(samtools faidx ${REF} -i ${gene1_region})
    out2=$(samtools faidx ${REF} ${gene2_region})
  
    echo -e "making target ${gene1_name}_${gene1_exon}_${gene2_name}_${gene2_exon}.fa"
    echo -e "${out1}\n${out2}" > ${OUT_DIR}/${gene1_name}_${gene1_exon}_${gene2_name}_${gene2_exon}.fa
  
    fi


  done < ${regionFile}

  # remove FusionRegions file
  rm ${regionFile}
  
done

echo -e "\nTarget Generation Complete\n"
