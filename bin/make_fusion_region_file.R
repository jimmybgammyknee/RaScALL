
####################################
#                                  #
#   Create custom fusion targets   #
#                                  #
####################################

## variation on generate_fusion_target_sequence_regions using rtracklayer to read gtf file

## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org"
       options(repos=r)
})

# Install stuff when needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
BiocManager::install(c("tidyverse","Biostrings", "rtracklayer"), ask = F, update = F)
}

# load required libraries
suppressMessages(library(tidyverse, warn.conflicts = F, quietly = T))
suppressMessages(library(rtracklayer, warn.conflicts = F, quietly = T))

# collect arguments
args <- commandArgs(trailingOnly = TRUE)

# args[1] text file containing gene and exon information for targets to be generated
# args[2] location of gtf file

# load target information file
target_file <- read_delim(file = args[1], delim = "\t", comment = "#", 
                          col_types = cols(
                            gene1_symbol = col_character(),
                            gene1_id = col_character(),
                            gene1_exons = col_character(),
                            gene2_symbol = col_character(),
                            gene2_id = col_character(),
                            gene2_exons = col_character(),
                            target_length = col_double()
                          ))

# extract information listing exon ranges to include in targets
exon_ranges <- suppressWarnings(tidyr::separate(target_file, 
                                                gene1_exons, 
                                                c("gene1_exon_start","gene1_exon_end"), sep = "[:]") %>% 
                                  separate(gene2_exons, 
                                           c("gene2_exon_start","gene2_exon_end"), sep = "[:]")) %>% 
  # convert character to numeric
  mutate(gene1_exon_start = as.numeric(gene1_exon_start),
         gene1_exon_end = as.numeric(gene1_exon_end),
         gene2_exon_start = as.numeric(gene2_exon_start),
         gene2_exon_end = as.numeric(gene2_exon_end))

# load gtf file
gtf_file <- args[2]
gtf_obj <- import(gtf_file)


#############################
#                           #
#    Gene Fusion targets    #
#                           #
#############################

## Run the following function for each line in the target_file
for(i in 1:nrow(target_file)){
  
  ## extract transcript information for given gene pair, limiting to exonic regions
  # gene1
  gene1_exon_data <- data.frame(gtf_obj[gtf_obj$gene_name == target_file$gene1_symbol[i] & gtf_obj$type == "exon" & gtf_obj$transcript_id == target_file$gene1_id[i], ])
  # gene2
  gene2_exon_data <- data.frame(gtf_obj[gtf_obj$gene_name == target_file$gene2_symbol[i] & gtf_obj$type == "exon" & gtf_obj$transcript_id == target_file$gene2_id[i], ])
  
  ## limit regions to protein_coding transcripts
  ## exclude any duplicate exonic regions from separate transcripts
  
  ### Need an Error Check to confirm there is gene data & print message if not ###

  ## Exon numbers are not always listed correctly (e.g. in UCSC hg19.refGene.gtf exons on neg strand genes in reverse order)
  ## manually add exon number depending on whether gene appears on positive or negative strand
  if(gene1_exon_data$strand[1] == "-"){
    
    # remove exon number
    gene1_exon_data <- gene1_exon_data %>% 
      dplyr::select(-exon_number) %>%
      # sort exons by desc(start)
      arrange(desc(start)) %>% 
      # add exon_number using row_number
      mutate(exon_number = row_number())
    
  } else { 
    
    # gene1
    gene1_exon_data <- gene1_exon_data %>% 
      # ensure exon_number is numeric variable
      mutate(exon_number = as.numeric(exon_number))
    
  }
  
  ## Repeat editing of exon number for gene 2
  if(gene2_exon_data$strand[1] == "-"){
    
    # remove exon number
    gene2_exon_data <- gene2_exon_data %>% 
      dplyr::select(-exon_number) %>%
      # sort exons by desc(start)
      arrange(desc(start)) %>% 
      # add exon_number using row_number
      mutate(exon_number = row_number())
    
  } else {
    
    # gene1
    gene2_exon_data <- gene2_exon_data %>% 
      # ensure exon_number is numeric variable
      mutate(exon_number = as.numeric(exon_number))
    
  }
  
  # # gene1
  # gene1_exon_data <- gene1_exon_data %>% 
  #   # filter(source == "protein_coding") %>% 
  #   distinct(start, end, .keep_all = TRUE) %>% 
  #   mutate(exon_number = as.numeric(exon_number))
  # 
  # # gene2
  # gene2_exon_data <- gene2_exon_data %>% 
  #   # filter(source == "protein_coding") %>% 
  #   distinct(start, end, .keep_all = TRUE) %>% 
  #   mutate(exon_number = as.numeric(exon_number))
  
  
  ########## calculating target regions ##########
  
  ## gene1 (+) strand OR gene2 (-) strand
  # mutate(gene1_region = paste(chromosome_name,":",(exon_chrom_end - length_gene_region),"-", exon_chrom_end, sep = ""))
  ## gene2 (+) strand OR gene1 (-) strand
  # mutate(gene2_region = paste(chromosome_name,":",exon_chrom_start,"-",(exon_chrom_start + length_gene_region), sep = ""))
  
  ################################################
  
  
  ## Gene 1
  gene1_exon_data <- gene1_exon_data %>% 
    # generate field containing target region
    mutate(gene1_region = ifelse(strand == "+", 
                                 paste(seqnames, ":", (end - (target_file$target_length[i]/2 -1)), "-", (end), sep = ""), 
                                 paste(seqnames, ":", start, "-", (start + (target_file$target_length[i]/2 -1)), sep = ""))) %>% 
    # limit target regions to requested exons
    filter(exon_number >= exon_ranges$gene1_exon_start[i] & exon_number <= exon_ranges$gene1_exon_end[i]) %>% 
    # generate field with exon_number
    mutate(exon_number = paste("exon", exon_number, sep = "")) %>% 
    # select required fields
    dplyr::select(gene1_region, gene_name, strand, exon_number) %>% 
    dplyr::rename(gene1_name = gene_name, gene1_strand = strand, gene1_exon = exon_number)
  
  ## Gene 2
  gene2_exon_data <- gene2_exon_data %>% 
    # generate field containing target region
    mutate(gene2_region = ifelse(strand == "+", 
                                 paste(seqnames, ":", start, "-", (start + (target_file$target_length[i]/2 -1)), sep = ""), 
                                 paste(seqnames, ":", (end - (target_file$target_length[i]/2 -1)) ,"-", end, sep = ""))) %>% 
    # limit target regions to requested exons
    filter(exon_number >= exon_ranges$gene2_exon_start[i] & exon_number <= exon_ranges$gene2_exon_end[i]) %>% 
    # generate field with exon_number
    mutate(exon_number = paste("exon", exon_number, sep = "")) %>% 
    # select required fields
    dplyr::select(gene2_region, gene_name, strand, exon_number) %>% 
    dplyr::rename(gene2_name = gene_name, gene2_strand = strand, gene2_exon = exon_number)
  
  
  # specify output fileName
  file = paste("FusionRegions", "_", exon_ranges$gene1_symbol[i], "_exons", 
               exon_ranges$gene1_exon_start[i], "-", exon_ranges$gene1_exon_end[i], "_",
               exon_ranges$gene2_symbol[i], "_exons",
               exon_ranges$gene2_exon_start[i], "-", exon_ranges$gene2_exon_end[i],".txt", sep = "")
  # outLocation = paste("~/km_ALL_targets/custom_targets/", sep = "")        
  # outFile  = paste(outLocation, file, sep = "")
  
  
  for(x in 1:nrow(gene1_exon_data)){
    
    for(y in 1:nrow(gene2_exon_data)){
      
      bind_cols(gene1_exon_data[x,], gene2_exon_data[y,]) %>% 
        write_tsv(file = file, 
                  append = TRUE)
      
    }
  }
  
}    

    
