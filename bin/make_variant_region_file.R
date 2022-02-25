
#################################
#                               #
#   Create custom SNV targets   #
#                               #
#################################

## using rtracklayer to read gtf file & extract regions of interest

# load required libraries
suppressMessages(library(tidyverse, warn.conflicts = F, quietly = T))
suppressMessages(library(rtracklayer, warn.conflicts = F, quietly = T))

# collect arguments
args <- commandArgs(trailingOnly = TRUE)


# args[1] text file containing gene and amino acid position information for targets to be generated
# args[2] location of gtf file

# load target information file
target_file <- read_delim(file = args[1], delim = "\t", comment = "#", 
                          col_types = cols(
                            gene_symbol = col_character(),
                            transcript_id = col_character(),
                            pos1 = col_double(),
                            pos2 = col_double()
                          ))

# load gtf file
gtf_file <- args[2]
gtf_obj <- import(gtf_file)




## Run the following function for each line in the target_file

for(i in 1:nrow(target_file)){
  
  # extract gene info: if target within CDS extract just the CDS; else extract CDS + exon information
  if(target_file$pos1[i] >= 12){
    
    # extract CDS regions for given gene from gtf file
    CDS_info <- data.frame(gtf_obj[gtf_obj$gene_name == target_file$gene_symbol[i] & gtf_obj$transcript_id == target_file$transcript_id[i] & gtf_obj$type == "CDS"])
    
  } else {
    
    print(paste("WARNING: Target for ", paste(target_file$gene_symbol[i], target_file$pos1[i], sep = " amino acid position "),
                " cannot be generated as it requires bases from the 5' UTR.", sep = "")) 
    print(paste("Target will need to be created manually"))
    
    next
    
  }
  
  ## Exon numbers are not always listed correctly (e.g. in UCSC hg19.refGene.gtf exons on neg strand genes in reverse order)
  ## manually add exon number depending on whether gene appears on positive or negative strand
  if(CDS_info$strand[1] == "-"){
    
    # remove exon number
    CDS_info <- CDS_info %>% 
      dplyr::select(-exon_number) %>%
      # sort exons by desc(start)
      arrange(desc(start)) %>% 
      # add exon_number using row_number
      mutate(exon_number = row_number())
    
  } else { 
    
    # gene1
    CDS_info <- CDS_info %>% 
      # ensure exon_number is numeric variable
      mutate(exon_number = as.numeric(exon_number))
    
  }
  
  # determine AA positions encoded by each exon of the CDS
  CDS_info <- CDS_info %>% 
    # need to be able to add AA.position as move across exons
    # therefore must ensure exons for a given transcript_id/transcript_name are in order
    mutate(exon_number = as.numeric(exon_number)) %>% 
    dplyr::select(seqnames,start,end,width,strand,phase,gene_name,transcript_id,exon_number) %>% 
    group_by(transcript_id) %>% 
    arrange(transcript_id) %>% 
    # number rows to indicate if first CDS exon, second CDS exon ...
    mutate(transcript_exon = row_number()) %>% 
    # determine the number of complete amino acids encoded within the exon
    mutate(remainder = ifelse(phase == 0, 0, ifelse(phase == 2, 1, 2))) %>% 
    mutate(no_AA_in_exon = (width + remainder)%/%3) %>% 
    mutate(AA_end = accumulate(no_AA_in_exon, `+`)) %>% 
    mutate(AA_start = (AA_end - no_AA_in_exon + 1)) %>% 
    # calculate number extra bases at end of exon (not complete AA)
    mutate(extra_bases = (width%%3)) %>% 
    ungroup()
  
  # determine which exons the target AA.positions fall within
  first_target_AA <- target_file$pos1[i] - 11
  final_target_AA <- target_file$pos2[i] + 11
  strand <- CDS_info %>% pull(strand) %>% unique() %>% as.character()
  
  # subset CDS for target exons
  exons_target_AA <- filter(CDS_info,
                            (AA_start <= first_target_AA & AA_end >= first_target_AA) | 
                              (AA_end >= final_target_AA & AA_start <= final_target_AA))
  
  
  # are there multiple rows for a single transcript_id? if so target spans exons.
  # split targets within exon and targets spanning exons into separate data.frames
  
  # target region within single exon
  single_exon_target_info <- exons_target_AA %>% 
    dplyr::add_count(transcript_id) %>% 
    filter(n == 1) %>% 
    dplyr::select(seqnames,start,end,width,strand,phase,
                  gene_name,transcript_id,transcript_exon,
                  remainder,no_AA_in_exon,AA_start,AA_end,extra_bases)
  
  # # remove additional transcripts that have identical start & end ranges
  # single_exon_target_info <- single_exon_target_info %>% 
  #   distinct(start, end, AA_start, AA_end, .keep_all = TRUE)
  
  # target region spans multiple exons
  multi_exon_target_info <- exons_target_AA %>% 
    dplyr::add_count(transcript_id) %>% 
    filter(n > 1) %>% 
    dplyr::select(seqnames,start,end,width,strand,phase,
                  gene_name,transcript_id,remainder,
                  no_AA_in_exon,AA_start,AA_end,extra_bases)
  
  ### extract target regions for targets occuring within single exon ###
  
  if(nrow(single_exon_target_info) > 0){
    
    # initialise empty data.frame for storing target region data
    single_exon_targets <- data.frame(gene_name=character(),
                                      transcript_id=character(), 
                                      strand=character(), 
                                      target_range1=character(),
                                      target_range2=character(), 
                                      first_aa_pos=character(),
                                      final_aa_pos=character(),
                                      target_aa_pos=character(),
                                      stringsAsFactors=FALSE) 
    
    # loop over each row of data.frame
    for(j in 1:nrow(single_exon_target_info)){
      
      # Vary target region calculation depending on strand (+ / -)
      if(strand == "+") {
        
        # Vary target region calculation based on position of first aa in target
        if(first_target_AA == single_exon_target_info$AA_start[j] & single_exon_target_info$remainder[j] == 0) {
          
          single_exon_targets[j,] <- single_exon_target_info[j,] %>% 
            mutate(start_range = start, 
                   end_range = start + (((final_target_AA - first_target_AA + 1) * 3) - 1), 
                   target_range1 = paste(seqnames, ":", start_range, "-", end_range, sep = ""), 
                   target_range2 = "NA", 
                   strand = as.character(strand)) %>% 
            # select require cols
            dplyr::select(gene_name, transcript_id, strand, target_range1, target_range2) %>% 
            # add extra info regarding target
            mutate(first_aa_pos = paste(first_target_AA),
                   final_aa_pos = paste(final_target_AA), 
                   target_aa_pos = ifelse(pos1 == pos2, paste(pos1), paste(pos1,"-",pos2,sep = "")))
          
          
          
        } else if(first_target_AA == single_exon_target_info$AA_start[j] & single_exon_target_info$remainder[j] > 0) {
          
          # the first AA of the target has some bases in previous exon. 
          # must therefore create an extra target range including these bases
          
          extra_target <- filter(CDS_info, AA_end == (first_target_AA - 1), transcript_id == single_exon_target_info$transcript_id[j]) %>% 
            mutate(end_range1 = end, 
                   start_range1 = end - extra_bases + 1, 
                   target_range1 = paste(seqnames, ":", start_range1, "-", end_range1, sep = "")) %>% 
            dplyr::select(transcript_id, target_range1)
          
          single_exon_targets[j,] <- single_exon_target_info[j,] %>% 
            mutate(start_range = start, 
                   end_range = start - remainder + (((final_target_AA - first_target_AA + 1) * 3) - 1), 
                   target_range2 = paste(seqnames, ":", start_range, "-", end_range, sep = ""), 
                   strand = as.character(strand)) %>% 
            # join this information with the additional bases
            left_join(extra_target) %>% 
            # select require cols
            dplyr::select(gene_name, transcript_id, strand,target_range1, target_range2) %>% 
            # add extra info regarding target
            mutate(first_aa_pos = paste(first_target_AA),
                   final_aa_pos = paste(final_target_AA), 
                   target_aa_pos = ifelse(target_file$pos1[i] == target_file$pos2[i], paste(target_file$pos1[i]), paste(target_file$pos1[i],"-",target_file$pos2[i],sep = "")))
          
          
        } else if(first_target_AA > single_exon_target_info$AA_start[j]) {
          
          single_exon_targets[j,] <- single_exon_target_info[j,] %>% 
            mutate(start_range = start - remainder + ((first_target_AA - AA_start) * 3), 
                   end_range = start - remainder + ((final_target_AA - AA_start + 1) * 3) - 1, 
                   target_range1 = paste(seqnames, ":", start_range, "-", end_range, sep = ""), 
                   target_range2 = "NA",
                   strand = as.character(strand)) %>% 
            # select require cols
            dplyr::select(gene_name, transcript_id, strand,target_range1, target_range2) %>% 
            # add extra info regarding target
            mutate(first_aa_pos = paste(first_target_AA),
                   final_aa_pos = paste(final_target_AA), 
                   target_aa_pos = ifelse(target_file$pos1[i] == target_file$pos2[i], paste(target_file$pos1[i]), paste(target_file$pos1[i],"-",target_file$pos2[i],sep = "")))
          
        }
        
        
      } else {
        
        # Vary target region calculation based on position of first aa in target
        if(first_target_AA == single_exon_target_info$AA_start[j] & single_exon_target_info$remainder[j] == 0) {
          
          single_exon_targets[j,] <- single_exon_target_info[j,] %>% 
            mutate(end_range = (end), 
                   start_range = (end - ((final_target_AA) * 3) + 1),
                   target_range1 = paste(seqnames, ":", start_range, "-", end_range, sep = ""), 
                   target_range2 = "NA",
                   strand = as.character(strand)) %>% 
            # select require cols
            dplyr::select(gene_name, transcript_id, strand, target_range1, target_range2) %>% 
            # add extra info regarding target
            mutate(first_aa_pos = paste(first_target_AA),
                   final_aa_pos = paste(final_target_AA), 
                   target_aa_pos = ifelse(target_file$pos1[i] == target_file$pos2[i], paste(target_file$pos1[i]), paste(target_file$pos1[i],"-",target_file$pos2[i],sep = "")))
          
          
        } else if(first_target_AA == single_exon_target_info$AA_start[j] & single_exon_target_info$remainder[j] > 0) {
          
          # the first AA of the target has some bases in previous exon. 
          # must therefore create an extra target range including these bases
          
          extra_target <- filter(CDS_info, AA_end == (first_target_AA - 1)) %>% 
            mutate(start_range1 = start, 
                   end_range1 = start + extra_bases - 1, 
                   target_range1 = paste(seqnames, ":", start_range1, "-", end_range1, sep = "")) %>% 
            dplyr::select(transcript_id, target_range1)
          
          # now extract range from main portion of exon
          single_exon_targets[j,] <- single_exon_target_info[j,] %>% 
            mutate(end_range = end,
                   start_range = (end + remainder - ((final_target_AA - AA_start + 1) * 3) + 1),
                   target_range2 = paste(seqnames, ":", start_range, "-", end_range, sep = ""), 
                   strand = as.character(strand)) %>% 
            # join this information with the additional bases
            left_join(extra_target) %>% 
            # select require cols
            dplyr::select(gene_name, transcript_id, strand,target_range1, target_range2) %>% 
            # add extra info regarding target
            mutate(first_aa_pos = paste(first_target_AA),
                   final_aa_pos = paste(final_target_AA), 
                   target_aa_pos = ifelse(target_file$pos1[i] == target_file$pos2[i], paste(target_file$pos1[i]), paste(target_file$pos1[i],"-",target_file$pos2[i],sep = "")))
          
          
        } else if(first_target_AA > single_exon_target_info$AA_start[j]) {
          
          single_exon_targets[j,] <- single_exon_target_info[j,] %>% 
            mutate(end_range = (end + remainder - ((first_target_AA - AA_start) * 3)),
                   start_range = (end + remainder - ((final_target_AA - AA_start + 1) * 3) + 1),
                   target_range1 = paste(seqnames, ":", start_range, "-", end_range, sep = ""), 
                   target_range2 = "NA",
                   strand = as.character(strand)) %>% 
            # select require cols
            dplyr::select(gene_name, transcript_id, strand,target_range1, target_range2) %>% 
            # add extra info regarding target
            mutate(first_aa_pos = paste(first_target_AA),
                   final_aa_pos = paste(final_target_AA), 
                   target_aa_pos = ifelse(target_file$pos1[i] == target_file$pos2[i], paste(target_file$pos1[i]), paste(target_file$pos1[i],"-",target_file$pos2[i],sep = "")))
          
          
        }
        
      } 
      
    }  
    
  }
  
  ### extract target regions for targets spanning multiple exons ###
  
  if(nrow(multi_exon_target_info) > 0){
    
    # function to collapse columns
    coalesce_all_columns <- function(df) {
      
      tibble(
        gene_name = df$gene_name[1],
        transcript_name = df$transcript_id[1],
        strand = df$strand[1],
        seqnames = df$seqnames[1],
        start_range1 = na.omit(df$start_range1),
        end_range1 = na.omit(df$end_range1), 
        start_range2 = na.omit(df$start_range2),
        end_range2 = na.omit(df$end_range2)
      )
    }
    
    
    # Vary target region calculation depending on strand (+ / -)
    if(strand == "+") {
      
      multi_exon_targets <- multi_exon_target_info %>% 
        group_by(transcript_id) %>% 
        # calculate start and end positions for each range
        mutate(start_range1 = ifelse(first_target_AA >= AA_start, (start - remainder + ((first_target_AA - AA_start) * 3)), NA), 
               end_range1 = ifelse(first_target_AA >= AA_start, end, NA), 
               start_range2 = ifelse(final_target_AA <= AA_end, start, NA),
               end_range2 = ifelse(final_target_AA <= AA_end, (start - remainder + ((final_target_AA - AA_start + 1) * 3) - 1), NA)) %>% 
        # collapse columns
        dplyr::select(gene_name,transcript_id,strand,seqnames,
                      start_range1,end_range1,start_range2,end_range2) %>% 
        do(coalesce_all_columns(.)) %>% 
        ungroup() %>% 
        # create final target ranges
        mutate(target_range1 = paste(seqnames, ":", start_range1, "-", end_range1, sep = ""), 
               target_range2 = paste(seqnames, ":", start_range2, "-", end_range2, sep = ""), 
               strand = as.character(strand)) %>% 
        # add general information regarding the SNV target
        mutate(first_aa_pos = paste(first_target_AA),
               final_aa_pos = paste(final_target_AA), 
               target_aa_pos = ifelse(target_file$pos1[i] == target_file$pos2[i], paste(target_file$pos1[i]), paste(target_file$pos1[i],"-",target_file$pos2[i],sep = ""))) %>% 
        # select required columns in desired order
        dplyr::select(gene_name,transcript_id,strand,
                      target_range1,target_range2,
                      first_aa_pos,final_aa_pos,target_aa_pos)
      
      
    } else {
      
      multi_exon_targets <- multi_exon_target_info %>% 
        group_by(transcript_id) %>% 
        # calculate start and end positions for each range
        mutate(start_range1 = ifelse(first_target_AA >= AA_start, start, NA), 
               end_range1 = ifelse(first_target_AA >= AA_start, (end + remainder - ((first_target_AA - AA_start) * 3)), NA), 
               start_range2 = ifelse(final_target_AA <= AA_end, (end + remainder - ((final_target_AA - AA_start + 1) * 3) + 1), NA), 
               end_range2 = ifelse(final_target_AA <= AA_end, end, NA)) %>% 
        # collapse columns
        dplyr::select(gene_name,transcript_id,strand,seqnames,
                      start_range1,end_range1,start_range2,end_range2) %>% 
        do(coalesce_all_columns(.)) %>% 
        ungroup() %>% 
        # create final target ranges
        mutate(target_range1 = paste(seqnames, ":", start_range1, "-", end_range1, sep = ""), 
               target_range2 = paste(seqnames, ":", start_range2, "-", end_range2, sep = ""), 
               strand = as.character(strand)) %>% 
        # add general information regarding the SNV target
        mutate(first_aa_pos = paste(first_target_AA),
               final_aa_pos = paste(final_target_AA), 
               target_aa_pos = ifelse(target_file$pos1[i] == target_file$pos2[i], paste(target_file$pos1[i]), paste(target_file$pos1[i],"-",target_file$pos2[i],sep = ""))) %>% 
        # select required columns in desired order
        dplyr::select(gene_name,transcript_id,strand,
                      target_range1,target_range2,
                      first_aa_pos,final_aa_pos,target_aa_pos)
      
      
      
    }
    
  }
  
  ## combine target information into single data.frame ##
  if(exists("single_exon_targets") == TRUE & exists("multi_exon_targets") == TRUE){
    
    targets <- bind_rows(single_exon_targets, multi_exon_targets)
    
  } else if(exists("single_exon_targets") == FALSE) {
    
    targets <- multi_exon_targets
    
  } else if (exists("multi_exon_targets") == FALSE) {
    
    targets <- single_exon_targets
    
  }
  
  ### convert NA col to "Nil"
  targets <- targets %>% mutate(target_range2 = gsub("NA","Nil",target_range2))
  
  ### write target region information to file ###
  
  target_aa_pos <- ifelse(target_file$pos1[i] == target_file$pos2[i], paste(target_file$pos1[i]), paste(target_file$pos1[i],"-",target_file$pos2[i],sep = ""))
  
  # specify output fileName & write target region info to file
  file = paste("variantRegions", "_", target_file$gene_symbol[i], "_pos", target_aa_pos, ".txt", sep = "")
  write_tsv(targets, file = file, col_names = FALSE)
  
  ## clean-up before next run of loop
  rm(list = setdiff(ls(), c("gtf_obj","gtf_file","target_file")))
  
  ## extract region information for targets to be created 
  # current_SNV_anno <- targets %>%
  #   mutate(Query = paste(transcript_name, "_", target_aa_pos, sep = ""),
  #          ref_genome = "unk") %>%
  #   dplyr::select(Query,ref_genome,strand,target_range1,target_range2,target_aa_pos,first_aa_pos) %>%
  #   dplyr::rename(range1 = target_range1, range2 = target_range2, target_aa = target_aa_pos)
  # 
  # SNV_anno <- rbind(get0("current_SNV_anno"),
  #                   get0("SNV_anno"))
  
}

