
###################################################################

#     Collate and filter km output for ALL clinical samples       #

###################################################################

## a sample run of seq files was downloaded, converted to k-mer count table and run through km for fusions, variants and IKZF1 del
## this script is designed to collate and filter output from each sample in the run, to make reporting at MTB easier

## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})

# Install stuff when needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install(c("tidyverse", "Biostrings", "rtracklayer"), ask = F, update = F)
}

# Load sneakily
suppressMessages(library(tidyverse, warn.conflicts = F, quietly = T))
suppressMessages(library(Biostrings, warn.conflicts = F, quietly = T))

# collect arguments
args <- commandArgs(trailingOnly = TRUE)

# args[1] output dir for sample e.g. "output/sampleFile"

# catch error when no arguments passed
if(length(args) <1) {
  
  stop("Argument (input file dir) must be supplied")
  
}

### Import the data

# list all output.csv files in results directory
resultFiles <- list.files(path = args[1], pattern = ".txt$", 
                          full.names = TRUE, recursive = TRUE)

# perform check to ensure there are result files to analyse
if (length(resultFiles) > 0) { 

  # split results according to whether they are for fusions, variants or IKZF1
  fusionFile <- resultFiles[grepl("_Fusion.txt", resultFiles)]
  SNVFile <- resultFiles[grepl("_SNV.txt", resultFiles)]
  focDelFile <- resultFiles[grepl("_focal_deletions.txt", resultFiles)]
  DUX4File <- resultFiles[grepl("_DUX4.txt", resultFiles)]
  IGHFile <- resultFiles[grepl("_IGH_fusion.txt",resultFiles)]

} else {
  
  print("No files for analysis. Check path to output files.")
  
}


  
############# Fusions #############

# check Fusion.txt file exists
if (isTRUE(file.exists(fusionFile))) {

  # read-in km fusion output
  km_fusions <- read_delim(fusionFile, delim = "\t",
                           comment = "#",
                           col_types = cols(Database = col_character(),
                                            Query = col_character(),
                                            Type = col_character(),
                                            Variant_name = col_character(),
                                            rVAF = col_double(),
                                            Expression = col_double(),
                                            Min_coverage = col_double(),
                                            Start_offset = col_double(),
                                            Sequence = col_character(),
                                            Reference_expression = col_double(),
                                            Reference_sequence = col_character(),
                                            Info = col_character()))

  # Filter for positive fusion calls
  km_fusions <- km_fusions %>%
    # exclude results for info != vs_ref or Type is ITD
    filter(Info == "vs_ref", Type != "ITD") %>%
    # select required fields
    dplyr::select(Query,Type,Variant_name,rVAF,Expression,
                  Min_coverage,Sequence,Reference_sequence) %>%
    # filter to exclude results when fusion not called
    filter(Min_coverage > 0) %>%
    # select query with greatest support (if > 1 path detected for a single fusion)
    dplyr::arrange(Query,desc(Min_coverage)) %>%
    dplyr::distinct(Query, .keep_all=TRUE) %>% 
    dplyr::arrange(desc(Min_coverage)) %>%
    # add fileName as a field in data.frame
    tibble::add_column(File = basename(fusionFile))
  
  
  # Edit df to ensure consistent reporting of all detected mutations
  km_fusions <- km_fusions %>% 
    separate(Query, c("Gene1",NA,"Gene2",NA), sep = "_", remove = FALSE) %>% 
    unite("Alteration", c(Gene1,Gene2), sep = "-") %>% 
    # arrange output by most likely breakpoint
    mutate(Type = factor(Type, levels = c("Reference","Substitution",
                                          "Insertion","Deletion","Indel"))) %>% 
    arrange(Alteration, Type, desc(Min_coverage))
  
  ### where fusion called multiple times extract most likely breakpoint(s) ###
  
  # group by Alteration and for each group extract the 'Sequence' for Type == "Reference"
  topBreakpoint <- km_fusions %>% 
    # group by alteration for instances where multiple fusions were called
    group_by(Alteration) %>% 
    # select top row for each fusion
    filter(row_number()==1) %>% 
    ungroup()
  
  ## for each fusion detected:
  for(i in 1:nrow(topBreakpoint)){
    
    # perform pattern matching and generate vector indicating if pattern found & therefore row/output should be dropped
    dropRows <- grepl(topBreakpoint$Sequence[i], km_fusions$Sequence)
    # remove rows where the reference sequence was a substring of the other breakpoints identified
    km_fusions <- km_fusions[!dropRows,]
    
  }
  
  # combine any remaining rows with the most likely breakpoints
  km_fusions <- topBreakpoint %>% 
    bind_rows(km_fusions) %>% 
    # arrange by fusion
    arrange(Alteration, desc(Min_coverage))
  

}

############# SNVs #############

# check SNV.txt file exists
if (isTRUE(file.exists(SNVFile))) {

  # read-in km SNV output
  km_SNV <- read_delim(SNVFile, delim = "\t", comment = "#",
                       col_types = cols(Database = col_character(),
                                        Query = col_character(),
                                        Type = col_character(),
                                        Variant_name = col_character(),
                                        rVAF = col_double(),
                                        Expression = col_double(),
                                        Min_coverage = col_double(),
                                        Start_offset = col_double(),
                                        Sequence = col_character(),
                                        Reference_expression = col_double(),
                                        Reference_sequence = col_character(),
                                        Info = col_character()))

  
  # Filter for queries where SNV detected
  km_SNV <- km_SNV %>%
    # exclude results for info != vs_ref
    filter(Info == "vs_ref") %>%
    # exclude results where the target sequence without variant was detected or where SNV ratio < 0.1
    filter(Type != "Reference", rVAF >= 0.1) %>%
    # select required cols
    dplyr::select(Query,Type,Variant_name,rVAF,Expression,
                  Min_coverage,Sequence,Reference_sequence) %>%
    # add fileName as a field in data.frame
    tibble::add_column(File = basename(SNVFile))
  
  ## Build SNV_anno file for targets tested against
  SNV_anno <- km_SNV %>% 
    dplyr::select(Query,Reference_sequence) %>% 
    distinct(Query, .keep_all = TRUE) %>% 
    separate(Query, c("Gene","Target_AA"), sep = "_", remove = FALSE) %>% 
    ### determine the first_aa in string
    # extract first target aa
    mutate(Target_AA = gsub("-.*","",Target_AA)) %>% 
    # remove characters leaving just the aa number
    mutate(Target_AA = gsub("[^0-9.-]", "", Target_AA)) %>% 
    # convert Target_AA from character to numeric
    mutate(Target_AA = as.numeric(Target_AA)) %>% 
    # subtract 11 to obtain the first aa position
    mutate(first_aa_pos = (Target_AA - 11)) 
    
  # use biostrings to convert Ref_sequence for RNA to amino acid
  SNV_anno <- SNV_anno %>% 
    dplyr::select(Query, Reference_sequence) %>% 
    # use tibble::deframe to convert to a named vector
    tibble::deframe() %>% 
    # convert to a DNAStringSet
    DNAStringSet() %>% 
    # translate to Amino acid sequence
    translate(no.init.codon = TRUE) %>% 
    # convert to a data.frame so easier to read
    as.vector() %>% 
    tibble::enframe() %>% 
    dplyr::rename(Query = name, Ref_AAseq = value) %>% 
    # combine with other information (i.e. first_aa)
    left_join(SNV_anno, by = "Query")



  ### Annotate with AA change ###

  # only perform annotation step if variants were detected

  ## Function to confirm that fusions were called by MetaCaller and Exit if not
  if(nrow(km_SNV) > 0) {

    # Annotate amino acid change
    AAseqs <- km_SNV %>%
      # first generate and ID consisting of Query + Variant_name
      unite("variantID", c(Query,Variant_name), sep = "_") %>% 
      # select variantID + Sequence and convert to DNAStringSet
      dplyr::select(variantID,Sequence) %>% 
      # use tibble::deframe to convert to a named vector
      tibble::deframe() %>% 
      # convert to a DNAStringSet
      DNAStringSet() %>% 
      # translate to Amino acid sequence
      translate(no.init.codon = TRUE) %>% 
      # convert to a data.frame so easier to read
      as.vector() %>% 
      tibble::enframe() %>% 
      # separate columns
      separate(name, c("Gene","GATK_AAchange","Variant_name"), sep = "_", remove = FALSE) %>% 
      unite("Query", c(Gene, GATK_AAchange), sep = "_") %>% 
      dplyr::rename(Alt_AAseq = value) %>% 
      # join this information with SNV_anno
      left_join(SNV_anno, by = "Query") %>% 
      # select just the info needed to determine AAchange
      dplyr::select(name, Alt_AAseq,Ref_AAseq,first_aa_pos)


    ## Function to annotate each row of AAseqs

    # check that there is data to be annotated
    if(nrow(AAseqs) > 0) {

      ### Loop through each varaint and get the AAChange
      for(i in 1:nrow(AAseqs)){

        # take Ref_AAseq and create a single col data fame
        Ref_AA <- strsplit(AAseqs$Ref_AAseq[i], "") %>% data.frame()
        # add name
        names(Ref_AA) <- "Ref_AA"
        # convert from factor to character
        Ref_AA$Ref_AA <- as.character(Ref_AA$Ref_AA)
        ## Add AA.position ##
        start <- AAseqs$first_aa_pos[i]
        Ref_AA <- Ref_AA %>% mutate(AA.pos = (start + 1:nrow(.) - 1))

        # repeat for the altered AA sequence
        Alt_AA <- strsplit(AAseqs$Alt_AAseq[i], "") %>% data.frame()
        # add name
        names(Alt_AA) <- "Alt_AA"
        # convert from factor to character
        Alt_AA$Alt_AA <- as.character(Alt_AA$Alt_AA)
        ## Add AA.position ##
        Alt_AA <- Alt_AA %>% mutate(AA.pos = (start + 1:nrow(.) - 1))

        ## determine the AA.Change ##
        # join the two data frames
        AA.Change <- left_join(Ref_AA, Alt_AA, by = "AA.pos") %>% 
          # exclude all AA.pos where Ref_AA matches Alt_AA
          dplyr::filter(Ref_AA != Alt_AA) %>% 
          # combine the 3 cols (Ref_AA, AA.pos, Alt_AA) to indicate the AA.Change
          unite("AA.Change", c(Ref_AA, AA.pos, Alt_AA), sep = "") %>%
          # extract the column
          dplyr::pull(AA.Change) %>%
          # if multiple AA.changes collapse into single string
          paste(., collapse = "_")
        
        ## Must account for posibility of a synonymous mutations 
        ## which does not result in an AA.Change
        
        if(AA.Change != "" ) {
          
          ## add this information to AAseqs
          AAseqs$AA.Change[i] <- AA.Change
          
        } else {
          
          ## add SynonymousVariant to AA.Change
          AAseqs$AA.Change[i] <- "SynonymousVariant"
          
        }
        
          
      }

    }


    ## re-create Query and Variant_name columns
    AAseqs <- AAseqs %>%
      tidyr::separate(name, c("Gene","Variant","Variant_name"), sep = "_") %>%
      tidyr::unite("Query", c(Gene,Variant), sep = "_")

    # Add AA.Change to km_SNV
    km_SNV <- km_SNV %>%
      # join the two data.frames
      left_join(AAseqs, by = c("Query", "Variant_name"))
    
    # Edit df to ensure consistent reporting of all detected mutations
    km_SNV <- km_SNV %>% 
      # extract gene name from Query
      separate(Query, c("Gene",NA), sep = "_", remove = FALSE) %>% 
      # combine gene name with detected AA.Change
      unite("Alteration", c(Gene,AA.Change), sep = " ") %>% 
      # drop additional cols added when annotating SNV
      dplyr::select(-Alt_AAseq,-Ref_AAseq,-first_aa_pos)

  }
    
}

########### IKZF1/ERG intgragenic deletion #############

# check FocalDel.txt file exists
if (isTRUE(file.exists(focDelFile))) {

  # read-in km output
  km_FocalDel <- read_delim(focDelFile, delim = "\t", comment = "#",
                            col_types = cols(Database = col_character(),
                                             Query = col_character(),
                                             Type = col_character(),
                                             Variant_name = col_character(),
                                             rVAF = col_double(),
                                             Expression = col_double(),
                                             Min_coverage = col_double(),
                                             Start_offset = col_double(),
                                             Sequence = col_character(),
                                             Reference_expression = col_double(),
                                             Reference_sequence = col_character(),
                                             Info = col_character()))
  
  # Filter for positive results
  km_FocalDel <- km_FocalDel %>%
    # exclude results for info != vs_ref
    filter(Info == "vs_ref") %>%
    # filter only for exact matches to target sequence
    filter(Type == "Reference") %>%
    # Results with Min_coverage < 10 represent low confidence prediction
    # Likely spurious mis-spliced transcript; remove these results
    filter(Min_coverage >= 10) %>%
    # select required cols
    dplyr::select(Query,Type,Variant_name,rVAF,Expression,
                  Min_coverage,Sequence,Reference_sequence) %>%
    # re-order results to list break-point with greatest support at top
    dplyr::arrange(desc(rVAF),desc(Min_coverage),desc(Expression)) %>%
    # add fileName as a field in data.frame
    tibble::add_column(File = basename(focDelFile))
  
  # Edit df to ensure consistent reporting of all detected mutations
  km_FocalDel <- km_FocalDel %>% 
    # extract gene name from Query
    separate(Query, c("Gene","exon1","exon2"), sep = "_", remove = FALSE) %>% 
    mutate(exon1 = gsub("exon1"," 2-",exon1), 
           exon1 = gsub("exon2"," 3-",exon1), 
           exon1 = gsub("exon3"," 4-",exon1), 
           exon1 = gsub("exon4"," 5-",exon1), 
           exon2 = gsub("exon4","3 del",exon2), 
           exon2 = gsub("exon5","4 del",exon2), 
           exon2 = gsub("exon6","5 del",exon2), 
           exon2 = gsub("exon7","6 del",exon2), 
           exon2 = gsub("exon8","7 del",exon2), 
           exon2 = gsub("exon9","8 del",exon2), 
           exon2 = gsub("exon10","9 del",exon2), 
           exon2 = gsub("exon11","10 del",exon2),
           exon2 = gsub("exon12","11 del",exon2)) %>% 
    # combine gene name with detected exon deletion
    unite("Alteration", c(Gene,exon1,exon2), sep = "")
  
}


############# DUX4 expression #############

# read-in km output files for DUX4 target sequences
if (isTRUE(file.exists(DUX4File))) {
  
  # read-in km output
  km_DUX4 <- read_delim(DUX4File, delim = "\t", comment = "#",
                        col_types = cols(Database = col_character(),
                                         Query = col_character(),
                                         Type = col_character(),
                                         Variant_name = col_character(),
                                         rVAF = col_double(),
                                         Expression = col_double(),
                                         Min_coverage = col_double(),
                                         Start_offset = col_double(),
                                         Sequence = col_character(),
                                         Reference_expression = col_double(),
                                         Reference_sequence = col_character(),
                                         Info = col_character()))


  # Filter for positive results
  km_DUX4 <- km_DUX4 %>%
    # exclude results for info != vs_ref
    filter(Info == "vs_ref") %>%
    # filter only for exact matches to target sequence with min_cov > 10
    filter(Min_coverage > 10) %>%
    # select only required columns
    dplyr::select(Query,Type,Variant_name,rVAF,Expression,
                  Min_coverage,Sequence,Reference_sequence) %>% 
    # add col 'Alteration' 
    mutate(Alteration = "") %>% 
    # add fileName as a field in data.frame
    tibble::add_column(File = basename(DUX4File)) 



####### Expression of DUX4 is limited to patients with DUX4-rearragement #######
  ####### Reasonable expression over 4+ targets would indicate this #######

##### Summarise number of targets detected and mean coverage #####

  if(nrow(km_DUX4) > 3) {

    ## Confirm sample has DUX4 over >= 4 targets and if so report back # targets and mean_coverage ##
    DUX4_summary <- km_DUX4 %>%
      # separate Query into 3 components
      tidyr::separate(Query, c("Gene","Chr","range"), sep = "_") %>%
      mutate(range = gsub("to"," - ",range)) %>%
      ## if sample has 2 hits for same range, select hit with highest coverage
      arrange(range, desc(Min_coverage),desc(Expression)) %>%
      dplyr::distinct(range, .keep_all = TRUE) %>%
      # count number targets detected in sample
      mutate(n_targets = nrow(.)) %>%
      # calculate the mean_coverage across detected targets
      mutate(mean_coverage = mean(Min_coverage)) %>%
      # select rows wish to report
      dplyr::select(File,n_targets,mean_coverage) %>%
      dplyr::distinct(., .keep_all = TRUE)
    
    ## Edit table to fit format of other results
    DUX4_summary <- DUX4_summary %>% 
      mutate(Query = as.character(NA),
             Alteration = "DUX4-rearrangement", 
             Type = as.character(NA), 
             Variant_name = as.character(NA),
             rVAF = as.integer(NA), 
             Expression = paste(n_targets, "of 9 targets", sep = " "), 
             Min_coverage = paste("mean",mean_coverage, sep = " "), 
             Sequence = as.character(NA), 
             Reference_sequence = as.character(NA)) %>% 
      # re-order / drop extra cols
      dplyr::select(Query,Alteration,Type,Variant_name,rVAF,Expression,
                    Min_coverage,Sequence,Reference_sequence,File)
      

  }

}

############# IGH fusions #############

# read-in km output files for DUX4 target sequences
if (isTRUE(file.exists(IGHFile))) {

  # read-in km output files for IGH
  km_IGH <- read_delim(IGHFile, delim = "\t", comment = "#",
                       col_types = cols(Database = col_character(),
                                        Query = col_character(),
                                        Type = col_character(),
                                        Variant_name = col_character(),
                                        rVAF = col_double(),
                                        Expression = col_double(),
                                        Min_coverage = col_double(),
                                        Start_offset = col_double(),
                                        Sequence = col_character(),
                                        Reference_expression = col_double(),
                                        Reference_sequence = col_character(),
                                        Info = col_character()))


  # Filter for positive results
  km_IGH <- km_IGH %>%
    # exclude results for info != vs_ref
    filter(Info == "vs_ref") %>%
    # filter only for matches to target sequence
    filter(Min_coverage > 0) %>%
    # select only required columns
    dplyr::select(Query,Type,Variant_name,rVAF,Expression,
                  Min_coverage,Sequence,Reference_sequence) %>%
    # select query with greatest support (if > 1 path detected for a single fusion/break-point)
    dplyr::arrange(Query,desc(Min_coverage,Expression)) %>%
    dplyr::distinct(Query, .keep_all=TRUE) %>%
    # re-order results to list break-point with greatest support at top
    dplyr::arrange(desc(rVAF),desc(Min_coverage),desc(Expression)) %>%
    # add fileName as a field in data.frame
    tibble::add_column(File = basename(IGHFile))

  # Edit df to ensure consistent reporting of all detected mutations
  km_IGH <- km_IGH %>% 
    separate(Query, c("Gene1",NA,"Gene2",NA), sep = "_", remove = FALSE) %>% 
    unite("Alteration", c(Gene1,Gene2), sep = "-")
  
}

############# Collate results into a single output csv #############

# only add DUX4 results if DUX4_summary exists (i.e. Expression of 4+ targets)

if (exists("DUX4_summary")) {
  
  # Join together all filtered data.frames
  km_results <- rbind(get0("km_fusions"),
                      get0("DUX4_summary"),
                      get0("km_DUX4"),
                      get0("km_IGH"), 
                      get0("km_SNV"), 
                      get0("km_FocalDel"))
  
} else {
  
  # Join together filtered data.frames except DUX4
  km_results <- rbind(get0("km_fusions"),
                      get0("km_IGH"), 
                      get0("km_SNV"), 
                      get0("km_FocalDel"))
  
  
}




# Extract variant type from 'File' column and add as initial column
km_results <- km_results %>%
  # remove '.txt' from end of fileName
  mutate(File = gsub(".txt","", File)) %>% 
  mutate(File = gsub("_R1.fastq.gz","",File)) %>% 
  mutate(File = gsub("focal_deletions","focalDeletion",File)) %>% 
  mutate(File = gsub("IGH_fusion","IGH",File)) %>% 
  # extract Target type from fileName
  tidyr::separate(File, c(NA, "Target_Type"), sep = "_") %>%
  # reorder columns to place Target_Type first
  select(Target_Type, Alteration, everything())


# specify output location
outputFile <- resultFiles[1] %>% gsub("_.+", "_final_variants.csv", .)

# write filtered results to csv file
write.csv(km_results, file=outputFile, row.names = FALSE)


  
  
  

