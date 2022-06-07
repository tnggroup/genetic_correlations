## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      comment=NA,
                      prompt=FALSE,
                      cache=FALSE)


## ----Clear global environment-----------------------------------------------------------------------------------
remove(list = ls())


## ----Packages---------------------------------------------------------------------------------------------------
#devtools::install_github("johanzvrskovec/shru")
# library(shru)
# library(googlesheets4)

library(optparse)
library(data.table)
library(tidyverse)


## ----command line setup-----------------------------------------------------------------------------------------
column_parser <- OptionParser()

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-t", "--task"),
  type = "character",
  default ="clean",
  help = "The explicit task to run:\nclean: Clean GWAS sumstat [default %default]"
  )

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-f", "--file"),
  type = "character",
  help = "The file to process"
  )

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-c", "--code"),
  type = "character",
  help = "The code to use for the dataset from the specified file."
  )

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-l", "--label"),
  type = "character",
  help = "Label describing GWAS phenotype/sample"
  )

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-n", "--sample-size"),
  type = "numeric",
  help = "Number of individuals in the sample."
  )

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-p", "--population"),
  type = "character",
  help = "GWAS sample super population."
  )

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-o", "--output"),
  type = "character",
  help = "Path to gwas_sumstats group folder."
  )

column_options <- parse_args(column_parser)

# For testing
#/scratch/groups/gwas_sumstats/original/alcdep_afr_all.jul2017.min_n.gz
#/scratch/groups/gwas_sumstats/original/daner_PGC_BIP32b_mds7a.gz

column_options$file<-'/scratch/groups/gwas_sumstats/original/daner_PGC_BIP32b_mds7a.gz'
column_options$label<-'Test run'
column_options$code<-'TEST_XXXX'
column_options$`sample-size`<-'100000'
column_options$population<-'EUR'
column_options$output<-'/scratch/groups/gwas_sumstats'



## ----Recent date------------------------------------------------------------------------------------------------
date = Sys.Date()
date


## ----Paths------------------------------------------------------------------------------------------------------
filePaths <- column_options$file
filePaths


## ----Trait names------------------------------------------------------------------------------------------------
traitNames <- column_options$label
traitNames


## ----Code-------------------------------------------------------------------------------------------------------
traitCodes <- column_options$code
traitCodes


## ----Creating the sumstats_meta object--------------------------------------------------------------------------
sumstats_meta <- data.frame(
  name = traitNames
  )

row.names(sumstats_meta) <- traitCodes
sumstats_meta


## ----Sample size------------------------------------------------------------------------------------------------
cN <- column_options$`sample-size`


## ----Indels-----------------------------------------------------------------------------------------------------
keep_indel <- TRUE # or change to FALSE


## ----Standardise GWAS column names------------------------------------------------------------------------------
source(file = "../functions/standardise_GWAS_column_names.R")


## ----Chromosome numbering---------------------------------------------------------------------------------------
source(file = "../functions/parse_SNP_column_as_rs_number.R")


## ----Read in sumstats-------------------------------------------------------------------------------------------
cCode <- traitCodes[1]

cFilePath <- filePaths[1]

cSumstats <- data.table::fread(
  file = cFilePath,
  header = T,
  quote = "\"",
  na.strings = c(
    ".",
    NA,
    "NA",
    ""
    ),
  # For testing read in first 10k rows
  nrows=10000
  )

# Number of rows before touching the data
nRowsRaw <- nrow(cSumstats)
sumstats_meta[cCode,c("n_snp_raw")]<-nRowsRaw #add number of SNPs raw to sumstats_meta (meta data file)

sumstats_meta$n_snp_raw


## ----interpret columns------------------------------------------------------------------------------------------
newNames <- standardise_GWAS_column_names(
  column_names = colnames(cSumstats)
  )

print(newNames)

colnames(cSumstats) <- as.character(newNames$std)


## ----Transform to data table------------------------------------------------------------------------------------
cSumstats <- setDT(cSumstats)
str(cSumstats)


## ----Fix SNP column---------------------------------------------------------------------------------------------
cSumstats.keys <- c('SNP')
cSumstats.keys

cSumstats$SNP <- as.character(cSumstats$SNP)
str(cSumstats$SNP)


## ----Format A1 and A2 to upper case-----------------------------------------------------------------------------
cSumstats$A1 <- toupper(as.character(cSumstats$A1))
str(cSumstats$A1)

cSumstats$A2 <- toupper(as.character(cSumstats$A2))
str(cSumstats$A2)


## ----Parse SNP--------------------------------------------------------------------------------------------------
cSumstats$SNP <- tolower(parse_SNP_column_as_rs_number(cSumstats$SNP))

if(any(colnames(cSumstats)=="CHR")) {
  cSumstats$CHR <- as.integer(parse_CHR_column(as.character(cSumstats$CHR)))
  cSumstats.keys <- c(cSumstats.keys,'CHR')
  str(cSumstats$CHR)
}

if(any(colnames(cSumstats)=="BP")) {
  cSumstats$BP <- as.integer(cSumstats$BP)
  cSumstats.keys <- c(cSumstats.keys,'BP')
  str(cSumstats$BP)
}

cSumstats.keys


## ----Set data table index on selected keys----------------------------------------------------------------------
setkeyv(
  x = cSumstats,
  cols = cSumstats.keys,
  verbose = TRUE
  )


## ----Interpret BGENIE SNP format--------------------------------------------------------------------------------

# Future development



## ----Interpret daner like format--------------------------------------------------------------------------------
if(!any(colnames(cSumstats)=="FRQ") & 
   !any(colnames(cSumstats)=="N_CAS") & 
   !any(colnames(cSumstats)=="N_CON"))
  {
  danerNcas <- colnames(cSumstats)[
    startsWith(colnames(cSumstats),
               prefix = "FRQ_A_")
    ][1]
  
  danerNcon <- colnames(cSumstats)[
    startsWith(colnames(cSumstats),
               prefix = "FRQ_U_")
    ][1]
  
  if(
    !is.na(danerNcas)
    & !is.na(danerNcon))
    {
    
    # Add number of cases and number of controls
    cSumstats$N_CAS <- as.numeric(gsub('FRQ_A_','',danerNcas))
    cSumstats$N_CON <- as.numeric(gsub('FRQ_U_','',danerNcon))
    
    # Add case and control specific effect allele frequencies
    colnames(cSumstats)[colnames(cSumstats) == danerNcas] <- "FRQ_CASES"
    colnames(cSumstats)[colnames(cSumstats) == danerNcon] <- "FRQ_CONTROLS"
    
    # Add total FRQ column
    cSumstats$FRQ<-((cSumstats$FRQ_CASES * cSumstats$N_CAS) + (cSumstats$FRQ_CONTROLS * cSumstats$N_CON))/(cSumstats$N_CAS + cSumstats$N_CON)
    
    colnames(cSumstats)
  }
  }


## ----Convert to numeric columns---------------------------------------------------------------------------------

cSumstats$CHR<-as.numeric(cSumstats$CHR)



## ----add missing N columns--------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="N_CAS") 
  && any(colnames(cSumstats)=="N_CON"))
  {
  
  # Calculate total N from number of cases and number of controls if they are present. Overwrite any specific total N.
  cSumstats$N <- cSumstats$N_CAS + cSumstats$N_CON
  
  } else if(!is.null(cN)) {
    
  cSumstats$N[which(is.na(cSumstats$N))] <- cN #default behaviour: add explicit N only to missing values.
  
  } else if(!any(colnames(cSumstats)=="N")) {
  
    warning("\nNo N column detected!\n")
    
    sumstats_meta[cCode,c("noN")]<-T
    
    }


## ----add BETA column if missing---------------------------------------------------------------------------------
# If no BETA column, but OR column detected
if(
  any(colnames(cSumstats) == "BETA") == 0 
  & any(colnames(cSumstats) == "OR"))
  {
  
  # Calculate BETA from log of OR
  cSumstats$BETA <- log(cSumstats$OR) 
  
  cat("BETA column inserted based on log OR value")
  
  sumstats_meta[cCode,c("Effect")]<-"OR"
  
  } else if(
    any(colnames(cSumstats) == "BETA")
    ) {
    
    cat("BETA column detected")
    
    sumstats_meta[cCode,c("Effect")]<-"BETA"
    
    }


## ----insert Z column if missing---------------------------------------------------------------------------------
# if no Z column
if(
  any(colnames(cSumstats) == "Z") == 0 
  & any(colnames(cSumstats)=="BETA")
  )
  {
  
  # Calculate Z value using qchisq
  cSumstats$Z <- sign(cSumstats$BETA)*sqrt(qchisq(cSumstats$P,1,lower=F))
  
  cat("Z column inserted based on P value and sign(BETA)")
  
  sumstats_meta[cCode,c("Z")]<-"P,sign(BETA)"
  } else if(
    any(colnames(cSumstats)=="SE") 
    & any(colnames(cSumstats)=="BETA")
    ) {
    
    cSumstats$Z <- cSumstats$BETA/cSumstats$SE
    
    cat("Z column inserted based on BETA value and SE")
    
    sumstats_meta[cCode,c("Z")]<-"BETA,SE"
}


## ----insert SE column if missing--------------------------------------------------------------------------------
# if no SE column
if(
  any(colnames(cSumstats) == "SE") == 0
  )
  {
  
  # Take BETA if BETA column exists
  if(any(colnames(cSumstats) == "BETA") == 1)
    
    {
    
    # Calculate SE from BETA and Z
    cSumstats$SE <- abs(cSumstats$BETA/cSumstats$Z)
    
  }
  
  cat("SE column inserted based on BETA/OR and P")
  
  sumstats_meta[cCode, c("SE")] <- "BETA,Z"
  
  }


## ----add meta data for PRS & pathway analysis-------------------------------------------------------------------
# Sufficient columns for PRS
sumstats_meta[cCode, c("enough_columns_PRS")] <-
  any(colnames(cSumstats) == "SNP") &
  any(colnames(cSumstats) == "P") &
  (any(colnames(cSumstats) == "BETA") |
     any(colnames(cSumstats) == "OR") |
     any(colnames(cSumstats) == "Z") ) &
  any(colnames(cSumstats) == "A1") &
  any(colnames(cSumstats) == "A2")

# Sufficient columns for pathway analysis
sumstats_meta[cCode,c("enough_columns_pathway_analysis")] <-
  any(colnames(cSumstats)=="SNP") &
  any(colnames(cSumstats)=="P")

sumstats_meta


## ----clean missing p--------------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="P")
  )
  {
  # If P is na
  cond <- is.na(cSumstats$P)
  
  # Create column n_removed_missing_p in meta file
  sumstats_meta[cCode, c("n_removed_missing_p")] <- sum(cond)
  
  # Subsets sum stats that are not na in the P column
  cSumstats <- cSumstats[which(!cond),]
}

sumstats_meta$n_removed_missing_p


## ----clean missing effect---------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="BETA") # if BETA present
  ){
  
  # If NA for BETA value, remove variant
  cond <- is.na(cSumstats$BETA) 
  sumstats_meta[cCode,c("n_removed_missing_effect")] <- sum(cond)
  cSumstats<-cSumstats[which(!cond),]
  
} else if(
  any(colnames(cSumstats)=="OR")
  ) {
  
  #If NA for OR value, remove variant
  cond <- is.na(cSumstats$OR) 
  sumstats_meta[cCode,c("n_removed_missing_effect")] <- sum(cond)
  cSumstats<-cSumstats[which(!cond),]
  
} else if(
  any(colnames(cSumstats)=="Z")
  ) {
  
  cond <- is.na(cSumstats$Z) #if NA for Z value remove variant
  sumstats_meta[cCode,c("n_removed_missing_effect")] <- sum(cond)
  cSumstats<-cSumstats[which(!cond),]
  
}

sumstats_meta$n_removed_missing_effect


## ----clean duplicate data---------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="SNP") 
  & any(colnames(cSumstats)=="A1") 
  & any(colnames(cSumstats)=="A2")
  ) {
  
  cSumstats.n <- nrow(cSumstats)
  
  cSumstats <- unique(cSumstats, by = c("SNP","A1","A2"))
  
  sumstats_meta[cCode,c("n_removed_duplicates")] <- cSumstats.n-nrow(cSumstats)
  
}

sumstats_meta$n_removed_duplicates


## ----P filter---------------------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="P")
  ){
  
  rm <- (!is.na(cSumstats$P) & (cSumstats$P > 1 | cSumstats$P < 0))
  
  cSumstats <- cSumstats[!rm, ]
  
  cat("Removing", sum(rm), "SNPs with P < 0 or P > 1;", nrow(cSumstats), "remain")
  
  sumstats_meta[cCode,c("n_removed_p")]<-sum(rm)
  
  } else {
    
    cat("Warning: The dataset does not contain a P column to apply the specified filter on.")
    
    }


## ----FRQ filter-------------------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="FRQ")
  ){
  
  frq.filter<-0.005
  
  # If outside of bounds 
  rm <- (!is.na(cSumstats$FRQ) 
         & (cSumstats$FRQ<frq.filter | cSumstats$FRQ>(1-frq.filter)))
  
  cSumstats <- cSumstats[!rm, ]
  
  cat("Removing", sum(rm), "SNPs with MAF <", frq.filter, ";", nrow(cSumstats), "remain")
  
  sumstats_meta[cCode,c("n_removed_frq")] <- sum(rm)
  
  } else {
  
    cat("Warning: The dataset does not contain a FRQ or MAF column to apply the specified filter on.")
    
    }


## ----INFO filter------------------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="INFO")
  ){
  
  rm <- (!is.na(cSumstats$INFO) & cSumstats$INFO<0.6)
  
  cSumstats <- cSumstats[!rm, ]
  
  cat("Removing", sum(rm), "SNPs with INFO <0.6;", nrow(cSumstats), "remain")
  
  sumstats_meta[cCode,c("n_removed_info")] <- sum(rm)
  
  } else {
    
    cat("Warning: The dataset does not contain an INFO column to apply the specified filter on.")
}


## ----OR filter--------------------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="OR")
  ){
  
  rm <- (!is.na(cSumstats$OR) & cSumstats$OR>10000)
  
  cSumstats <- cSumstats[!rm, ]
  
  cat("Removing", sum(rm), "SNPs with OR > 10000;", nrow(cSumstats), "remain")
  
  sumstats_meta[cCode,c("n_removed_or")] <- sum(rm)
  
  } else {
    
    cat("Warning: The dataset does not contain an INFO column to apply the specified filter on.")
    
    }


## ----N filter---------------------------------------------------------------------------------------------------
if(any(names(cSumstats) == 'NEF')){
  # NEF is present in sumstats so variants will be filtered by provided NEF
  N_sd <- sd(cSumstats$NEF)
  N_median <- median(cSumstats$NEF)
  
  cSumstats$N_outlier<-cSumstats$NEF > N_median+(3*N_sd) | cSumstats$NEF < N_median-(3*N_sd)
  
  cat(sum(cSumstats$N_outlier, na.rm=T), "SNPs have reported NEF outside median(N) ± 3SD(N).\n", sep='')
  sumstats_meta[cCode,c("n_outlier_n")] <- sum(cSumstats$N_outlier, na.rm=T)

} else {
  if(any(names(cSumstats) == 'N_CAS') & any(names(cSumstats) == 'N_CON')){
    # NEF isn't present, but N_CAS and N_CON are present, so we will calculate NEF from N_CAS and N_CON.
    cSumstats$NEF_est<-4/(1/cSumstats$N_CAS+1/cSumstats$N_CON)
    
    N_sd <- sd(cSumstats$NEF_est)
    N_median <- median(cSumstats$NEF_est)
    
    cSumstats$N_outlier<-cSumstats$NEF_est > N_median+(3*N_sd) | cSumstats$NEF < N_median-(3*N_sd)
    
    cat(sum(cSumstats$N_outlier, na.rm=T), "SNPs have estimated NEF outside median(N) ± 3SD(N).\n", sep='')
    sumstats_meta[cCode,c("n_outlier_n")] <- sum(cSumstats$N_outlier, na.rm=T)

  } else {
    # NEF nor N_CAS/N_CON columns are present. We will therefore filter by N
    if(length(unique(cSumstats$N)) > 1){
      # There is variation in the N column
      N_sd <- sd(cSumstats$N)
      N_median <- median(cSumstats$N)
      
      cSumstats$N_outlier<-cSumstats$N > N_median+(3*N_sd) | cSumstats$N < N_median-(3*N_sd)
      
      cat(sum(cSumstats$N_outlier, na.rm=T), "SNPs have N outside median(N) ± 3SD(N).\n", sep='')
      sumstats_meta[cCode,c("n_outlier_n")] <- sum(cSumstats$N_outlier, na.rm=T)
    } else {
      # Per variant sample size information is not available. Set column indicating outliers to NA
      cSumstats$N_outlier<-NA
      cat("Per variant sample size is not avialable.\n", sep='')
      sumstats_meta[cCode,c("n_outlier_n")] <- NA

    }
  }
}



## ----SE filter--------------------------------------------------------------------------------------------------
if(
  any(colnames(cSumstats)=="SE")
  ){
  
  # If SE is 0 or NA, SNP gets removed
  rm <- (!is.na(cSumstats$SE) & cSumstats$SE == 0)
  
  cSumstats <- cSumstats[!rm, ]
  
  cat("Removing", sum(rm), "SNPs with SE = 0 or NA;", nrow(cSumstats), "remain")
  
  sumstats_meta[cCode,c("n_removed_se")] <- sum(rm)
  
  } else {
    
    cat("Warning: The dataset does not contain an SE column to apply the specified filter on.")
    
    }


## ---------------------------------------------------------------------------------------------------------------
#if a Zscore already exists, use this to calculate SE from the BETA
if(
  any(colnames(cSumstats)=="Z")
  ){
  
  cSumstats_se_check <- cSumstats %>%
    mutate(
      se_check = 
        BETA / Z
    )
  
  #if there is only SE column, calculate the Zscore first and then the SE from the BETA
} else if(
  any(colnames(cSumstats)=="SE")
  ){
  
  cSumstats_se_check <- cSumstats %>%
    mutate(
      z_check =
        sign(BETA) * abs(qnorm(P/2))
    )
  
  cSumstats_se_check <- cSumstats %>%
    mutate(
      se_check = 
        BETA / z_check
    )
}

#calculate the relative difference between the two SE
mean_rel_diff <- all.equal(target = cSumstats_se_check$se_check, 
                           current = cSumstats$SE)

#print output message
if(
    mean_rel_diff > 0.1
    ){
    
    cat("Mean relative difference in SE is greater than 0.1, please investigate")
      
      } else if(
        mean_rel_diff < 0.1
      ){
        
        cat("Mean relative difference in SE is smaller than 0.1 - passed check")
        
        #remove the dataframe used in the checking process
        rm(cSumstats_se_check)
        
      }



## ----indels filter----------------------------------------------------------------------------------------------
if(keep_indel == FALSE)
  { 
  cSumstats$A1 <- as.character(toupper(cSumstats$A1), c("A", "C", "G", "T"))
  
  cSumstats$A2 <- as.character(toupper(cSumstats$A2), c("A", "C", "G", "T"))
  
  rm <- (is.na(cSumstats$A1) | is.na(cSumstats$A2))
  
  sumstats_meta[cCode,c("n_removed_indels")] <- sum(cond)
  
  cSumstats<-cSumstats[!rm,]
}

sumstats_meta$n_removed_indels


## ----Check and correct for genomic control----------------------------------------------------------------------
if(
  any(colnames(cSumstats) == "SE") == 1)
  {
  if(
    any(colnames(cSumstats) == "BETA") == 1
    )
    {

    cSumstats$P_check <- 2*pnorm(-abs(cSumstats$Z))
    
    if(
      abs(mean(cSumstats$P[!is.na(cSumstats$P_check)]) - mean(cSumstats$P_check[!is.na(cSumstats$P_check)])) > 0.01)
      {
      
      cSumstats$P <- cSumstats$P_check
      
      cSumstats$P_check <- NULL
      
      cat("Genomic control detected. P-value recomputed using BETA and SE.")
      
      sumstats_meta[cCode,c("GC")] <- TRUE
     
    } else {
     
      cat("Genomic control was not detected.")
      
      sumstats_meta[cCode,c("GC")] <- FALSE
    
      cSumstats$P_check <- NULL
    }
    
  }
} else {
  
  cat("SE column is not present, genomic control cannot be detected.")
 
}


## ----Calculate lambda GC----------------------------------------------------------------------------------------
sumstats_meta[cCode,c("lambdaGC")] <- median(cSumstats$Z^2)/qchisq(0.5,1)

cat(
  "\nThe genomic inflation factor (lambda GC) was calculated as:",
  round(sumstats_meta[cCode,c("lambdaGC")],
        digits = 2)
  )


## ----Read in functions for manipulating allele codes------------------------------------------------------------
source(file = "../functions/snp_allele.R")


## ----insert IUPAC codes-----------------------------------------------------------------------------------------
# Insert IUPAC codes into target
cSumstats$IUPAC<-snp_iupac(cSumstats$A1, cSumstats$A2)



## ----define super population of GWAS sample---------------------------------------------------------------------
# Define super population
super_pop <- column_options$population



## ----check presence of CHR, BP and RSID information-------------------------------------------------------------
# Initially target should be merged with the reference based on CHR, BP and IUPAC. RSIDs will be inserted for all SNPs, and reference MAF will be inserted for all non-ambiguous SNPs if the super population is specified.
# If CHR and BP are unavailable, but RSID is present, target will be merged with the reference by RSID and IUPAC. CHR and BP will be inserted for all SNPs, and reference MAF will be inserted for all non-ambiguous SNPs if the super population is specified.

# Check whether chromosome and base pair position information is present
chr_bp_avail<-sum(c('CHR','ORIGBP') %in% names(cSumstats)) == 2 

# Check whether RSIDs are available for majority of SNPs in GWAS
rsid_avail<-(sum(grepl('rs', cSumstats$SNP)) > 0.9*length(cSumstats$SNP))



## ----merge with reference---------------------------------------------------------------------------------------
if(chr_bp_avail){
  ###
  # Determine build
  ###
  # Read in random chromosome of reference data present in GWAS
  i<-unique(cSumstats$CHR)[1]
  
  ref<-readRDS(file = paste0('/scratch/groups/gwas_sumstats/1kg_ref/1kg_ref_chr',i,'.rds'))
  
  # Check target-ref condordance of BP across builds
  ref$CHR<-as.numeric(ref$CHR)
  matched<-list()
  matched[['GRCh37']]<-merge(cSumstats, ref, by.x=c('CHR','ORIGBP','IUPAC'), by.y=c('CHR','BP_GRCh37','IUPAC'))
  matched[['GRCh38']]<-merge(cSumstats, ref, by.x=c('CHR','ORIGBP','IUPAC'), by.y=c('CHR','BP_GRCh38','IUPAC'))
  
  cat('GRCh37 match: ',round(nrow(matched[['GRCh37']])/sum(cSumstats$CHR == i)*100, 2),'%\n',sep='')
  cat('GRCh38 match: ',round(nrow(matched[['GRCh38']])/sum(cSumstats$CHR == i)*100, 2),'%\n',sep='')

  target_build<-NA
  if((nrow(matched[['GRCh37']])/sum(cSumstats$CHR == i)) > 0.7 & (nrow(matched[['GRCh37']])/sum(cSumstats$CHR == i)) > (nrow(matched[['GRCh38']])/sum(cSumstats$CHR == i))){
    target_build<-'GRCh37'
  }
  
  if((nrow(matched[['GRCh38']])/sum(cSumstats$CHR == i)) > 0.7 & (nrow(matched[['GRCh38']])/sum(cSumstats$CHR == i)) > (nrow(matched[['GRCh37']])/sum(cSumstats$CHR == i))){
    target_build<-'GRCh38'
  }
  
  rm(matched,ref)
  
  if(!is.na(target_build)){
    # Build detected, so can continue reference harmonisation
    # Run per chromosome to reduce memory requirements
    chrs<-unique(cSumstats$CHR)

    cSumstats_matched<-NULL
    cSumstats_unmatched<-NULL
    for(i in chrs){
      print(i)
      
      # Read reference data
      tmp<-readRDS(file = paste0('/scratch/groups/gwas_sumstats/1kg_ref/1kg_ref_chr',i,'.rds'))
      
      # Subset relevent data
      if(super_pop %in% c('AFR','AMR','EAS','EUR','SAS')){
        tmp<-tmp[,c("CHR","SNP","BP_GRCh37","BP_GRCh38","A1","A2","IUPAC",paste0('FREQ_',super_pop)),with=F]
        names(tmp)[names(tmp) == 'BP_GRCh37']<-'REF.BP_GRCh37'
        names(tmp)[names(tmp) == 'BP_GRCh38']<-'REF.BP_GRCh38'
        names(tmp)[names(tmp) == paste0('FREQ_',super_pop)]<-'REF.FRQ'
      } else {
        tmp<-tmp[,c("CHR","SNP","BP_GRCh37","BP_GRCh38","A1","A2","IUPAC"),with=F]
      }

      tmp<-tmp[nchar(tmp$A1) == 1 & nchar(tmp$A2) == 1,]
      tmp$CHR<-NULL
      
      names(tmp)[names(tmp) == 'SNP']<-'REF.SNP'

      # Merge target and reference by BP
      cSumstats_chr<-cSumstats[cSumstats$CHR == i,]
      ref_target<-merge(cSumstats_chr, tmp, by.x='ORIGBP', by.y=paste0('REF.BP_',target_build))

      
      # Identify SNPs that are opposite strands
      flipped<-ref_target[(ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'Y') | 
                            (ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'R') | 
                            (ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'M') |
                            (ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'K'),]
      

      # Change target alleles to compliment for flipped variants
      flipped$A1.x<-snp_allele_comp(flipped$A1.x)
      flipped$A2.x<-snp_allele_comp(flipped$A2.x)
      
      # Update IUPAC codes
      flipped$IUPAC.x<-snp_iupac(flipped$A1.x, flipped$A2.x)

      
      # Identify SNPs that have matched alleles
      matched<-ref_target[ref_target$IUPAC.x == ref_target$IUPAC.y,]
      matched<-rbind(matched, flipped)
      
      # Flip REF.FRQ if alleles are swapped
      matched$REF.FRQ[matched$A1.x != matched$A1.y]<-1-matched$REF.FRQ[matched$A1.x != matched$A1.y]
      
      # Remove REF.FRQ for ambiguous SNPs
      if(!is.na(super_pop)){
        matched$REF.FRQ[matched$IUPAC.x %in% c('W','S')]<-NA
      } else {
        matched$REF.FRQ<-NA
      }
      
      # Retain reference CHR, BP, but target A1, and A2 information
      matched$A1<-matched$A1.x
      matched$A1.y<-NULL
      matched$A1.x<-NULL
      matched$A2<-matched$A2.x
      matched$A2.y<-NULL
      matched$A2.x<-NULL
      matched$IUPAC<-matched$IUPAC.x
      matched$IUPAC.y<-NULL
      matched$REF.CHR<-matched$CHR
      matched[[paste0('REF.BP_',target_build)]]<-matched$ORIGBP

      # Identify SNPs that are unmatched
      unmatched_chr<-cSumstats_chr[!(paste0(cSumstats_chr$ORIGBP,':',cSumstats_chr$IUPAC) %in% paste0(matched$ORIGBP,':',matched$IUPAC.x)),]
      matched$IUPAC.x<-NULL

      cSumstats_matched<-rbind(cSumstats_matched, matched)
      cSumstats_unmatched<-rbind(cSumstats_unmatched, unmatched_chr)
    }
    
    # Reinsert variants in GWAS that could not be matched to the reference
    # These unmatched variants include those with CHR:BP:IUPAC not in the reference, and INDELS.
    cSumstats_unmatched$REF.CHR<-NA
    cSumstats_unmatched$REF.BP_GRCh37<-NA
    cSumstats_unmatched$REF.BP_GRCh38<-NA
    cSumstats_unmatched$REF.SNP<-NA
    cSumstats_unmatched$REF.FRQ<-NA

    cSumstats_harm<-rbind(cSumstats_matched, cSumstats_unmatched)
  }
} else {
  if(rsid_avail){
    # Run per chromosome to reduce memory requirements
    chrs<-c(1:22,'X','Y','MT')
  
    cSumstats_matched<-NULL
    for(i in chrs){
      print(i)
      
      # Read reference data
      tmp<-readRDS(file = paste0('/scratch/groups/gwas_sumstats/1kg_ref/1kg_ref_chr',i,'.rds'))
      
      # Subset relevent data
      if(super_pop %in% c('AFR','AMR','EAS','EUR','SAS')){
        tmp<-tmp[,c("CHR","SNP","BP_GRCh37","BP_GRCh38","A1","A2","IUPAC",paste0('FREQ_',super_pop)),with=F]
        names(tmp)[names(tmp) == 'BP_GRCh37']<-'REF.BP_GRCh37'
        names(tmp)[names(tmp) == 'BP_GRCh38']<-'REF.BP_GRCh38'
        names(tmp)[names(tmp) == paste0('FREQ_',super_pop)]<-'REF.FRQ'
      } else {
        tmp<-tmp[,c("CHR","SNP","BP_GRCh37","BP_GRCh38","A1","A2","IUPAC"),with=F]
      }

      tmp<-tmp[nchar(tmp$A1) == 1 & nchar(tmp$A2) == 1,]
      tmp$CHR<-NULL
      
      # Merge target and reference by SNP ID
      ref_target<-merge(cSumstats, tmp, by='SNP')
      
            # Identify SNPs that are opposite strands
      flipped<-ref_target[(ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'Y') | 
                            (ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'R') | 
                            (ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'M') |
                            (ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'K'),]
      
      # Change target alleles to compliment for flipped variants
      flipped$A1.x<-snp_allele_comp(flipped$A1.x)
      flipped$A2.x<-snp_allele_comp(flipped$A2.x)
      
      # Update IUPAC codes
      flipped$IUPAC.x<-snp_iupac(flipped$A1.x, flipped$A2.x)
      
      # Identify SNPs that have matched alleles
      matched<-ref_target[ref_target$IUPAC.x == ref_target$IUPAC.y,]
      matched<-rbind(matched, flipped)
      
      # Flip REF.FRQ if alleles are swapped
      matched$REF.FRQ[matched$A1.x != matched$A1.y]<-1-matched$REF.FRQ[matched$A1.x != matched$A1.y]
      
      # Remove REF.FRQ for ambiguous SNPs
      if(!is.na(super_pop)){
        matched$REF.FRQ[matched$IUPAC.x %in% c('W','S')]<-NA
      } else {
        matched$REF.FRQ<-NA
      }

      # Retain reference CHR, BP, but target A1, and A2 information
      matched$A1<-matched$A1.x
      matched$A1.y<-NULL
      matched$A1.x<-NULL
      matched$A2<-matched$A2.x
      matched$A2.y<-NULL
      matched$A2.x<-NULL
      matched$IUPAC<-matched$IUPAC.x
      matched$IUPAC.y<-NULL
      matched$IUPAC.x<-NULL
      matched$REF.SNP<-matched$SNP
  
      cSumstats_matched_chr<-rbind(matched)
      cSumstats_matched<-rbind(cSumstats_matched, cSumstats_matched_chr)
    }
    
    # Identify SNPs that are unmatched
    unmatched<-cSumstats[!(cSumstats$SNP %in% cSumstats_matched$SNP),]
    
    # Reinsert variants in GWAS that could not be matched to the reference
    # These unmatched variants include those with RSIDs not in the reference, any ambiguous SNPs, and INDELS.
    unmatched$REF.SNP<-NA
    unmatched$REF.FRQ<-NA
    unmatched$REF.CHR<-NA
    unmatched$REF.BP_GRCh37<-NA
    unmatched$REF.BP_GRCh38<-NA
  
    cSumstats_harm<-rbind(cSumstats_matched, unmatched)

  }
}



## ----allele Frequency columns-----------------------------------------------------------------------------------
# Check FRQ column is valid
if('FRQ' %in% names(cSumstats_harm)){
  if(any(cSumstats_harm$FRQ>1) || any(cSumstats_harm$FRQ<0)) {
    
    stop(
        paste0('\nThere are FRQ values larger than 1 (',sum(cSumstats_harm$FRQ>1),') or less than 0 (',sum(cSumstats_harm$FRQ<0),') which is outside of the possible FRQ range.')
    )
  }
}



## ----allele Frequency comparison with reference-----------------------------------------------------------------
# Check FRQ column is valid
if('FRQ' %in% names(cSumstats_harm) & !is.na(super_pop)){
  cSumstats_harm$diff<-abs(cSumstats_harm$FRQ-cSumstats_harm$REF.FRQ)
  cSumstats_harm$FRQ_bad<-F
  cSumstats_harm$FRQ_bad[is.na(cSumstats_harm$diff) | cSumstats_harm$diff > 0.2]<-T
  cSumstats_harm$diff<-NULL
}



## ----save cleaned file as .gz file in cleaned folder------------------------------------------------------------

fwrite(cSumstats_harm, paste0(column_options$output,"/cleaned/", cCode, ".gz"), sep='\t', na='NA', quote=F)



## ----create log file and save in log file folder----------------------------------------------------------------
sumstats_meta %>%
  write_tsv(
    path = paste0(column_options$output,"/cleaned_logs/", cCode, "_log.txt"),
    col_names = TRUE
    )

