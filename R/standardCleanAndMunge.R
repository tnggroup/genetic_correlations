#you will need to have these packages installed
# Standardised pipeline as an R function


#needs the shru package
#devtools::install_github("johanzvrskovec/shru")


# #Test
# filePaths = c(
#  file.path("/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/raw/","PGC3_ED_2022","daner_BENARROW.gz"),
#  file.path("/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/raw/","Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
# )
# traitCodes = c("BEN","BMI")
# traitNames = c("Binge Eating (Narrow)","BMI 2018")
# referenceFilePath = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/variant_lists/combined.hm3_1kg.snplist.vanilla.jz2020.gz"
# n_threads=5
# keep_indel=T
# maf_filter=0.01
# info_filter=0.6
# serviceAccountTokenPath=normalizePath("/scratch/prj/gwas_sumstats/tngpipeline/tngpipeline-8130dbd7d58a.json",mustWork = T)
# groupFolderPath = normalizePath("/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/gwas_sumstats_test",mustWork = T)
# munge="opmunge"
# N = c(830917,681275)
# ancestrySetting =c("EUR")

#
# filePaths = c(
#   file.path("new","bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"),
#   file.path("new","ADHD2023.gz")
# )
# traitCodes = NA_character_
# traitNames = c("BMI", "ADHD")
# referenceFilePath = "variant_lists/combined.hm3_1kg.snplist.vanilla.jz2020.gz"
# n_threads=5
# keep_indel=T
# maf_filter=0.01
# info_filter=0.6
# serviceAccountTokenPath=normalizePath("/scratch/prj/gwas_sumstats/tngpipeline/tngpipeline-8130dbd7d58a.json",mustWork = T)
# groupFolderPath = normalizePath("/scratch/prj/gwas_sumstats",mustWork = T)
# munge="opmunge"
# N = c(806834,NA)
# ancestrySetting =c("EUR")

standardPipelineCleaningAndMunging <- function(
    filePaths,
    traitCodes=NA_character_, #set explicit code(s) here
    traitNames,
    referenceFilePath,
    n_threads=5,
    keep_indel=T,
    maf_filter=0.01,
    info_filter=0.6,
    serviceAccountTokenPath=normalizePath("/scratch/prj/gwas_sumstats/tngpipeline/tngpipeline-8130dbd7d58a.json",mustWork = T),
    groupFolderPath = normalizePath("/scratch/prj/gwas_sumstats",mustWork = T),
    munge="opmunge", #alt supermunge
    filter.mhc=NULL, #can be either 37 or 38 for filtering the MHC region according to either grch37 or grch38
    N=NA_integer_,
    ancestrySetting
    ){


#set up metadata df
  sumstats_meta <- data.frame(
    name = traitNames
  )

  #row.names(sumstats_meta) <- traitCodes
  sumstats_meta$code<-traitCodes
  sumstats_meta$path_orig<-filePaths
  sumstats_meta$N<-N
  sumstats_meta$ancestrySetting<-ancestrySetting
  sumstats_meta$n_case<-NA_character_
  sumstats_meta$n_control<-NA_character_


  #configure gs4 for non-browser login #https://gargle.r-lib.org/articles/non-interactive-auth.html
  drive_auth(email = "johan.kallberg_zvrskovec@kcl.ac.uk", path = serviceAccountTokenPath)
  #drive_auth(email = "jane_doe@example.com") # gets a suitably scoped token
  # and stashes for googledrive use
  gs4_auth(token = drive_token())            # registers token with googlesheets4

  #used by supermunge - may be harmonised later so all steps use the same format of variant lists (summary level reference panel).
  varlist<-fread(file = referenceFilePath, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = n_threads, showProgress = F)

  #testTrait<-readFile(filePath = filePaths[iTrait])

  for(iTrait in 1:filePaths){
    #iTrait <-1 #test

    #read trait metadata if exists
    metaFilePath <- paste0(filePaths[iTrait],".txt")
    #metaFilePath <- "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/raw/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz.txt"
    if(file.exists(metaFilePath)){
      nMetaList <- readMetadata(sumstats_meta_list = unlist(sumstats_meta[iTrait,]), filePath = metaFilePath)
      nMetaList.names<-names(nMetaList)

      sumstats_meta[iTrait,nMetaList.names]<-nMetaList
    }

    #set new code using the spreadsheet metadata
    if(is.na(sumstats_meta[iTrait,]$code)){
      cSort<-sumstats_meta[iTrait,]$sort
      if(is.na(cSort)) cSort<-"UNST" #unsorted
      nCode <- assignCodeFromSortAndSpreadsheet(sort = cSort)
      sumstats_meta[iTrait,c("code")]<-nCode
    }


    #edit metadata after reading the file-based metadata
    if(is.na(sumstats_meta[iTrait,]$N)) sumstats_meta[iTrait,c("N")] <- sum(as.integer(sumstats_meta[iTrait,c("n_case")]), as.integer(sumstats_meta[iTrait,c("n_control")]),na.rm = T)

    sumstats_meta[iTrait,c("dependent_variable")]<-ifelse(!is.na(sumstats_meta[iTrait,]$n_case) & !is.na(sumstats_meta[iTrait,]$n_control), "binary", "continuous")

    #clean using shru::supermunge - implements most of the cleaning and parsing steps in the previous implementation, plus some additions and fixes. this may be harmonised later to either increase or reduce the dependency on the shru package.
    ref_df_arg <-NULL
    if(munge=="supermunge") ref_df_arg<-varlist

    smungeResults <- shru::supermunge(
      filePaths = filePaths[iTrait],
      ref_df = ref_df_arg,
      traitNames = traitCodes[iTrait],
      ancestrySetting = ancestrySetting[iTrait],
      N = N[iTrait],
      keepIndel = keep_indel,
      writeOutput = F,
      filter.maf = maf_filter,
      filter.info = info_filter,
      lossless = T)

    procResults <- standardPipelineExplicitSumstatProcessing(
      cSumstats = smungeResults$last,
      sumstats_meta = sumstats_meta,
      cCode = traitCodes[iTrait],
      super_pop = ancestrySetting[iTrait],
      munge = (munge=="opmunge")
      )

    cSumstats <- procResults$procResults
    sumstats_meta <-procResults$sumstats_meta

    if(!is.na(filter.mhc)){
      smungeResults2 <- shru::supermunge(
        list_df = list(procResults),
        filter.mhc=38,
        process = F
        )
      cSumstats <- smungeResults2$last
    }



  }

}

readMetadata <- function(sumstats_meta_list,filePath){
  #filePath <- "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/raw/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz.txt"
  f <- readLines(con = file(filePath,open="r"),encoding = "UTF-8")
  for(iLine in 1:length(f)){
    if(nchar(f[iLine])>2){
      #iLine<-29
      #strsplit(x = f[iLine], split = "\\s")
      m<-gregexpr(pattern = "\\s*(\\w+)\\s*",text = f[iLine])
      k<-NA
      v<-NA
      if(length(m)>0){
        if(length(m[[1]])>0) k<- trimws(substr(x =f[iLine], start = m[[1]][[1]],  stop = m[[1]][[1]]+(attr(x = m[[1]], which = "match.length")[[1]]-1)),which = "both")

        v <- trimws(substr(f[iLine],start = m[[1]][[1]]+(attr(x = m[[1]], which = "match.length")[[1]]),stop = nchar(f[iLine])), which = "both")
      }

      if(!is.na(k) & !is.na(v)) sumstats_meta_list[k]<-v

    }
  }

  return(sumstats_meta_list)

}

standardPipelineExplicitSumstatProcessing <- function(cSumstats, sumstats_meta, cCode, super_pop, munge=T, refPath='/scratch/prj/gwas_sumstats/1kg_ref'){
  # cSumstats <- smungeResults$last
  # cCode <- traitCodes[iTrait]
  # super_pop <- "EUR"

  # Meta data about sufficient columns for analyses
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


  #explicit hard coded FRQ filter
  if(
    any(colnames(cSumstats)=="FRQ")
  ){

    frq.filter<-0.005

    # If outside of bounds
    rm <- (!is.na(cSumstats$FRQ)
           & (cSumstats$FRQ<frq.filter | cSumstats$FRQ>(1-frq.filter)))

    cSumstats <- cSumstats[!rm, ]

    cat("Removing", sum(rm), "SNPs with FRQ <", frq.filter, " or (1-",frq.filter,");", nrow(cSumstats), "remain")

    sumstats_meta[cCode,c("n_removed_frq")] <- sum(rm)

  } else {

    cat("Warning: The dataset does not contain a FRQ or MAF column to apply the specified filter on.")

  }


  #explicit OR filter
  if(
    any(colnames(cSumstats)=="OR")
  ){

    rm <- (!is.na(cSumstats$OR) & cSumstats$OR>10000)

    cSumstats <- cSumstats[!rm, ]

    cat("Removing", sum(rm), "SNPs with OR > 10000;", nrow(cSumstats), "remain")

    sumstats_meta[cCode,c("n_removed_or")] <- sum(rm)

  } else {

    cat("Warning: The dataset does not contain an OR column to apply the specified filter on.")

  }


  #explicit N statistics
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

  #explicit SE filter
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

  #explicit SE checks
  #if a Zscore already exists, use this to calculate SE from the BETA
  if(any(colnames(cSumstats)=="BETA") & any(colnames(cSumstats)=="SE")){
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

    cSumstats_se_check <- as.data.frame(cSumstats_se_check)
    #calculate the relative difference between the two SE
    mean_rel_diff <- mean(cSumstats_se_check$se_check - cSumstats$SE, na.rm = TRUE)


    print(mean_rel_diff)
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
  }


  # Check and correct for Genomic Control
  # +++JZ: Would this rather detect logistic SEs and not corresponding ln(OR) effect?
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


  # Genomic Control calculation: Lambda
  sumstats_meta[cCode,c("lambdaGC")] <- median(cSumstats$Z^2,na.rm = T)/qchisq(0.5,1)

  cat(
    "\nThe genomic inflation factor (lambda GC) was calculated as:",
    round(sumstats_meta[cCode,c("lambdaGC")],
          digits = 2)
  )

  # Harmonise with the reference
  # +++JZ: This is munging IMO, so I will regulate it with a munge argument
  if(munge){

    # Insert IUPAC codes into target
    cSumstats$IUPAC<-snp_iupac(cSumstats$A1, cSumstats$A2)

    # Check whether chromosome and base pair position information is present
    chr_bp_avail<-sum(c('CHR','BP') %in% names(cSumstats)) == 2

    # Check whether RSIDs are available for majority of SNPs in GWAS
    rsid_avail<-(sum(grepl('rs', cSumstats$SNP)) > 0.9*length(cSumstats$SNP))

    if(chr_bp_avail){
      ###
      # Determine build
      ###
      # Read in random chromosome of reference data present in GWAS
      i<-unique(cSumstats$CHR)[1]

      ref<-readRDS(file = file.path(refPath,'1kg_ref_chr',i,'.rds'))
      #ref<-readRDS(file = paste0('/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/gwas_sumstats_test/1kg_ref/1kg_ref_chr4.rds')) #test

      # Check target-ref condordance of BP across builds
      ref$CHR<-as.numeric(ref$CHR)
      matched<-list()
      matched[['GRCh37']]<-merge(cSumstats, ref, by.x=c('CHR','BP','IUPAC'), by.y=c('CHR','BP_GRCh37','IUPAC'))
      matched[['GRCh38']]<-merge(cSumstats, ref, by.x=c('CHR','BP','IUPAC'), by.y=c('CHR','BP_GRCh38','IUPAC'))

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
          #i<-4 #test
          print(i)

          # Read reference data
          tmp<-readRDS(file = file.path(refPath,'1kg_ref_chr',i,'.rds'))
          #tmp<-readRDS(file = '/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/gwas_sumstats_test/1kg_ref/1kg_ref_chr4.rds') #test

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
          ref_target<-merge(cSumstats_chr, tmp, by.x='BP', by.y=paste0('REF.BP_',target_build))


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
          matched[[paste0('REF.BP_',target_build)]]<-matched$BP

          # Identify SNPs that are unmatched
          unmatched_chr<-cSumstats_chr[!(paste0(cSumstats_chr$BP,':',cSumstats_chr$IUPAC) %in% paste0(matched$BP,':',matched$IUPAC.x)),]
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

        cSumstats_harm<-rbindlist(l = list(cSumstats_matched, cSumstats_unmatched), use.names = T, fill = T)
      }
    } else {
      if(rsid_avail){
        # Run per chromosome to reduce memory requirements
        chrs<-c(1:25)

        cSumstats_matched<-NULL
        for(i in chrs){
          #i<-4
          print(i)

          tmp<-readRDS(file = file.path(refPath,'1kg_ref_chr',i,'.rds'))
          #tmp<-readRDS(file = '/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/gwas_sumstats_test/1kg_ref/1kg_ref_chr4.rds') #test

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
          cSumstats_matched<-rbindlist(l = list(cSumstats_matched, cSumstats_matched_chr), use.names = T, fill = T)
        }

        cSumstats_matched[,BP:=REF.BP_GRCh38] #+++JZ: Harmonise BP to reference BP when merging on SNP. Allows for coordinate based filters later.

        # Identify SNPs that are unmatched
        unmatched<-cSumstats[!(cSumstats$SNP %in% cSumstats_matched$SNP),]

        # Reinsert variants in GWAS that could not be matched to the reference
        # These unmatched variants include those with RSIDs not in the reference, any ambiguous SNPs, and INDELS.
        unmatched$REF.SNP<-NA
        unmatched$REF.FRQ<-NA
        unmatched$REF.CHR<-NA
        unmatched$REF.BP_GRCh37<-NA
        unmatched$REF.BP_GRCh38<-NA

        cSumstats_harm<-rbindlist(l = list(cSumstats_matched, unmatched), use.names = T, fill = T)

      }
    }

    if(!is.null(cSumstats_matched)) cSumstats<-cSumstats_matched

    # Check FRQ column is valid
    if('FRQ' %in% names(cSumstats) & !is.na(super_pop)){
      cSumstats$diff<-abs(cSumstats$FRQ-cSumstats$REF.FRQ)
      cSumstats$FRQ_bad<-F
      cSumstats$FRQ_bad[is.na(cSumstats$diff) | cSumstats$diff > 0.2]<-T
      cSumstats$diff<-NULL
    }

    sumstats_meta[cCode,c("frq_bad")] <- sum(cSumstats$FRQ_bad)
    sumstats_meta[cCode,c("frq_bad")]

  }




  return(list(cSumstats=cSumstats,sumstats_meta=sumstats_meta))

}

assignCodeFromSortAndSpreadsheet <- function(
    sort,
    sheetLink="https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"
    ){

  #sheetLink <- "https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"

  currentSheet <- as.data.frame(read_sheet(ss = sheetLink, col_names = T, range="SGDP_GWASLIST_EDITTHIS"))
  sorts<-currentSheet[grep(pattern = paste0("^",toupper(sort)), x = currentSheet$code),]
  if(nrow(sorts)>0){
    indexesLengths<-regexec(pattern = "\\d+", text=sorts$code)
    matches<-regmatches(sorts$code,indexesLengths)
    nums<-as.integer(unique(unlist(matches)))
    return(paste0(sort,shru::padStringLeft(s = paste0("",(max(nums)+1)), padding = "0", targetLength = 2)))
  } else return(paste0(sort,"01"))

}



updateSpreadsheet <- function(sheetLink,sumstats_meta){
  #sheetLink <- "https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"

  currentSheet <- as.data.frame(read_sheet(ss = sheetLink, col_names = T, range="SGDP_GWASLIST_EDITTHIS"))
  cols<-colnames(currentSheet)

  if(nrow(currentSheet)>0) currentSheet$x_row<-1:nrow(currentSheet)
  setDT(currentSheet)
  setkeyv(currentSheet,cols = c("code"))

  # sumstats_meta_dummy <- data.frame(
  #   code=c("ANXI04","DEPR05","NEUR02"),
  #   name=c("Generalised anxiety symptoms","MDD, narrow","Neuroticism"),
  #   ancestry=c("EUR","EUR","EUR"),
  #   n_ca=c(19012,58113,449484),
  #   n_co=c(58113,25632,NA_integer_)
  #   )
  # sumstats_meta_dummy$n_tot<-sumstats_meta_dummy$n_ca + ifelse(is.na(sumstats_meta_dummy$n_co),0,sumstats_meta_dummy$n_co)

  dfToInsert<-as.data.frame(matrix(data=NA,nrow = 0, ncol = length(colnames(currentSheet))))
  colnames(dfToInsert)<-colnames(currentSheet)

  setDT(sumstats_meta)
  sumstats_meta<-sumstats_meta[,.(
    name,
    code,
    N,
    n_case,
    n_control,
    doi,
    pmid,
    permission,
    phenotype_type,
    ancestry,
    sex,
    phenotype,
    category,
    notes,
    dependent_variable
  )]
  setkeyv(sumstats_meta,cols = c("code"))


  dfToInsert<-rbindlist(list(
    dfToInsert,
    sumstats_meta
  ), fill = T)


  dfToInsert <- dfToInsert[,..cols]
  for(iRec in 1:nrow(dfToInsert)){
    #iRec<-1
    cRec <- dfToInsert[iRec,]
    if(!is.na(cRec$code)){
      if(any(currentSheet$code==cRec$code)){
        #inactivated update for now
        # toUpdate<-currentSheet[currentSheet$code==cRec$code,]$x_row
        # range_write(ss = sheetLink,data = cRec,cell_rows(toUpdate[[1]]:toUpdate[[1]]))
      } else {
        sheet_append(ss = sheetLink,data = cRec, sheet="SGDP_GWASLIST_EDITTHIS") #this does not recognise some columns
      }
    }
  }
}




