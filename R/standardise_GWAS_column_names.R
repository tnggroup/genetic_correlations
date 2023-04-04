#Function to standardise column names
#Johan Zvrskovec, 2021
#The function interprets the columns of the original raw GWAS summary statistics and replaces them with standard column names. The standard column names are the first item in the vector:
#   SNP	CHR	ORIGBP	A1	A2	FREQ	FREQ_CASES	P	INFO	OR/BETA/Z	SE	N	Ncas	Ncon
# We have chosen ORIGBP as output name, because the summary statistics are based on different genome builds which means SNPs will be located at different BPs
# To avoid follow-up software using the BP column (cave: different genome builds), we named it ORIGBP.
# The function stops if the essential columns SNP, A1, A2, P, BETA or OR or Z 
# 
# _Future_
# +++CH: FREQ_CASES is missing when reading in, NEF is missing when reading in
# +++CH: Ollie has code that detects the genome build of sumstats; we could add the chunk for lifting over. (not top priority).
# +++CH: Johan, can you try to find the most recent version of Helena's script on the virtual machine?
standardise_GWAS_column_names <- function(column_names, stop_on_missing_essential_column = T,
                                          c.SNP = c(
                                            "SNP",
                                            "PREDICTOR",
                                            "SNPID",
                                            "MARKERNAME",
                                            "MARKER_NAME",
                                            "SNPTESTID",
                                            "ID_DBSNP49",
                                            "RSID",
                                            "ID",
                                            "RS_NUMBER",
                                            "MARKER",
                                            "RS",
                                            "RSNUMBER",
                                            "RS_NUMBERS",
                                            "SNP.NAME",
                                            "SNP ID",
                                            "SNP_ID",
                                            "LOCATIONALID",
                                            "ASSAY_NAME"
                                          ),
                                          c.A1 = c(
                                            "A1",
                                            "ALLELE1",
                                            "ALLELE_1",
                                            "INC_ALLELE",
                                            "EA",
                                            "A1_EFFECT",
                                            "REF",
                                            "EFFECT_ALLELE",
                                            "RISK_ALLELE",
                                            "EFFECTALLELE",
                                            "EFFECT_ALL",
                                            "REFERENCE_ALLELE",
                                            "REF_ALLELE",
                                            "REFERENCEALLELE",
                                            "EA",
                                            "ALLELE_1",
                                            "INC_ALLELE",
                                            "ALLELE1",
                                            "A",
                                            "A_1",
                                            "CODED_ALLELE",
                                            "TESTED_ALLELE"
                                          ),
                                          c.A2 = c(
                                            "A2",
                                            "ALLELE2",
                                            "ALLELE_2",
                                            "OTHER_ALLELE",
                                            "NON_EFFECT_ALLELE",
                                            "DEC_ALLELE",
                                            "OA",
                                            "NEA",
                                            "ALT",
                                            "A2_OTHER",
                                            "NONREF_ALLELE",
                                            "NEFFECT_ALLELE",
                                            "NEFFECTALLELE",
                                            "NONEFFECT_ALLELE",
                                            "OTHER_ALL",
                                            "OTHERALLELE",
                                            "NONEFFECTALLELE",
                                            "ALLELE0",
                                            "ALLELE_0",
                                            "ALT_ALLELE",
                                            "A_0",
                                            "NONCODED_ALLELE"
                                          ),
                                          #c.EFFECT = c(
                                          #"EFFECT",
                                          #"OR",
                                          #"B",
                                          #"BETA",
                                          #"LOG_ODDS",
                                          #"EFFECTS",
                                          #"SIGNED_SUMSTAT",
                                          #"EST"
                                          #),
                                          c.BETA = c(
                                            "BETA",
                                            "B",
                                            "EFFECT_BETA",
                                            "EFFECT",
                                            "EFFECTS",
                                            "SIGNED_SUMSTAT",
                                            "EST",
                                            "GWAS_BETA",
                                            "EFFECT_A1",
                                            "EFFECTA1",
                                            "EFFECT_NW"
                                          ),
                                          c.OR = c(
                                            "OR",
                                            "LOG_ODDS",
                                            "OR",
                                            "ODDS-RATIO",
                                            "ODDS_RATIO",
                                            "ODDSRATIO",
                                            "OR(MINALLELE)",
                                            "OR.LOGISTIC",
                                            "OR_RAN",
                                            "OR(A1)"
                                          ),
                                          c.SE = c(
                                            "SE",
                                            "STDER",
                                            "STDERR",
                                            "STD",
                                            "STANDARD_ERROR",
                                            "OR_SE",
                                            "STANDARDERROR",
                                            "STDERR_NW",
                                            "META.SE",
                                            "SE_DGC",
                                            "SE.2GC"
                                          ),
                                          c.Z = c(
                                            "Z",
                                            "ZSCORE",
                                            "Z-SCORE",
                                            "ZSTAT",
                                            "ZSTATISTIC",
                                            "GC_ZSCORE",
                                            "BETAZSCALE"
                                          ),
                                          c.INFO = c(
                                            "INFO",
                                            "IMPINFO",
                                            "IMPQUALITY",
                                            "INFO.PLINK",
                                            "INFO_UKBB",
                                            "INFO_UKB"
                                          ),
                                          c.P = c(
                                            "P",
                                            "PVALUE",
                                            "PVAL",
                                            "P_VALUE",
                                            "GC_PVALUE",
                                            "WALD_P",
                                            "P.VAL",
                                            "GWAS_P",
                                            "P-VALUE",
                                            "P-VAL",
                                            "FREQUENTIST_ADD_PVALUE",
                                            "P.VALUE",
                                            "P_VAL",
                                            "SCAN-P",
                                            "P.LMM",
                                            "META.PVAL",
                                            "P_RAN",
                                            "P.ADD",
                                            "P_BOLT_LMM"
                                          ),
                                          c.N = c(
                                            "N",
                                            "WEIGHT",
                                            "NCOMPLETESAMPLES",
                                            "TOTALSAMPLESIZE",
                                            "TOTALN",
                                            "TOTAL_N",
                                            "N_COMPLETE_SAMPLES",
                                            "N_TOTAL",
                                            "N_SAMPLES",
                                            "N_ANALYZED",
                                            "NSAMPLES",
                                            "SAMPLESIZE",
                                            "SAMPLE_SIZE",
                                            "TOTAL_SAMPLE_SIZE",
                                            "TOTALSAMPLESIZE"
                                          ),
                                          c.N_CAS = c(
                                            "N_CAS",
                                            "NCASE",
                                            "N_CASE",
                                            "N_CASES",
                                            "NCAS",
                                            "NCA",
                                            "NCASES",
                                            "CASES",
                                            "CASES_N",
                                            "FRQ_A" # +++CH: Is this correct here? +++JZ: Yes, this is the use of 'FRQ' in different standards. FRQ can refer to frequency of cases as well. I believe this interpretation is taken from Helena G's scripts. FRQ_A I believe is frequencey affected while FRQ_U is frequency unaffected. Not the same as FREQ_A and FREQ_U which refers to allele frequencies, maybe stratified per affected and unaffected.
                                          ),
                                          c.N_CON = c(
                                            "N_CON",
                                            "NCONTROL",
                                            "N_CONTROL",
                                            "N_CONTROLS",
                                            "NCON",
                                            "NCO",
                                            "N_CON",
                                            "NCONTROLS",
                                            "CONTROLS",
                                            "CONTROLS_N",
                                            "FRQ_U" # +++CH: Is this correct here? +++JZ: See above.
                                          ),
                                          c.NEF = c(
                                            "NEF",
                                            "NEFF",
                                            "NEFFECTIVE",
                                            "NE"
                                          ),
                                          #include FRQ_A? # +++CH: Who is this question by?
                                          c.FRQ = c(
                                            "FRQ", #+++JZ: Do we prefer FRQ being the standard column name for effect allele frequency?
                                            "FREQ",
                                            "MAF",
                                            "AF",
                                            "CEUAF",
                                            "FREQ1",
                                            "EAF",
                                            "FREQ1.HAPMAP",
                                            "FREQALLELE1HAPMAPCEU",
                                            "FREQ.ALLELE1.HAPMAPCEU",
                                            "EFFECT_ALLELE_FREQ",
                                            "FREQ.A1",
                                            "F_A",
                                            "F_U",
                                            "FREQ_A",
                                            "FREQ_U",
                                            "MA_FREQ",
                                            "MAF_NW",
                                            "FREQ_A1",
                                            "A1FREQ",
                                            "CODED_ALLELE_FREQUENCY",
                                            "FREQ_TESTED_ALLELE_IN_HRS",
                                            "EAF_HRC",
                                            "EAF_UKB"
                                          ),
                                          c.CHR = c(
                                            "CHR",
                                            "CH",
                                            "CHROMOSOME",
                                            "CHROM",
                                            "CHR_BUILD38",
                                            "CHR_BUILD37",
                                            "CHR_BUILD36",
                                            "CHR_B38",
                                            "CHR_B37",
                                            "CHR_B36",
                                            "CHR_ID",
                                            "SCAFFOLD",
                                            "HG19CHR",
                                            "CHR.HG19",
                                            "CHR_HG19",
                                            "HG18CHR",
                                            "CHR.HG18",
                                            "CHR_HG18",
                                            "CHR_BP_HG19B37",
                                            "HG19CHRC"
                                          ),
                                          c.BP = c(
                                            "ORIGBP",
                                            "BP",
                                            "POS",
                                            "POSITION",
                                            "LOCATION",
                                            "PHYSPOS",
                                            "GENPOS",
                                            "CHR_POSITION",
                                            "POS_B38",
                                            "POS_BUILD38",
                                            "POS_B37",
                                            "POS_BUILD37",
                                            "BP_HG19B37",
                                            "POS_B36",
                                            "POS_BUILD36",
                                            "POS.HG19",
                                            "POS.HG18",
                                            "POS_HG19",
                                            "POS_HG18",
                                            "BP_HG19",
                                            "BP_HG18",
                                            "BP.GRCH38",
                                            "BP.GRCH37",
                                            "POSITION(HG19)",
                                            "POSITION(HG18)",
                                            "POS(B38)",
                                            "POS(B37)"
                                          )
){
  
  column_names.upper <- toupper(column_names)
  column_names.orig <- column_names
  
  #+++CH: What does this bit do?
  #+++JZ: This part sets the names of all columns of each set to the first value of the set; the standard column name. Multiple columns which correspond to the same set may then end up with the same column name. Warnings are thrown for those below.
  column_names[column_names.upper %in% c.SNP] <- c.SNP[1]
  column_names[column_names.upper %in% c.A1] <- c.A1[1]
  column_names[column_names.upper %in% c.A2] <- c.A2[1]
  column_names[column_names.upper %in% c.BETA] <- c.BETA[1]
  column_names[column_names.upper %in% c.OR] <- c.OR[1] 
  column_names[column_names.upper %in% c.Z] <- c.Z[1] 
  column_names[column_names.upper %in% c.SE] <- c.SE[1]
  column_names[column_names.upper %in% c.INFO] <- c.INFO[1]
  column_names[column_names.upper %in% c.P] <- c.P[1]
  column_names[column_names.upper %in% c.N] <- c.N[1]
  column_names[column_names.upper %in% c.N_CAS] <- c.N_CAS[1]
  column_names[column_names.upper %in% c.N_CON] <- c.N_CON[1]
  column_names[column_names.upper %in% c.NEF] <- c.NEF[1]
  column_names[column_names.upper %in% c.FRQ] <- c.FRQ[1]
  column_names[column_names.upper %in% c.CHR] <- c.CHR[1]
  column_names[column_names.upper %in% c.BP] <- c.BP[1]
  
  if(stop_on_missing_essential_column){
    # Stop if any of these columns are not found
    if(!any(column_names=="SNP")) stop("\nCould not find the 'SNP' column.\n")
    if(!any(column_names=="A1")) stop("\nCould not find the 'A1' column.\n")
    if(!any(column_names=="A2")) stop("\nCould not find the 'A2' column.\n")
    if(!any(column_names=="BETA") & !any(column_names=="OR") & !any(column_names=="Z")) stop("Could not find any effect column.\n")  
  }
  
  if(!any(column_names=="P")) warning("\nCould not find the P-value column. Standard is 'P'.\n")
  if(!any(column_names=="FRQ")) warning("\nCould not find the 'FRQ' column.\n")
  
  # Warn if multiple of these columns are found #+++CH: Do we want to add more checks for duplicates?
  if(sum(column_names=="SNP")>1) warning("\nMultiple 'SNP' columns found!\n")
  if(sum(column_names=="P")>1) warning("\nMultiple 'P' columns found!\n")
  if(sum(column_names=="A1")>1) warning("\nMultiple 'A1' columns found!\n")
  if(sum(column_names=="A2")>1) warning("\nMultiple 'A2' columns found!\n")
  if(sum(column_names=="BETA")>1) warning("\nMultiple 'BETA' columns found!\n")
  if(sum(column_names=="OR")>1) warning("\nMultiple 'OR' columns found!\n")
  if(sum(column_names=="Z")>1) warning("\nMultiple 'Z' columns found!\n")
  if(sum(column_names=="SE")>1) warning("\nMultiple 'SE' columns found!\n")
  if(sum(column_names=="FRQ")>1) warning("\nMultiple 'FRQ' columns found!\n")
  if(sum(column_names=="INFO")>1) warning("\nMultiple 'INFO' columns found!\n")
  
  #+++JZ: Moved this part here to de-clutter the main script
  #rename (some) duplicate column names
  ##SNP
  iDup<-grep(pattern = "^SNP$",column_names)
  if(length(iDup)>1){
    iDup<-iDup[2:length(iDup)]
    column_names[iDup]<-"XSNP"
  }
  
  ##BP
  if('BP' %in% names(cSumstats)){
    iDup<-grep(pattern = "^BP$",column_names)
    if(length(iDup)>1){
      iDup<-iDup[2:length(iDup)]
      column_names[iDup]<-"XBP"
    }
  }
  
  ##FRQ
  if('FRQ' %in% names(cSumstats)){
    iDup<-grep(pattern = "^FRQ$",column_names)
    if(length(iDup)>1){
      iDup<-iDup[2:length(iDup)]
      column_names[iDup]<-"XFRQ"
    }
  }
  
  return(
    data.frame(
      std = column_names,
      orig = column_names.orig
    )
  )
}