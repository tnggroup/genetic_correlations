library(shru)
library(optparse)

clParser <- OptionParser()
clParser <- add_option(clParser, c("-t", "--task"), type="character", default="0",
                       help="Task - used to specify trait code [default %default]")

clParser <- add_option(clParser, c("-a", "--task_argument"), type="character", default=NA,
                       help="Argument - used for specifying the level string [default %default]")

clOptions <- parse_args(clParser)

cTraitString <- clOptions$task
cLevelString <- clOptions$task_argument

#settings
filepathSNPReference <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/data/combined.hm3_1kg.snplist.vanilla.jz2020.txt", mustWork = T)
folderpathLDscores <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/data/eur_w_ld_chr.1KG_Phase3", mustWork = T)
folderpathEvaluationSumstats <- normalizePath("/users/k1806347/brc_scratch/Analyses/sumstat_imputation/filtered", mustWork = T)
folderpathEvaluationOutput <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/working_directory", mustWork = T)
folderpathEvaluationOutputLIMP <- file.path(folderpathEvaluationOutput,"LIMP")
dir.create(path = folderpathEvaluationOutputLIMP)

#traitList <- c("ADHD05","ANXI02","COAD01","OBES01","SMOK04")
#levelList <- c("05","1","15","2","25","3")

#for(cTraitString in traitList){
  #cTraitString <- "ADHD05"
  cFolderpathEvaluationSumstats <- file.path(folderpathEvaluationSumstats,cTraitString)
  cFileList <- file.path(cFolderpathEvaluationSumstats,paste0(cTraitString,".cleaned.",cLevelString,"_missing.gz"))
  cNameList <- paste0(cTraitString,".",cLevelString,".imputed")
  cOutputDir <- file.path(folderpathEvaluationOutputLIMP,cTraitString)
  dir.create(path = cOutputDir)
  supermunge(
    filePaths = cFileList,
    refFilePath = filepathSNPReference,
    ldDirPath=folderpathLDscores,
    traitNames = cNameList,
    imputeFromLD=T,
    pathDirOutput = cOutputDir
  ) 
#}