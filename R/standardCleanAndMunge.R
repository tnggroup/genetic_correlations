#you will need to have these packages installed
# Standardised pipeline as an R function


#needs the shru package
#devtools::install_github("johanzvrskovec/shru")

standardPipelineCleaningAndMunging <- function(traitCodes,traitNames,referenceFilePath,n_threads=5,keep_indel=T,maf_filter=0.01,info_filter=0.6, outputFolderPath){

  sumstats_meta <- data.frame(
    name = traitNames
  )

  row.names(sumstats_meta) <- traitCodes
  sumstats_meta$code<-traitCodes

  varlist<-fread(file = referenceFilePath, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = n_threads, showProgress = F)

  #testTrait<-readFile(filePath = filePaths[iTrait])

  for(iTrait in 1:filePaths){
    #iTrait <-1 #test
    smungeResults <- shru::supermunge(filePaths = filePaths[iTrait], ref_df = varlist, traitNames = traitCodes[iTrait], ancestrySetting = ancestrySetting[iTrait], N = N[iTrait], keepIndel = keep_indel, writeOutput = F,filter.maf = maf_filter, filter.info = info_filter)
  }

}

standardPipelineExplicitProcessing <- function(traitCodes,traitNames,referenceFilePath,n_threads=5){

}




