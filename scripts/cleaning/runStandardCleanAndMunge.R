#New cleaning/munging routine to be run on CREATE
# Will run the standardised pipeline as a command line program

#install the pipeline package
#devtools::install_github("tnggroup/genetic_correlations",ref = 'dev_jz', auth_token = "YOUR_PAT")

library(tidyverse)
library(googlesheets4)
library(data.table)
library(optparse)

#command line options
column_parser <- OptionParser()

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
  opt_str = c("-r", "--reference-file-path"),
  type = "character",
  help = "Path to the reference variant list file."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-s", "--rs-synonyms-file-path"),
  type = "character",
  help = "Path to the RS synonyms file."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-o", "--output"),
  type = "character",
  help = "Path to where you want tp place the output files."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("--filter.maf"),
  type = "numeric",
  help = "MAF filter lower bound to include."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("--filter.info"),
  type = "numeric",
  help = "INFO filter lower bound to include."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("--process"),
  type = "character",
  help = "Standardised edits of file content when possible. True or false (as parseable by R). Default TRUE."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("--doPipelineSpecific"),
  type = "character",
  help = "Run pipeline specific parts (besides supermunge). True or false (as parseable by R). Default TRUE."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("--format"),
  type = "character",
  help = "Output format. Default supermunge default format. Other options are ldsc: ldsc compartible sumstat format, cojo: GCTA COJO sumstat format."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("--setntoneff"),
  type = "character",
  help = "True or false (as parseable by R), per dataset. Set N=NEFF in the output. This will be capped in case of backed out NEFF (formula version, rather than previously provided)."
)



column_options <- parse_args(column_parser)

cat("\n***Run standard clean and munge***\n")

# # # test
# column_options<-c()
# #projectFolderpath <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1",mustWork = T)
# projectFolderpath <- normalizePath("/scratch/prj/gwas_sumstats",mustWork = T)
# # column_options$file <-paste0(
# #   file.path(projectFolderpath,"cleaned","ACCU01.gz"),
# #   ",",
# #   file.path(projectFolderpath,"cleaned","ALCD03.gz")
# # )
# #column_options$label <- "Binge Eating (Narrow),BMI 2018"
# column_options$code <- "SMOK10"
# #column_options$code <- "ACCU01,ALCD03"
# #column_options$`sample-size`<-"830917,681275"
# #column_options$population<-"EUR"
# column_options$`reference-file-path`<- file.path(projectFolderpath,"variant_lists","hc1kgp3.b38.mix.l2.jz2023.gz")
# column_options$output<-'/scratch/groups/gwas_sumstats/munged_hc1kg'



#settings
filePaths <- NULL
if(!is.null(column_options$file)) filePaths <-  unlist(strsplit(column_options$file,split = ",",fixed = T))
traitCodes <- NULL
if(!is.null(column_options$code)) traitCodes <- unlist(strsplit(column_options$code,split = ",",fixed = T))
#traitCodes <- c("BEN","BMI")
traitNames<-NULL
if(!is.null(column_options$label)) {
  traitNames <- unlist(strsplit(column_options$label,split = ",",fixed = T))
} else if (length(traitCodes)>0){
  traitNames <- traitCodes
}

#traitNames <- c("Binge Eating (Narrow)","BMI 2018")
ancestrySetting<-NULL
if(!is.null(column_options$population)) ancestrySetting <- unlist(strsplit(column_options$population,split = ",",fixed = T))
referenceFilePath<-NULL
if(!is.null(column_options$`reference-file-path`)) referenceFilePath <- column_options$`reference-file-path`
arg.rsSynonymsFilePath<-NULL
if(!is.null(column_options$`rs-synonyms-file-path`)) arg.rsSynonymsFilePath <- column_options$`rs-synonyms-file-path`
N<-NULL
if(!is.null(column_options$`sample-size`)) N <- unlist(strsplit(column_options$`sample-size`,split = ",",fixed = T))
pathDirOutput<-NULL
if(!is.null(column_options$output)) pathDirOutput <- column_options$output
maf_filter<-NULL
if(!is.null(column_options$`filter.maf`)) maf_filter <- column_options$`filter.maf`
info_filter<-NULL
if(!is.null(column_options$`filter.info`)) info_filter <- column_options$`filter.info`

arg.process <- TRUE
if(!is.null(column_options$process)) arg.process <-  as.logical(column_options$process)

arg.doPipelineSpecific <- TRUE
if(!is.null(column_options$doPipelineSpecific)) arg.doPipelineSpecific <-  as.logical(column_options$doPipelineSpecific)

arg.outputFormat <- NULL
if(!is.null(column_options$format)) arg.outputFormat <-  column_options$format

arg.setNToNEFF <- NULL
if(!is.null(column_options$setntoneff)) arg.setNToNEFF <-  as.logical(unlist(strsplit(column_options$setntoneff,split = ",",fixed = T)))

#hard coded options
n_threads <- 5
keep_indel <- TRUE
#maf_filter <- 0.01
#info_filter <- 0.6

#groupFolderPath <- normalizePath("/scratch/prj/gwas_sumstats",mustWork = T)
#groupFolderPath <- normalizePath("/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/gwas_sumstats_test",mustWork = T) #for test


cat("\nFile paths:")
print(filePaths)

cat("\nTrait names:")
print(traitNames)

cat("\nTrait codes:")
print(traitCodes)

cat("\nOutput path:")
print(pathDirOutput)

tngpipeline::standardPipelineCleanAndMunge(filePaths = filePaths, traitCodes = traitCodes, traitNames = traitNames,referenceFilePath = referenceFilePath, rsSynonymsFilePath=arg.rsSynonymsFilePath, n_threads = n_threads, keep_indel = keep_indel, maf_filter = maf_filter,info_filter = info_filter, process=arg.process, doPipelineSpecific=arg.doPipelineSpecific, pathDirOutput = pathDirOutput, setNtoNEFF=arg.setNToNEFF, outputFormat = arg.outputFormat)
