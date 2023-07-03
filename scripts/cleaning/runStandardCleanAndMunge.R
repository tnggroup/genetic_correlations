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

column_options <- parse_args(column_parser)

cat("\n***Run standard clean and munge***\n")

# # # test
# column_options<-c()
# #projectFolderpath <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1",mustWork = T)
# projectFolderpath <- normalizePath("/scratch/prj/gwas_sumstats",mustWork = T)
# column_options$file <-paste0(
#   file.path(projectFolderpath,"cleaned","ACCU01.gz"),
#   ",",
#   file.path(projectFolderpath,"cleaned","ALCD03.gz")
# )
# #column_options$label <- "Binge Eating (Narrow),BMI 2018"
# column_options$code <- "ACCU01,ALCD03"
# #column_options$`sample-size`<-"830917,681275"
# #column_options$population<-"EUR"
# column_options$`reference-file-path`<- file.path(projectFolderpath,"variant_lists","hc1kgp3.b38.mix.l2.jz2023.gz")
# column_options$output<-'/scratch/groups/gwas_sumstats/munged_hc1kg'
#
#

#settings
filePaths <-  unlist(strsplit(column_options$file,split = ",",fixed = T))
traitCodes <- unlist(strsplit(column_options$code,split = ",",fixed = T))
#traitCodes <- c("BEN","BMI")
traitNames<-NULL
if(!is.null(column_options$label)) {
  traitNames <- unlist(strsplit(column_options$label,split = ",",fixed = T))
} else {
  traitNames <- traitCodes
  }
#traitNames <- c("Binge Eating (Narrow)","BMI 2018")
ancestrySetting<-NULL
if(!is.null(column_options$population)) ancestrySetting <- unlist(strsplit(column_options$population,split = ",",fixed = T))
referenceFilePath<-NULL
if(!is.null(column_options$`reference-file-path`)) referenceFilePath <- column_options$`reference-file-path`
N<-NULL
if(!is.null(column_options$`sample-size`)) N <- unlist(strsplit(column_options$`sample-size`,split = ",",fixed = T))
pathDirOutput<-NULL
if(!is.null(column_options$output)) pathDirOutput <- column_options$output
maf_filter<-NULL
if(!is.null(column_options$`filter.maf`)) maf_filter <- column_options$`filter.maf`
info_filter<-NULL
if(!is.null(column_options$`filter.info`)) info_filter <- column_options$`filter.info`

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

tngpipeline::standardPipelineCleanAndMunge(filePaths = filePaths, traitCodes = traitCodes, traitNames = traitNames,referenceFilePath = referenceFilePath, n_threads = n_threads, keep_indel = keep_indel, maf_filter = maf_filter,info_filter = info_filter, pathDirOutput = pathDirOutput)
