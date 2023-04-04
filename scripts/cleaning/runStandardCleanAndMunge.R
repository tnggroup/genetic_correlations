#New cleaning/munging routine to be run on CREATE
# Will run the standardised pipeline as a command line program

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
  opt_str = c("-ref", "--reference-file-path"),
  type = "character",
  help = "Path to the reference variant list file."
)

column_parser <- add_option(
  object = column_parser,
  opt_str = c("-o", "--output"),
  type = "character",
  help = "Path to gwas_sumstats group folder."
)


# # test
# projectFolderpath <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1",mustWork = T)
# column_options$file <-paste0(
#   file.path(projectFolderpath,"data","gwas_sumstats","raw","PGC3_ED_2022","daner_BENARROW.gz"),
#   ",",
#   file.path(projectFolderpath,"data","gwas_sumstats","raw","Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
# )
# column_options$label <- "Binge Eating (Narrow),BMI 2018"
# column_options$code <- "BEN,BMI"
# column_options$`sample-size`<-"830917,681275"
# column_options$population<-"EUR"
# column_options$`reference-file-path`<- file.path(projectFolderpath,"data","variant_lists","combined.hm3_1kg.snplist.vanilla.jz2020.gz")
# column_options$output<-'/scratch/groups/gwas_sumstats'

column_options <- parse_args(column_parser)

#settings
filePaths <-  unlist(strsplit(column_options$file,split = ",",fixed = T))
traitNames <- unlist(strsplit(column_options$label,split = ",",fixed = T)) #not working now for some reason
#traitNames <- c("Binge Eating (Narrow)","BMI 2018")
traitCodes <- unlist(strsplit(column_options$code,split = ",",fixed = T)) #not working now for some reason
#traitCodes <- c("BEN","BMI")
ancestrySetting <- unlist(strsplit(column_options$population,split = ",",fixed = T))
referenceFilePath <- column_options$`reference-file-path`
N <- unlist(strsplit(column_options$`sample-size`,split = ",",fixed = T))




#hard coded options
n_threads <- 5
keep_indel <- TRUE
maf_filter <- 0.01
info_filter <- 0.6

groupFolderPath <- normalizePath("/scratch/prj/gwas_sumstats",mustWork = T)




tngpipeline::standardPipelineCleaningAndMunging(traitCodes = traitCodes, traitNames = traitNames,)
