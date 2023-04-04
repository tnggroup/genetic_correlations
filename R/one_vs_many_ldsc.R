#Function to calculate genetic correlations between a trait GWAS dataset and a set of other traits. The first specified trait will be used as the main trait to be used for computing genetic correlations against the rest of the traits, as per original ldsc.
#Johan Zvrskovec, 2022
#Based on Helena's previous script [sumstats_path]/scripts/runcorrelations_name.sh

one_vs_many_ldsc <- function(
  trait_codes=NA, #A list of codes to use for your traits. If not specified, the codes will be set to the filenames of your gwas sumstat files as specified below, not including file sufix.
  trait_filepaths, #A list of paths to your trait gwas sumstat files to compute genetic correlations against.
  ldsc_filepath=NA_character_, #Path to ldsc to use.
  ld_folderpath=NA_character_, #Path to a folder containing LD-scores as used by ldsc.
  output_file_path=NA_character_ #Path to the log file where you want the mandatory ldsc log output to be saved. Defaults to the main trait code and the file suffix 'ldsc' in the current folder.
){
  
  #Set the trait codes if not specified in the argument.
  if(length(trait_codes<1)){
    trait_codes <- basename(trait_filepaths)
  }
  
  #Errors when not providing paths to ldsc or to ld-score folder
  if(is.na(ldsc_filepath)) stop("You need to set the path to your installation of ldsc through the 'ldsc_filepath' argument!")
  if(is.na(ld_folderpath)) stop("You need to set the path to your ld-scores (folder) through the 'ld_folderpath' argument!")
  
  #Set the output file path if not specified.
  if(is.na(output_file_path)) output_file_path <- paste0(current_trait_code,".ldsc")
  
  #Append the mandatory folder (file?) separators to the end of the ld folder path if they are missing.
  if(!substr(ld_folderpath,nchar(ld_folderpath),nchar(ld_folderpath))==base::.Platform$file.sep) ld_folderpath <- paste0(ld_folderpath,base::.Platform$file.sep)

  #Set up arguments for ldsc in a list.
	args <- c(
			paste0("--n-blocks 200"),
			paste0("--rg ", paste(trait_filepaths,collapse = ",",sep = ",")),
			paste0("--ref-ld-chr ",ld_folderpath),
			paste0("--w-ld-chr ",ld_folderpath),
			paste0("--out ",output_file_path),
			paste0("--no-check-alleles")
	)
	
	#Print the current arguments to the screen, as a reminder
	print(paste(ldsc_filepath,paste(args, collapse=" ")))
	
	#Run ldsc with the configured arguments and printing the output to the standard out stream (the screen).
  output <- system2(
        command=normalizePath(ldsc_filepath,mustWork=T),
        args=args,
        stdout = T
  		)
  
  #Set up the output object that will be returned. Will be composed of multiple objects in a list.
  toreturn <- c()
  toreturn$text.rows <- output #Store the output from ldsc
  
  #Get the index from where to parse the ldsc results from the ldsc text log output
  itab = grep(pattern = "Summary of Genetic Correlation Results",x=output)+1
  
  #Extract the text from the ldsc log outoput that is to be parsed to a dataframe
  dfstring <- paste(output[itab:(itab+length(other_trait_filepaths)+1)],collapse = "\n")
  
  #Convert the text to a dataframe
  toreturn$df <- read.table(text = dfstring, header = T)
  
  return(toreturn)
}

#test - this does not work - you need to add your paths to ldsc and ld scores
# sumstats_path <- normalizePath("/scratch/groups/gwas_sumstats",mustWork = T)
# output <- one_vs_many_ldsc(
# 	trait_codes=c("ANXI03","DEPR05","DEPR08"),
# 	trait_filepaths=file.path(sumstats_path,"munged",c("ANXI03.sumstats.gz","DEPR05.sumstats.gz","DEPR08.sumstats.gz"))
# )
# output$df #access the ldsc data in a dataframe
