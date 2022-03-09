#Function to calculate genetic correlations between a trait GWAS dataset and a set of other traits
#Johan Zvrskovec, 2022
#Based on Helena's previous script [sumstats_path]/scripts/runcorrelations_name.sh

one_vs_many_ldsc <- function(
  current_trait_code=NA_character_,
  current_trait_filepath,
  other_trait_codes=NA,
  other_trait_filepaths,
  ldsc_filepath="/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc/ldsc.py",
  ld_folderpath="/scratch/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/data/eur_w_ld_chr.1KG_Phase3",
  output_file_path=NA_character_

){
  
  if(length(current_trait_code<1)){
    current_trait_code <- basename(current_trait_filepath)[1]
  }
  
  if(length(other_trait_codes<1)){
    other_trait_codes <- basename(other_trait_filepaths)
  }
  
  if(is.na(output_file_path)) output_file_path <- paste0(current_trait_code,".ldsc")
  
  other.df <- data.frame(code=other_trait_codes,path=other_trait_filepaths)
  
  
  if(!substr(ld_folderpath,nchar(ld_folderpath),nchar(ld_folderpath))==base::.Platform$file.sep) ld_folderpath <- paste0(ld_folderpath,base::.Platform$file.sep)

  
	args <- c(
			paste0("--n-blocks 200"),
			paste0("--rg ", paste(c(current_trait_filepath,other_trait_filepaths),collapse = ",",sep = ",")),
			paste0("--ref-ld-chr ",ld_folderpath),
			paste0("--w-ld-chr ",ld_folderpath),
			paste0("--out ",output_file_path),
			paste0("--no-check-alleles")
	)
	print(paste(ldsc_filepath,paste(args, collapse=" ")))
  output <- system2(
        command=normalizePath(ldsc_filepath,mustWork=T),
        args=args,
        stdout = T
  		)
  
  toreturn <- c()
  toreturn$text.rows <- output
  
  itab = grep(pattern = "Summary of Genetic Correlation Results",x=output)+1
  
  dfstring <- paste(output[itab:(itab+length(other_trait_filepaths)+1)],collapse = "\n")
  
  toreturn$df <- read.table(text = dfstring, header = T)
  
  return(toreturn)
}

#test
# sumstats_path <- normalizePath("/scratch/groups/gwas_sumstats",mustWork = T)
# output <- one_vs_many_ldsc(
# 	current_trait_code="ANXI03",
# 	current_trait_filepath=file.path(sumstats_path,"munged","ANXI03.sumstats.gz"),
# 	other_trait_codes=c("DEPR05","DEPR08"),
# 	other_trait_filepaths=file.path(sumstats_path,"munged",c("DEPR05.sumstats.gz","DEPR08.sumstats.gz"))
# )
