# munge all repository sumstats
# devtools::install_github("tnggroup/genetic_correlations",ref = 'dev_jz', auth_token = "YOUR_PAT")

library(data.table)
library(googledrive)
library(googlesheets4)

#settings

folderpath.raw<-normalizePath("/scratch/prj/gwas_sumstats/original",mustWork = T)
folderpath.cleaned<-normalizePath("/scratch/prj/gwas_sumstats/cleaned",mustWork = T)
folderpath.munged2<-normalizePath("/scratch/prj/gwas_sumstats/munged_hc1kg",mustWork = T)

filepath.varlist<-normalizePath("/scratch/prj/gwas_sumstats/variant_lists/hc1kgp3.b38.mix.l2.jz2023.gz",mustWork = T)

filepath.commands.out<-"mungeAll.txt"

sheetLink <- "https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"
serviceAccountTokenPath=normalizePath("/scratch/prj/gwas_sumstats/tngpipeline/tngpipeline-8130dbd7d58a.json",mustWork = T)


#setup gsheets
drive_auth(email = "johan.kallberg_zvrskovec@kcl.ac.uk", path = serviceAccountTokenPath)
gs4_auth(token = drive_token())

#read sheet
currentSheet <- as.data.frame(read_sheet(ss = sheetLink, col_names = T, range="SGDP_GWASLIST_EDITTHIS"))
currentSheet.cols<-colnames(currentSheet)
if(nrow(currentSheet)>0) currentSheet$x_row<-1:nrow(currentSheet)
setDT(currentSheet)
setkeyv(currentSheet,cols = c("code"))

if(file.exists(filepath.commands.out)){
    file.remove(filepath.commands.out)
    file.create(filepath.commands.out)
}

for(iTrait in 1:nrow(currentSheet)){
  #iTrait<-1
  cTrait <- currentSheet[iTrait,]

  filepath.raw<-file.path(folderpath.raw,cTrait$file_name)
  filepath.cleaned<-file.path(folderpath.cleaned,paste0(cTrait$code,".gz"))
  if(file.exists(filepath.raw)) {
    filepath.touse <-filepath.raw
  } else {
    filepath.touse <-filepath.cleaned
  }

  wrapCommand<-paste0("sleep 30; Rscript ~/project/genetic_correlations/scripts/cleaning/runStandardCleanAndMunge.R -f '",filepath.touse,"' -c '",cTrait$code,"' -r '",filepath.varlist,"' -p '",cTrait$ancestry,"' --filter.maf 0.01 --filter.info 0.6;")

  args <- c(
    paste0("--time 01:00:00"),
    paste0("--partition cpu"),
    paste0("--job-name=\"tng",iTrait,"\""),
    paste0("--ntasks 1"),
    paste0("--cpus-per-task 5"),
    paste0("--mem 64G"),
    paste0("--wrap \"",wrapCommand,"\""),
    paste0("--output \"",cTrait$code,".$(date +%Y%m%d).out.txt\"") #this probably does not work
  )

  # output <- system2(
  #   command="sbatch",
  #   args=args,
  #   stdout = T
  # )

  fullTextCommand<-paste("sbatch",paste(args, collapse=" "),";")

  data.table::fwrite(x = list(fullTextCommand), append = T, file = filepath.commands.out, quote = F)

}
