# munge all repository sumstats

library(data.table)
library(googledrive)
library(googlesheets4)

#settings

folderpath.raw<-normalizePath("/scratch/prj/gwas_sumstats/original",mustWork = T)
folderpath.cleaned<-normalizePath("/scratch/prj/gwas_sumstats/cleaned",mustWork = T)
folderpath.munged2<-normalizePath("/scratch/prj/gwas_sumstats/munged_hc1kg",mustWork = T)

filepath.varlist<-normalizePath("/scratch/prj/gwas_sumstats/variant_lists/hc1kgp3.b38.mix.l2.jz2023.gz",mustWork = T)

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

for(iTrait in 1:nrow(currentSheet)){
  #iTrait<-1
  cTrait <- currentSheet[iTrait,]

  filepath.raw<-file.path(folderpath.raw,cTrait$file_name)
  filepath.cleaned<-file.path(folderpath.cleaned,paste0(cTrait$code,".gz"))


}
