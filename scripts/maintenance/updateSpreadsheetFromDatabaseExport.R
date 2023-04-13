#Update spreadsheet with data from the latest database export

library(data.table)
library(googlesheets4)

filepathDBExport <- normalizePath("data_gwas_sumstats_repository_2021/database_export_gwas_reference_category_phenotype.Rds",mustWork = T)
sheetLink <- "https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"

dbExport<-readRDS(file=filepathDBExport)
setDT(dbExport)
setkeyv(dbExport,cols = c("code"))

currentSheet <- as.data.frame(read_sheet(ss = sheetLink, col_names = T, range="SGDP_GWASLIST_EDITTHIS"))
currentSheet.cols<-colnames(currentSheet)
if(nrow(currentSheet)>0) currentSheet$x_row<-1:nrow(currentSheet)
setDT(currentSheet)
setkeyv(currentSheet,cols = c("code"))

resultsTable<-as.data.frame(matrix(data = NA,nrow = 0,ncol = 0))

for(iGWAS in 1:nrow(dbExport)){
  #iGWAS<-1
  cGWAS<-dbExport[iGWAS,]
  cat("\n",cGWAS$code[[1]])
  #check duplicate ID
  match<-currentSheet[code==eval(cGWAS$code[[1]]),]
  cNrow<-nrow(match)
  if(cNrow>0){
    if(cNrow>1) {
      cat("\tMulti")
      resultsTable[cGWAS$code[[1]],c("w.multi")]<-T
    } else {
      #check consistency of match
        if(nrow(currentSheet[
          code==eval(cGWAS$code[[1]])
          & phenotype_type==eval(cGWAS$phenotype_type[[1]])
          & sex==eval(cGWAS$sex[[1]])
          & ((is.na(pmid) & is.na(eval(cGWAS$pmid[[1]]))) | pmid==eval(cGWAS$pmid[[1]]))
          ,])==0){
            #inconsistent match - warn
          cat("\tInconsistent")
          resultsTable[cGWAS$code[[1]],c("w.inconsistent")]<-T
        }
    }
  } else {
    #add missing GWAS - use the function for this
    cat("\tAdded")
    resultsTable[cGWAS$code[[1]],c("a")]<-T
    tngpipeline::updateSpreadsheet(sheetLink = sheetLink, sumstats_meta = cGWAS)
  }
}

#inconsistent
print(rownames(resultsTable[which(resultsTable$w.inconsistent==T),]))

#qdded
print(rownames(resultsTable[which(resultsTable$a==T),]))

