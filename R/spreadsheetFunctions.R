authenticateSpreadsheet <- function(
    serviceAccountTokenPath="/scratch/prj/gwas_sumstats/tngpipeline/tngpipeline-8130dbd7d58a.json",
    serviceAccountEmail="johan.kallberg_zvrskovec@kcl.ac.uk"){
  #configure gs4 for non-browser login #https://gargle.r-lib.org/articles/non-interactive-auth.html
  drive_auth(email = serviceAccountEmail, path = serviceAccountTokenPath)
  #drive_auth(email = "jane_doe@example.com") # gets a suitably scoped token
  # and stashes for googledrive use
  gs4_auth(token = drive_token())            # registers token with googlesheets4
}

readSpreadsheet <- function(
    sheetLink = "https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"
    ){
  currentSheet <- as.data.frame(read_sheet(ss = sheetLink, col_names = T, range="SGDP_GWASLIST_EDITTHIS"))
  currentSheet.cols<-colnames(currentSheet)
  if(nrow(currentSheet)>0) currentSheet$x_row<-1:nrow(currentSheet)
  setDT(currentSheet)
  setkeyv(currentSheet,cols = c("code"))
  return(currentSheet)
}


assignCodeFromSortAndSpreadsheet <- function(
    sort,
    sheetLink="https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"
){

  #sheetLink <- "https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"

  currentSheet <- as.data.frame(read_sheet(ss = sheetLink, col_names = T, range="SGDP_GWASLIST_EDITTHIS"))
  sorts<-currentSheet[grep(pattern = paste0("^",toupper(sort)), x = currentSheet$code),]
  if(nrow(sorts)>0){
    indexesLengths<-regexec(pattern = "\\d+", text=sorts$code)
    matches<-regmatches(sorts$code,indexesLengths)
    nums<-as.integer(unique(unlist(matches)))
    return(paste0(sort,shru::padStringLeft(s = paste0("",(max(nums)+1)), padding = "0", targetLength = 2)))
  } else return(paste0(sort,"01"))

}



updateSpreadsheet <- function(sheetLink,sumstats_meta){
  #sheetLink <- "https://docs.google.com/spreadsheets/d/1gjKI0OmYUxK66-HoXY9gG4d_OjiPJ58t7cl-OsLK8vU/edit?usp=sharing"
  #sumstats_meta <- cGWAS

  currentSheet <- as.data.frame(read_sheet(ss = sheetLink, col_names = T, range="SGDP_GWASLIST_EDITTHIS"))
  cols<-colnames(currentSheet)

  if(nrow(currentSheet)>0) currentSheet$x_row<-1:nrow(currentSheet)
  setDT(currentSheet)
  setkeyv(currentSheet,cols = c("code"))

  # sumstats_meta_dummy <- data.frame(
  #   code=c("ANXI04","DEPR05","NEUR02"),
  #   name=c("Generalised anxiety symptoms","MDD, narrow","Neuroticism"),
  #   ancestry=c("EUR","EUR","EUR"),
  #   n_ca=c(19012,58113,449484),
  #   n_co=c(58113,25632,NA_integer_)
  #   )
  # sumstats_meta_dummy$n_tot<-sumstats_meta_dummy$n_ca + ifelse(is.na(sumstats_meta_dummy$n_co),0,sumstats_meta_dummy$n_co)

  #dfToInsert<-as.data.frame(matrix(data=NA,nrow = 0, ncol = length(colnames(currentSheet))))
  #colnames(dfToInsert)<-colnames(currentSheet)
  dfToInsert <- currentSheet
  dfToInsert<-dfToInsert[FALSE,]

  setDT(sumstats_meta)
  #add missing columns
  if(!any(colnames(sumstats_meta)=="doi")) sumstats_meta[,doi:=NA_character_]
  if(!any(colnames(sumstats_meta)=="permissions")) sumstats_meta[,permissions:=NA_character_]
  sumstats_meta$phenotype_type<-as.character(sumstats_meta$phenotype_type)
  sumstats_meta$assembly<-as.character(sumstats_meta$assembly)
  sumstats_meta$ancestry<-as.character(sumstats_meta$ancestry)
  sumstats_meta$sex<-as.character(sumstats_meta$sex)
  sumstats_meta$permissions<-as.character(sumstats_meta$permissions)
  sumstats_meta$dependent_variable<-as.character(sumstats_meta$dependent_variable)
  sumstats_meta<-sumstats_meta[,.(
    name,
    code,
    n_cases,
    n_controls,
    sample_size_discovery=n_total,
    download_link,
    filename,
    platform,
    n_details,
    ancestry_details,
    permissions,
    assembly,
    uk_biobank,
    ancestry,
    doi,
    pmid,
    trait_detail=phenotype,
    phenotype_type,
    ancestry,
    sex,
    phenotype,
    phenotype_category=category_name,
    notes,
    dependent_variable,
    year,
    consortium
  )]
  setkeyv(sumstats_meta,cols = c("code"))

  print(colnames(sumstats_meta))
  dfToInsert<-rbindlist(list(
    dfToInsert,
    sumstats_meta
  ), fill = T)


  dfToInsert <- dfToInsert[,..cols]
  for(iRec in 1:nrow(dfToInsert)){
    #iRec<-1
    cRec <- dfToInsert[iRec,]
    if(!is.na(cRec$code)){
      if(any(currentSheet$code==cRec$code)){
        #inactivated update for now
        # toUpdate<-currentSheet[currentSheet$code==cRec$code,]$x_row
        # range_write(ss = sheetLink,data = cRec,cell_rows(toUpdate[[1]]:toUpdate[[1]]))
      } else {
        sheet_append(ss = sheetLink,data = cRec, sheet="SGDP_GWASLIST_EDITTHIS") #this does not recognise some columns
      }
    }
  }
}
