
#Johan Zvrskovec, 2021
# _Future_
# +++CH: These functions need a sentence explaining what they are doing. Who is going to do this?
#   +++CH: What is this comment about the number of SNPs? Can we move this outside of the chunk including explanation?
#   +++JKZ: It only checks the top 100,000 SNPs in a file. We don't know what happens if it detects less than 100,000

#Function to parse a SNP column as the typical RS-number that we expect it to contain
#Can deal with the BGENIE SNP column format
parse_SNP_column_as_rs_number <- function(text){
  #decide if BGENIE SNP format using top 100,000 SNPs
  #TODO this condition may be improved to not rely on the number of variants being >100,000
  #test
  #text<-files[[i]]$SNP
  if(sum(grepl(pattern = "^\\d+:\\w+_\\w+_\\w+", x= head(x = text, n=100000)))>90000){
    #extract and format rs-no
    indexesLengths<-regexec(pattern = "^\\d+:(\\w+)_\\w+_\\w+", text=text)
    matches<-regmatches(text,indexesLengths)
    return(lapply(X = matches, FUN = function(x)paste0("rs",x[2])))
  }
  
  text<-sub(pattern = "^chr",replacement = "",x = text, ignore.case = T)
  text<-sub(pattern = "^XY:",replacement = "25:",x = text, ignore.case = T)
  text<-sub(pattern = "^X:",replacement = "23:",x = text, ignore.case = T)
  text<-sub(pattern = "^Y:",replacement = "24:",x = text, ignore.case = T)
  text<-sub(pattern = "^MT:",replacement = "26:",x = text, ignore.case = T)
  text<-sub(pattern = "^M:",replacement = "26:",x = text, ignore.case = T)
  text<-sub(pattern = "^Un:",replacement = "0:",x = text, ignore.case = T)
  text<-sub(pattern = "_",replacement = ":",x = text, ignore.case = T)
  
  return(text)
}

#Function to parse a CHR column. it replaces chromosome character codes with numeric codes.
parse_CHR_column <- function(text){
  text<-trimws(text)
  text<-sub(pattern = "^chr",replacement = "",x = text, ignore.case = T)
  text<-sub(pattern = "^XY",replacement = "25",x = text, ignore.case = T)
  text<-sub(pattern = "^X",replacement = "23",x = text, ignore.case = T)
  text<-sub(pattern = "^Y",replacement = "24",x = text, ignore.case = T)
  text<-sub(pattern = "^MT",replacement = "26",x = text, ignore.case = T)
  text<-sub(pattern = "^M",replacement = "26",x = text, ignore.case = T)
  text<-sub(pattern = "^Un",replacement = "0",x = text, ignore.case = T)
  return(text)
}