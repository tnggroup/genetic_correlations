#function to parse text ancestry values from the google sheet

parseAncestryText <- function(l.ancestryText){
  sapply(l.ancestryText, FUN = function(ancestryText){
    if(is.na(ancestryText)) return("MIX")
    #ancestryText.orig<-ancestryText
    ancestryText<-toupper(as.character(ancestryText))
    if(any(ancestryText == c("MIX","TRANSETHNIC","WORLDWIDE","WORLD"))) ancestryText<-"MIX"
    if(any(ancestryText == c("EUR","EUROPEAN","EUROPE","MOSTLY EUROPEAN"))) ancestryText<-"EUR"
    if(any(ancestryText == c("AFR","AFRICAN","AFRICA"))) ancestryText<-"AFR"
    if(any(ancestryText == c("AMR","AD MIXED AMERICAN"))) ancestryText<-"AMR"
    if(any(ancestryText == c("EAS","EAST ASIAN","EAST ASIA"))) ancestryText<-"EAS"
    if(any(ancestryText == c("SAS","SOUTH ASIAN","SOUTH ASIA"))) ancestryText<-"SAS"

    #fallback
    if(!any(ancestryText==c("MIX","EUR","AFR","AMR","EAS","SAS"))) {
      warning("No ancestry could be matched to the provided string")
      ancestryText<-"MIX"
    }

    return(ancestryText)
  }
           )
}
