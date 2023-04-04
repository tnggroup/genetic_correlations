#Regression test function for the genetic correlations pipeline
#Johan Zvrskovec, 2022

#Put code here that will test your part of the main script. Running this script should verify that nothing breaks when implementing further changes into the main script. A test can be to run the main script with a defined input and check the results. This would work much easier if the main script was a function, or had sub-parts in functions.

regression_test <- function(){
  
  #test column name standardisation
  source(file = "scripts/functions/standardise_GWAS_column_names.R")
  testColumnNames<-c("SNP","A1","A2","MAF","BETA","SE","N","P") #a standard set of column names
  t<-standardise_GWAS_column_names(column_names = testColumnNames) #run with all defaults
  if(!any(t$std=="SNP") & any(t$std=="A1") & any(t$std=="A2") & any(t$std=="FRQ") & any(t$std=="BETA") & any(t$std=="SE") & any(t$std=="N") & any(t$std=="P")) stop("Column name standardisation is not working properly!")
  
  
  return(T)
}

#example, commented out - runs the whole test.
#regression_test()
