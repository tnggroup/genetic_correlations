# Create IUPAC code function
snp_iupac <- function(x = NA, y = NA){
  if(length(x) != length(y)){
    print('x and y are different lengths')
  } else {
    iupac <- rep(NA, length(x))
    iupac[x == 'A' & y =='T' | x == 'T' & y =='A'] <- 'W'
    iupac[x == 'C' & y =='G' | x == 'G' & y =='C'] <- 'S'
    iupac[x == 'A' & y =='G' | x == 'G' & y =='A'] <- 'R'
    iupac[x == 'C' & y =='T' | x == 'T' & y =='C'] <- 'Y'
    iupac[x == 'G' & y =='T' | x == 'T' & y =='G'] <- 'K'
    iupac[x == 'A' & y =='C' | x == 'C' & y =='A'] <- 'M'
    return(iupac)
  }
}

# Create function to change allele to complement
snp_allele_comp<-function(x=NA){
  x_new<-x
  x_new[x == 'A']<-'T'
  x_new[x == 'T']<-'A'
  x_new[x == 'G']<-'C'
  x_new[x == 'C']<-'G'
  x_new[!(x %in% c('A','T','G','C'))]<-NA
  return(x_new)
}
