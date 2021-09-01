#Based on the amazing work by Grotzinger, A. D. et al. Nat. Hum. Behav. 3, 513–525 (2019) and Bulik-Sullivan, B. K. et al. Nat. Genet. 47, 291–295 (2015).
#This is a modified copy of the multivariate ldsc function in GenomicSEM: https://github.com/GenomicSEM/GenomicSEM
#Forked and modified by Johan Zvrskovec at: https://github.com/johanzvrskovec/GenomicSEM

ldsc.mod <- function(traits, sample.prev, population.prev, ld, wld,
                trait.names = NULL, sep_weights = FALSE, chr = 22,
                n.blocks = 600, ldsc.log = NULL, stand = FALSE,select=FALSE,chisq.max = NA, info.filter = .6,maf.filter=0.01, N=NULL, forceN=FALSE, limited=FALSE) {
  
  # traits = project$sumstats.sel$mungedpath[c(15,16)]
  # sample.prev =  project$sumstats.sel$samplePrevalence[c(15,16)]
  # population.prev = project$sumstats.sel$populationPrevalence[c(15,16)]
  # trait.names = project$sumstats.sel$code[c(15,16)]
  # ld = project$folderpath.data.mvLDSC.ld
  # wld = project$folderpath.data.mvLDSC.ld
  # n.blocks = 600
  # info.filter = 0.6
  # maf.filter = 0.01
  # ldsc.log = project$setup.code.date
  # stand = TRUE
  # select=FALSE
  # chisq.max = NA
  # chr = 22
  # sep_weights = FALSE
  # N = project$sumstats.sel$n_total[c(15,16)]
  # forceN=FALSE
  
  
  LOG <- function(..., print = TRUE) {
    msg <- paste0(...)
    if (print) print(msg)
    cat(msg, file = log.file, sep = "\n", append = TRUE)
  }
  
  time <- proc.time()
  
  begin.time <- Sys.time()
  
  if(is.null(ldsc.log)){
    logtraits<-gsub(".*/","",traits)
    log2<-paste(logtraits,collapse="_")
    if(object.size(log2) > 200){
      log2<-substr(log2,1,100)
    }
    log.file <- file(paste0(log2, "_ldsc.log"),open="wt")
  }else{log.file<-file(paste0(ldsc.log, "_ldsc.log"),open="wt")}
  
  LOG("Multivariate ld-score regression of ", length(traits), " traits ",
      "(", paste(traits, collapse = " "), ")", " began at: ", begin.time)
  
  
  if(select == "ODD" | select == "EVEN"){
  odd<-seq(1,chr,2)
  even<-seq(2,chr,2)
  }
  
  # Dimensions
  n.traits <- length(traits)
  n.V <- n.traits * (n.traits + 1) / 2
  
  if(!(is.null(trait.names))){
    check_names<-str_detect(trait.names, "-")
    if(any(check_names))
      warning("Your trait names specified include mathematical arguments (e.g., + or -) that will be misread by lavaan. Please rename the traits using the trait.names argument.")
  }
  
  if(length(traits)==1)
    warning("Our version of ldsc requires 2 or more traits. Please include an additional trait.")
  
  
  # Storage:
  cov <- matrix(NA,nrow=n.traits,ncol=n.traits)
  V.hold <- matrix(NA,nrow=n.blocks,ncol=n.V)
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  Liab.S <- rep(1, n.traits)
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)
  
  
  #########  READ LD SCORES:
  LOG("Reading in LD scores")
 
  if(select == FALSE){
  x <- do.call("rbind", lapply(1:chr, function(i) {
    suppressMessages(read_delim(
      file.path(ld, paste0(i, ".l2.ldscore.gz")),
      delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
  }))
  }
  
  if(select == "ODD"){
    x <- do.call("rbind", lapply(odd, function(i) {
      suppressMessages(read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
  }
  
  if(select == "EVEN"){
    x <- do.call("rbind", lapply(even, function(i) {
      suppressMessages(read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
  }

  if(is.numeric(select)){
    x <- do.call("rbind", lapply(select, function(i) {
      suppressMessages(read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
  }
  
  #mod addition, correct SNP to lower case in case it is upper case
  x$SNP<-tolower(x$SNP)
  x$CM <- NULL
  x$MAF <- NULL
  
  
  ######### READ weights:
  if(sep_weights){
    if(select == FALSE){
    w <- do.call("rbind", lapply(1:chr, function(i) {
      suppressMessages(read_delim(
        file.path(wld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
    }
    if(select == "EVEN"){
      w <- do.call("rbind", lapply(even, function(i) {
        suppressMessages(read_delim(
          file.path(wld, paste0(i, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      }))
    }
    if(select == "ODD"){
      w <- do.call("rbind", lapply(even, function(i) {
        suppressMessages(read_delim(
          file.path(wld, paste0(i, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      }))
    }
      if(is.numeric(select)){
        w <- do.call("rbind", lapply(select, function(i) {
          suppressMessages(read_delim(
            file.path(wld, paste0(i, ".l2.ldscore.gz")),
            delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
        }))
    }
  }else{w<-x}
  
  #mod addition, correct SNP to lower case in case it is upper case
  w$SNP<-tolower(w$SNP)
  w$CM <- NULL
  w$MAF <- NULL
  
  colnames(w)[ncol(w)] <- "wLD"
  
  ### READ M
  
  if(select == FALSE){
  m <- do.call("rbind", lapply(1:chr, function(i) {
    suppressMessages(read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
  }))
  }
  
  if(select == "EVEN"){
    m <- do.call("rbind", lapply(even, function(i) {
      suppressMessages(read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
  }
  
  if(select == "ODD"){
    m <- do.call("rbind", lapply(odd, function(i) {
      suppressMessages(read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
  }
  
  if(is.numeric(select)){
    m <- do.call("rbind", lapply(select, function(i) {
      suppressMessages(read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
  }
  
  M.tot <- sum(m)
  m <- M.tot
  
  ### READ ALL CHI2 + MERGE WITH LDSC FILES
  s <- 0
  
  all_y <- lapply(traits, function(chi1) {
    #chi1<-traits[1]
    ## READ chi2
    ## mod change - use read.table
    y1 <- suppressMessages(read.table(
      chi1, header=T, quote="\"",fill=T,na.string=c(".",NA,"NA","")))
    #mod addition, harmonise character case
    y1$SNP<-tolower(as.character(y1$SNP))
    y1$A1<-toupper(as.character(y1$A1))
    y1$Z<-as.numeric(y1$Z)
    #mod addition - use FRQ - frequency of A1 rather than MAF. use MAF as FRQ in case there is no FRQ.
    if('FRQ' %in% names(y1)) y1$FRQ<-as.numeric(y1$FRQ)
    if('MAF' %in% names(y1)) y1$MAF<-as.numeric(y1$MAF)
    if('INFO' %in% names(y1)) y1$INFO<-as.numeric(y1$INFO)
    #mod addition - use explicit new N value
    if(!is.null(N) & length(N)>=length(traits)) {
      if(forceN) { 
        y1$N<-as.integer(N[s])
        LOG("Added explicit N=",N[s])
      } else {
        cond<-is.na(y1$N)
        y1$N<-ifelse(cond,N[s],y1$N)
        if(sum(cond)>0) LOG("Added explicit N=",N[s], " to ",sum(cond)," SNPs where it was missing.")
        }
    }
    
    LOG("Read in summary statistics [", s <<- s + 1, "/", n.traits, "] from: ", chi1)
    
    ##mod addition, check existence of FRQ, otherwise, use MAF
    if(!any(colnames(y1)=="FRQ") & any(colnames(y1)=="MAF")){
      y1$FRQ<-y1$MAF
    }
    
    ##mod addition
    ## REMOVE SNPs MAF<maf.filter and INFO<info.filter
    if(!is.na(maf.filter)){
      if("FRQ" %in% names(y1)){
        rm <- (!is.na(y1$FRQ) & (y1$FRQ<maf.filter | (1-y1$FRQ)<maf.filter))
        y1 <- y1[!rm, ]
        LOG("Removing ", sum(rm), " SNPs with MAF <", maf.filter, "; ", nrow(y1), " remain")
      } else {
        #warning("The dataset does not contain a MAF column to apply the specified filter on.")
        LOG("Warning: The dataset does not contain a FRQ or MAF column to apply the specified filter on.")
      }
    }
    if(!is.na(info.filter)){
      if("INFO" %in% names(y1)){
        rm <- (!is.na(y1$INFO) & y1$INFO<info.filter)
        y1 <- y1[!rm, ]
        LOG("Removing ", sum(rm), " SNPs with INFO <", info.filter, "; ", nrow(y1), " remain")
      } else {
        #warning("The dataset does not contain an INFO column to apply the specified filter on.")
        LOG("Warning: The dataset does not contain an INFO column to apply the specified filter on.")
      }
    }
    
    ## Merge files
    y1.columns.orig<-c("SNP", "N", "Z", "A1")
    y1.columns<-y1.columns.orig
    ##mod addition
    if('FRQ' %in% names(y1)) y1.columns<-c(y1.columns.orig,'FRQ')
    if('MAF' %in% names(y1)) y1.columns<-c(y1.columns.orig,'MAF')
    if('INFO' %in% names(y1)) y1.columns<-c(y1.columns.orig,'INFO')
    y1<-as_tibble(y1)
    
    merged <- merge(y1[, y1.columns.orig], w[, c("SNP", "wLD")], by = "SNP", sort = FALSE)
    merged <- merge(merged, x, by = "SNP", sort = FALSE)
    merged <- merged[with(merged, order(CHR, BP)), ]
    
    LOG("Out of ", nrow(y1), " SNPs, ", nrow(merged), " remain after merging with LD-score files")
    LOG("Columns after merge: ", paste(colnames(merged),collapse = " ")) ##mod addition
    if(all(is.na(merged$N))) LOG("Warning: The data has no N values!!") ##mod addition
    
    
    ## REMOVE SNPS with excess chi-square:
    
    if(is.na(chisq.max)){
    chisq.max <- max(0.001 * max(merged$N), 80)
    }
    rm <- (merged$Z^2 > chisq.max)
    merged <- merged[!rm, ]
    
    LOG("Removing ", sum(rm), " SNPs with Chi^2 > ", chisq.max, "; ", nrow(merged), " remain")
    
    ##mod addition - remove any rows with na values as was originally performed earlier before merge
    if(any(is.na(merged))){
      nBefore <- nrow(merged)
      merged <- na.omit(merged)
      LOG("Removing ", nBefore-nrow(merged), " SNPs with NA-values present; ", nrow(merged), " remain")
    }
    
    merged
  })
  
  # count the total nummer of runs, both loops
  s <- 1
  
  for(j in 1:n.traits){
    #j<-1
    chi1 <- traits[j]
    
    y1 <- all_y[[j]]
    y1$chi1 <- y1$Z^2
    
    for(k in j:length(traits)){
      #k<-1
      ##### HERITABILITY code
      
      if(j == k){
        
        LOG("     ", "     ", print = FALSE)
        LOG("Estimating heritability [", s, "/", n.V, "] for: ", chi1)
        
        samp.prev <- sample.prev[j]
        pop.prev <- population.prev[j]
        
        merged <- y1
        n.snps <- nrow(merged)
        
        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        
        
        #### MAKE WEIGHTS:
        
        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- pmax(merged$L2, 1)
        merged$w.ld <- pmax(merged$wLD, 1)
        merged$c <- tot.agg*merged$N/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        merged$weights <- merged$initial.w/sum(merged$initial.w)
        
        N.bar <- mean(merged$N)
        
        
        ## preweight LD and chi:
        
        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$chi1*merged$weights)
        
        
        ## Perfrom analysis:
        
        n.annot <- 1
        
        
        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)
        
        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocks),ncol =(n.annot+1))
        colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
        replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
        colnames(xtx)<- colnames(weighted.LD)
        for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))}
        
        reg <- solve(xtx, xty)
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete, xty.delete)
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[,s] <- pseudo.values[,1]
        N.vec[1,s] <- N.bar
        
        if(is.na(pop.prev)==F & is.na(samp.prev)==F){
          conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
          Liab.S[j] <- conversion.factor
          LOG("     ", print = FALSE)
          LOG("Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait ",
              chi1, " when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.")
          LOG("     ", print = FALSE)
        }
        
        cov[j,j] <- reg.tot
        I[j,j] <- intercept
        
        lambda.gc <- median(merged$chi1) / qchisq(0.5, df = 1)
        mean.Chi <- mean(merged$chi1)
        ratio <- (intercept - 1) / (mean.Chi - 1)
        ratio.se <- intercept.se / (mean.Chi - 1)
        
        LOG("Heritability Results for trait: ", chi1)
        LOG("Mean Chi^2 across remaining SNPs: ", round(mean.Chi, 4))
        LOG("Lambda GC: ", round(lambda.gc, 4))
        LOG("Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")")
        LOG("Ratio: ", round(ratio, 4), " (", round(ratio.se, 4), ")")
        LOG("Total Observed Scale h2: ", round(reg.tot, 4), " (", round(tot.se, 4), ")")
        LOG("h2 Z: ", format(reg.tot / tot.se), digits = 3)
      }
      
      
      ##### GENETIC COVARIANCE code
      
      if(j != k){
        
        LOG("     ", print = FALSE)
        
        chi2 <- traits[k]
        LOG("Calculating genetic covariance [", s, "/", n.V, "] for traits: ", chi1, " and ", chi2)
        
        # Reuse the data read in for heritability
        y2 <- all_y[[k]]
        y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")], by = "SNP", sort = FALSE)
        
        y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x)
        y$ZZ <- y$Z.y * y$Z.x
        y$chi2 <- y$Z.y^2
        merged <- na.omit(y)
        n.snps <- nrow(merged)
        
        LOG(n.snps, " SNPs remain after merging ", chi1, " and ", chi2, " summary statistics")
        
        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        
        
        #### MAKE WEIGHTS:
        
        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N.x)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- pmax(merged$L2, 1)
        merged$w.ld <- pmax(merged$wLD, 1)
        merged$c <- tot.agg*merged$N.x/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        
        tot.agg2 <- (M.tot*(mean(merged$chi2)-1))/mean(merged$L2*merged$N.y)
        tot.agg2 <- max(tot.agg2,0)
        tot.agg2 <- min(tot.agg2,1)
        merged$ld2 <- pmax(merged$L2, 1)
        merged$w.ld2 <- pmax(merged$wLD, 1)
        merged$c2 <- tot.agg2*merged$N.y/M.tot
        merged$het.w2 <- 1/(2*(1+(merged$c2*merged$ld))^2)
        merged$oc.w2 <- 1/merged$w.ld2
        merged$w2 <- merged$het.w2*merged$oc.w2
        merged$initial.w2 <- sqrt(merged$w2)
        
        
        merged$weights_cov <- (merged$initial.w + merged$initial.w2)/sum(merged$initial.w + merged$initial.w2 )
        
        N.bar <- sqrt(mean(merged$N.x)*mean(merged$N.y))
        
        ## preweight LD and chi:
        
        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$ZZ *merged$weights_cov)
        
        ## Perfrom analysis:
        
        
        n.annot <- 1
        
        
        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)
        
        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocks),ncol =(n.annot+1))
        colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
        replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
        colnames(xtx)<- colnames(weighted.LD)
        for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))}
        
        reg <- solve(xtx, xty)
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete, xty.delete)
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[, s] <- pseudo.values[, 1]
        N.vec[1, s] <- N.bar
        
        cov[k, j] <- cov[j, k] <- reg.tot
        I[k, j] <- I[j, k] <- intercept
        
        LOG("Results for genetic covariance between: ", chi1, " and ", chi2)
        LOG("Mean Z*Z: ", round(mean(merged$ZZ), 4))
        LOG("Cross trait Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")")
        LOG("Total Observed Scale Genetic Covariance (g_cov): ", round(reg.tot, 4), " (", round(tot.se, 4), ")")
        LOG("g_cov Z: ", format(reg.tot / tot.se), digits = 3)
        LOG("g_cov P-value: ", format(2 * pnorm(abs(reg.tot / tot.se), lower.tail = FALSE), digits = 5))
      }
      
      ### Total count
      s <- s + 1
    }
    
    #mod addition
    gc() #do garbage collect if this can help with out of memory issues.
    
    if(limited) break
  }
  
  
  ## Scale V to N per study (assume m constant)
  # /!\ crossprod instead of tcrossprod because N.vec is a one-row matrix
  v.out <- cov(V.hold) / crossprod(N.vec * (sqrt(n.blocks) / m))
  
  ### Scale S and V to liability:
  ratio <- tcrossprod(sqrt(Liab.S))
  S <- cov * ratio
  
  #calculate the ratio of the rescaled and original S matrices
  scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
  
  #rescale the sampling correlation matrix by the appropriate diagonals
  V <- v.out * tcrossprod(scaleO)
  
  
  #name traits according to trait.names argument
  #use general format of V1-VX if no names provided
  colnames(S) <- if (is.null(trait.names)) paste0("V", 1:ncol(S)) else trait.names
  
  if(mean(Liab.S)!=1){
    r<-nrow(S)
    SE<-matrix(0, r, r)
    SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(V))
    
    LOG(c("     ", "     "), print = FALSE)
    LOG("Liability Scale Results")
    
    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j == k){
          LOG("     ", print = FALSE)
          LOG("Liability scale results for: ", chi1)
          LOG("Total Liability Scale h2: ", round(S[j, j], 4), " (", round(SE[j, j], 4), ")")
        }
        
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          LOG("Total Liability Scale Genetic Covariance between ", chi1, " and ",
              chi2, ": ", round(S[k, j], 4), " (", round(SE[k, j], 4), ")")
          LOG("     ", print = FALSE)
        }
      }
    }
  }
  
  #mod additions - initialise S_Stand, set all NA element in S to 0 so not to have NA values
  S_Stand<-matrix(NA, nrow(S), nrow(S))
  V_Stand<-matrix(NA, nrow(V), nrow(V))
  S[is.na(S)]<-0
  if(all(diag(S) > 0)) {
  
    ##calculate standardized results to print genetic correlations to log and screen
    ratio <- tcrossprod(1 / sqrt(diag(S)))
    S_Stand <- S * ratio
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
    
    ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in
    # -> not possible because of 'all(diag(S) > 0)'
    # scaleO[is.nan(scaleO)] <- 0
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    V_Stand <- V * tcrossprod(scaleO)
    
    #enter SEs from diagonal of standardized V
    r<-nrow(S)
    SE_Stand<-matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand,diag=TRUE)] <-sqrt(diag(V_Stand))
    
    
    LOG(c("     ", "     "), print = FALSE)
    LOG("Genetic Correlation Results")
    
    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          LOG("Genetic Correlation between ", chi1, " and ", chi2, ": ",
              round(S_Stand[k, j], 4), " (", round(SE_Stand[k, j], 4), ")")
          LOG("     ", print = FALSE)
        }
      }
    }
  
  } else {
    warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
    LOG("Your genetic covariance matrix includes traits estimated to have a negative heritability.", print = FALSE)
    LOG("Genetic correlation results could not be computed due to negative heritability estimates.")
  }
  
  end.time <- Sys.time()
  
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- floor(total.time-mins*60)
  
  LOG("     ", print = FALSE)
  LOG("LDSC finished running at ", end.time)
  LOG("Running LDSC for all files took ", mins, " minutes and ", secs, " seconds")
  LOG("     ", print = FALSE)
  
  flush(log.file)
  close(log.file)
  
  # mod additions - added the suggested computations of standard error matrices from the website
  rownames(S)<-colnames(S)
  S.SE<-matrix(0, nrow(S), nrow(S))
  colnames(S.SE)<-colnames(S)
  rownames(S.SE)<-colnames(S)
  S.SE[lower.tri(S.SE,diag=TRUE)] <-sqrt(diag(V))
  S.SE[upper.tri(S.SE)]<-t(S.SE)[upper.tri(S.SE)]
  if(stand){
    colnames(S_Stand)<-colnames(S)
    rownames(S_Stand)<-colnames(S)
    S_Stand.SE<-matrix(0, nrow(S_Stand), nrow(S_Stand))
    colnames(S_Stand.SE)<-colnames(S)
    rownames(S_Stand.SE)<-colnames(S)
    S_Stand.SE[lower.tri(S_Stand.SE,diag=TRUE)] <-sqrt(diag(V_Stand))
    S_Stand.SE[upper.tri(S_Stand.SE)]<-t(S_Stand.SE)[upper.tri(S_Stand.SE)]
  }
  
  if(stand){
    list(V=V,S=S,S.SE=S.SE,I=I,N=N.vec,m=m,V_Stand=V_Stand,S_Stand=S_Stand,S_Stand.SE=S_Stand.SE)
  } else {
    list(V=V,S=S,S.SE,S.SE,I=I,N=N.vec,m=m)
  }
}
