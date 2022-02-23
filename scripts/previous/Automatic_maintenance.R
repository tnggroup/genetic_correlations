#!/usr/bin/Rscript
# Author: Oliver Pain (oliver.pain@kcl.ac.uk)
suppressMessages(library("optparse"))

option_list = list(
  make_option("--old_csv", action="store", default=NA, type='character',
  		help="Path to old .csv to be updated [optional]"),
  make_option("--output", action="store", default=NA, type='character',
  		help="File name of output [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(easyPubMed)))

# Tasks:
# 1. Download a csv from Helena's SQL database
# 2. Insert additional information:
#   - LDSC statistics (h2)
#   - Enough for PRS information
#   - Title, abstract, link to reference from PubMed
#   - Whether PRS in UK Biobank is available
# 3. Generate PRS in UK biobank if not present

####################

# Download Helena's .csv
#cat('Downloading up to date GWAS repository data...')
#system(paste0("rm ",opt$output_dir,"download_csv"), ignore.stdout = T, ignore.stderr = T)
#system(paste0("wget http://10.200.105.5/download_csv -P ",opt$output_dir), ignore.stdout = T, ignore.stderr = T)

cat('Reading up to date GWAS repository data...')
csv_new<-fread(paste0(opt$output_dir,"download_csv"))
#system(paste0("rm ",opt$output_dir,"download_csv"), ignore.stdout = T, ignore.stderr = T)

cat('Done!\n')

# While testing reduce number of entries
csv_new<-csv_new[1:30,]

if(is.na(opt$old_csv) == F){
  cat('Read in old GWAS respository data.\n')
  # Read in old csv to idenitfy new entries
  csv_old<-fread(opt$old_csv)
  new_gwas<-csv_new[!(csv_new$code %in% csv_old$code),]
  
  # Insert check for corresponding LDSC munged file (if not present the upload may be incomplete)
  if(sum(!(csv_new$code %in% csv_old$code)) == 0){
    cat('No new entries detected.\n')
    q()
  } else {
    cat(dim(new_gwas)[1],'new entries detected.\n')
    cat('Downloading PubMed information for new entries...\n')
  }
} else {
  # If no old csv provided, continue with all entries in csv_new
  new_gwas<-csv_new
  cat('Downloading PubMed information for all entries...\n')
}

# Download PubMed information for new entries
new_gwas$article_title<-NA
new_gwas$article_abstract<-NA

for(i in which(!is.na(new_gwas$pmid))){
  pmid_temp<-get_pubmed_ids(new_gwas$pmid[i], api_key = NULL)
  dat<-fetch_pubmed_data(pmid_temp, format='xml')
  new_gwas$article_title[i] <- custom_grep(dat, "ArticleTitle", "char")
  new_gwas$article_abstract[i] <- custom_grep(dat, "AbstractText", "char")
  cat(paste0(round(i/sum(!is.na(new_gwas$pmid)),2)*100,"% complete\n"))
}

# Check whether there is sufficient data is available for PRS
new_gwas$enough_for_PRS<-NA

for(i in 1:dim(new_gwas)[1]){    
  new_gwas_i_log<-read.fwf(paste0('/mnt/lustre/groups/ukbiobank/sumstats/cleaned_logs/',new_gwas$code[i],'.log'),widths = 1000000)
  if('ENOUGH_FOR_PRS' %in% new_gwas_i_log$V1){
    new_gwas$enough_for_PRS[i]<-T
  } else {
    new_gwas$enough_for_PRS[i]<-F
  }
}

# Check whether PRS has been calculated in UK biobank
new_gwas$PRS_created<-NA
PRS_calculated<-list.files(path='/scratch/groups/ukbiobank/sumstats/PRS/ukb18177_glanville/PRS_for_use', pattern='.all.score')
PRS_calculated<-gsub('_header.all.score','',PRS_calculated)

for(i in 1:dim(new_gwas)[1]){
  if(new_gwas$uk_biobank[i] == 'no'){
    if(!(new_gwas$code[i] %in% PRS_calculated)){
      new_gwas$PRS_created<-F
    } else {
      new_gwas$PRS_created<-T
    }
  } else {
    new_gwas$PRS_created<-F
  }
}

# Look up LDSC stats
new_gwas$ldsc_h2_observed<-NA
new_gwas$ldsc_h2_se<-NA
new_gwas$ldsc_lambda_GC<-NA
new_gwas$ldsc_mean_chi2<-NA
new_gwas$ldsc_intercept<-NA
new_gwas$ldsc_intercept_se<-NA
new_gwas$ldsc_ratio<-NA
new_gwas$ldsc_ratio_se<-NA

ldsc_logs<-list.files(path='/scratch/groups/ukbiobank/sumstats/munged', pattern='_herit.log')
ldsc_logs<-gsub('_herit.log','',ldsc_logs)

for(i in 1:dim(new_gwas)[1]){
  if((new_gwas$code[i] %in% ldsc_logs)){
    ldsc_logs_i<-read.fwf(paste0('/scratch/groups/ukbiobank/sumstats/munged/',new_gwas$code[i],'_herit.log'),widths = 1000000)
        
    new_gwas$ldsc_h2_observed[i]<-as.numeric(gsub(' .*','', gsub('Total Observed scale h2: ','',ldsc_logs_i[grepl('Total Observed scale h2: ', ldsc_logs_i$V1),])))
    new_gwas$ldsc_h2_se[i]<-as.numeric(gsub("\\)",'', gsub(".*\\(",'', gsub('Total Observed scale h2: ','',ldsc_logs_i[grepl('Total Observed scale h2: ', ldsc_logs_i$V1),]))))
    new_gwas$ldsc_lambda_GC[i]<-as.numeric(gsub('Lambda GC: ','',ldsc_logs_i[grepl('Lambda GC: ', ldsc_logs_i$V1),]))
    new_gwas$ldsc_mean_chi2[i]<-as.numeric(gsub("Mean Chi\\^2: ",'',ldsc_logs_i[grepl("Mean Chi", ldsc_logs_i$V1),]))
    new_gwas$ldsc_intercept[i]<-as.numeric(gsub(' .*','', gsub("Intercept: ",'',ldsc_logs_i[grepl("Intercept: ", ldsc_logs_i$V1),])))
    new_gwas$ldsc_intercept_se[i]<-as.numeric(gsub("\\)",'', gsub(".*\\(",'', gsub("Intercept: ",'',ldsc_logs_i[grepl("Intercept: ", ldsc_logs_i$V1),]))))
    
    if((new_gwas$ldsc_intercept[i]-1)/(new_gwas$ldsc_mean_chi2[i]-1) > 0){
      new_gwas$ldsc_ratio[i]<-as.numeric(gsub(' .*','', gsub("Ratio: ",'',ldsc_logs_i[grepl("Ratio: ", ldsc_logs_i$V1),])))
      new_gwas$ldsc_ratio_se[i]<-as.numeric(gsub("\\)",'', gsub(".*\\(",'', gsub("Ratio: ",'',ldsc_logs_i[grepl("Ratio: ", ldsc_logs_i$V1),]))))
    } else {
      new_gwas$ldsc_ratio[i]<-'<0'
      new_gwas$ldsc_ratio_se[i]<-'NA'
    }
  }
}

# Write out updated csv
if(is.na(opt$old_csv) == F){
  new_gwas<-rbind(csv_old, new_gwas)
}

cat('Done.\n')

# Save updated .csv
write.csv(new_gwas, paste0(opt$output,'.',Sys.Date(),'.csv'), row.names=F, quote=T)


