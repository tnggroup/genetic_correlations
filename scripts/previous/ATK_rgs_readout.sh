################################
############# Software to run to create csv file #######

###input paths
INPUTPATH=/filepath/rGs_LDSC/output/results/correlations/CTS
OUTPUTPATH=/filepath/rGs_LDSC/output/results/correlations/CTS/summary
CODE=*
OUTPUTPHENOTYPE=TRAITNAME #e.g. CTS
OUTPUTFILENAME=TRAITNAME_results

# fgrep accepts a fixed pattern (e.g. _noMHC.sumstats.gz  /scratch/groups/ukbiobank/sumstats/munged_noMHC/) to match, but then splits that pattern up into one search string per line. the pattern here is the end of phenotype1 and start of phenotype2 in the rG table output at the end of the log file.
LANG=C fgrep -h "CTS_noMHC.sumstats.gz  /scratch/groups/ukbiobank/sumstats/munged_noMHC/" ${INPUTPATH}/${CODE}_*.log > ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs_temp.txt

awk 'BEGIN {OFS=","} NR==1 {print "p1,p2,rg,se,z,p,h2_obs,h2_obs_se,h2_int,h2_int_se,gcov_int,gcov_int_se"} {print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs_temp.txt > ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs.csv

LANG=C fgrep -h "SNPs with valid alleles" ${INPUTPATH}/${CODE}_*.log /dev/null > ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs_n_snps.txt

rm ${OUTPUTPATH}/*temp*

# Delete the folder path in the file
sed -i 's/\/scratch\/groups\/ukbiobank\/sumstats\/munged_noMHC\///g' ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs.csv
head -2 ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs.csv

# Delete the sumstats end
sed -i 's/_noMHC.sumstats.gz//g' ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs.csv
head -2 ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs.csv | column -t

# Delete the folder path in the file for CTS
sed -i 's/\/scratch\/groups\/ukbiobank\/usr\/abi\/PhD\/projects\/CTS_gSEM\/data\/munged\///g' ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs.csv
head -2 ${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs.csv

# Open new terminal with local folder

# Repeat for local hardrive terminal
# Output folder
OUTPUTPATH=/filepath/rGs_LDSC/output/results/correlations/CTS/summary
echo $OUTPUTPATH
# output file and phentoype names
OUTPUTPHENOTYPE=TRAITNAME
OUTPUTFILENAME=TRAITNAME_results
echo $OUTPUTFILENAME

# find Folder on hard drive and pwd to copy filepath in next step
cd /filepath/results
pwd

# Copy file to hard drive from rosalind (insert your own k number)
scp Knumber@login3.rosalind.kcl.ac.uk:/${OUTPUTPATH}/${OUTPUTPHENOTYPE}_${OUTPUTFILENAME}_rgs.csv /filepathHarDrive/results
