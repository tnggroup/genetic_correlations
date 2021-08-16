##### SCRIPT TO REMOVE MHC REGION FROM GWAS SUMMARY STATISTICS
#### NB: If adding mhc removed traits to the sumstats/munged_noMHC folder, please check first that the trait you want the mhc region removed is not already there to avoid duplicates

### Removing mhc region for a single trait

name=TRAITNAME
echo $name

cd /scratch/groups/ukbiobank/sumstats
if [ -f munged/${name}.sumstats.gz ]; then zcat munged/${name}.sumstats.gz > tmp; gawk -F"\t" '{if(f==1){r[$1]}else if(!($1 in r)){print}}' f=1 scripts/hapmap_hmc_uniq.txt f=2 tmp > munged_noMHC/${name}_noMHC.sumstats; gzip munged_noMHC/${name}_noMHC.sumstats; fi;


### Loop to remove mhc region for a list of traits in the /ukbiobank/sumstats folder that have not yet had mhc region removed

varid=for_mhc_removal.txt #text file that contains list of traits in the sumstats folder that are missing from /scratch/groups/ukbiobank/sumstats/munged_noMHC folder. This can be created using echo command
echo $varid

pathvarid="/filepath/lists/${varid}"
echo $pathvarid

countvarid=$(wc -l ${pathvarid}|awk '{printf $1}')
echo $countvarid

cd /scratch/groups/ukbiobank/sumstats

for i in `seq 1 ${countvarid}`
do name=$(gawk -v myvar=$i 'NR==myvar{printf $1}' ${pathvarid})
if [ -f munged/${name}.sumstats.gz ]; then zcat munged/${name}.sumstats.gz > tmp; gawk -F"\t" '{if(f==1){r[$1]}else if(!($1 in r)){print}}' f=1 scripts/hapmap_hmc_uniq.txt f=2 tmp > munged_noMHC/${name}_noMHC.sumstats; gzip munged_noMHC/${name}_noMHC.sumstats; fi;

done
