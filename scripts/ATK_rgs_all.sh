### RUNNING GENETIC CORRELATIONS FOR ALL (or multiple) TRAITS IN THE UKB SUMSTATS RESPOSITORY EXAMPLE SCRIPT


# Step 1: MAKE LIST OF SUMSTATS FOR GENETIC CORRELATIONS
# this can be done in two ways
# a) create a list of all of the sumstats availble in a folder in the ukb sumstats repository. Use the munged_no_mhc folder for traits with the mhc region removed.
# b) create a specific list of traits if you only want to analyse certain traits in the repository

# a) create list with all traits in folder e.g. munged_noMHC
    ls munged_noMHC > filepath/munged_noMHC.txt

    #remove end of file text from the list e.g. '_noMHC.sumstats.gz'

    sed 's/_noMHC.sumstats.gz//g' munged_noMHC_list.txt > munged_noMHC_list_290320.txt

# b) create a list of specific traits from the repository, for example:
echo 'RISK01
RISK02
RISK03
SCBP01
SCBP02
SCHI01
SCHI02
' > filepath/list_rg_070721.txt

# MAKE RESULTS DIRECTORIES - INCLUDING ONE FOR DELETE VALUES:

mkdir -p filepath/rGs_LDSC/output/results/correlations/
mkdir -p filepath/rGs_LDSC/output/results/correlations/deletevalues/

#STEP 2: Create genetic correlations script 1 (which will be used in the second rG script below in Step 3)


cat <<'EOT'>> filepath/scripts_folder/correlFromFile.sh
#!/bin/bash -l

#SBATCH -p shared,brc
#SBATCH --mem 15G
#SBATCH --nodes=1
#SBATCH -t 24:00:00
#SBATCH -J rGs_all

#specified at end of ldsc.correlations.sh after correlFromFile.sh is called in, where $1 is munged_nomhc sumstats repo file for each trait, $2 specific single trait, $3 file path of output for each rg log
input1=$1
input2=$2
output=$3
ldsc=filepath/ldsc
python=/scratch/groups/ukbiobank/Edinburgh_Data/Software/anaconda2/bin/python2.7

$python $ldsc/ldsc.py \
--n-blocks 200 \
--print-delete-vals \
--rg ${input2},${input1} \ #note order of input1 and 2 is important here. The second trait listed here (p2) is from the list of traits (input1) to get heritability results for each trait from the list
--ref-ld-chr $ldsc/eur_w_ld_chr/ \
--w-ld-chr $ldsc/eur_w_ld_chr/ \
--out ${output}

gawk 'c&&c--;/^p1.*p2.*gcov_int_se$/{c=1}' ${output}.log

mv ${output}*delete /scratch/groups/ukbiobank/usr/abi/PhD/projects/CTS_gSEM/rGs_LDSC/output/results/correlations/deletevalues/
#moves the delete values files into delete vaulues folder

EOT

#STEP 3: Create Genetic correlations script 2
cat <<'EOT'>> filepath/scripts_folder/ldsc.correlations.sh

#!/bin/bash -l

#SBATCH -p shared,brc
#SBATCH --mem 15G
#SBATCH --nodes=1
#SBATCH -t 24:00:00
#SBATCH -J rGs_all_CTS

#name1 and 2 specified in command line when job submitted
#name1: munged and MHC cleaned trait from list (txt file with list of traits made in step 1)
#name2: munged and mhc cleaned trait of choice that will run rg with all traits from list (name of specific trait e.g. ADHD01)

name2=$2

mkdir -p filepath/rGs_LDSC/output/results/correlations/$name2

varid=$1
pathvarid="filepath/scripts_folder/${varid}"
countvarid=$(wc -l ${pathvarid}|awk '{printf $1}')

for i in `seq 1 ${countvarid}`
do name1=$(gawk -v myvar=$i 'NR==myvar{printf $1}' ${pathvarid})

sh filepath/scripts_folder/correlFromFile.sh \
/scratch/groups/ukbiobank/sumstats/munged_noMHC/${name1}_noMHC.sumstats.gz \
filepath/data/munged/${name2}_noMHC.sumstats.gz \
filepath/rGs_LDSC/output/results/correlations/${name2}/${name1}_${name2}

done

EOT

####STEP 4: SUBMIT SCRIPT 2 (WHICH WILL ALSO SUBMIT SCRIPT 1 USING SH):
cd filepath/scripts_folder/
sbatch -p brc ldsc.correlations.sh munged_noMHC_list_290320.txt TRAITNAME ## script 1, $1 = list of sumstats, $2 = trait name e.g. ADHD01
