Commands for manually running the pipeline

```bash

#munge file with existing info in spreadsheet
sbatch --time 01:00:00 --partition cpu --job-name="munge" --ntasks 1 --cpus-per-task 5 --mem 64G --wrap "Rscript ~/project/genetic_correlations/scripts/cleaning/runStandardCleanAndMunge.R -f '/scratch/prj/gwas_sumstats/new/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz' -c 'BODY14' -r '/scratch/prj/gwas_sumstats/variant_lists/w_hm3.snplist.flaskapp2018' -o '/scratch/prj/gwas_sumstats/munged'" --output "BODY14.$(date +%Y%m%d).out.txt" ;

```
