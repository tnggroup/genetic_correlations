# Commands for manually running the pipeline

Clean file with existing info in spreadsheet
```bash
sbatch --time 01:00:00 --partition cpu --job-name="clean" --ntasks 1 --cpus-per-task 5 --mem 64G --wrap "Rscript ~/project/genetic_correlations/scripts/cleaning/runStandardCleanAndMunge.R -f '/scratch/prj/gwas_sumstats/new/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz' -c 'BODY14' -o '/scratch/prj/gwas_sumstats/cleaned' --filter.maf 0.001 --filter.info 0.6" --output "BODY14.$(date +%Y%m%d).out.txt"
```

Munge file with existing info in spreadsheet
```bash
sbatch --time 01:00:00 --partition cpu --job-name="munge" --ntasks 1 --cpus-per-task 5 --mem 64G --wrap "Rscript ~/project/genetic_correlations/scripts/cleaning/runStandardCleanAndMunge.R -f '/scratch/prj/gwas_sumstats/new/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz' -c 'BODY14' -r '/scratch/prj/gwas_sumstats/variant_lists/w_hm3.snplist.flaskapp2018' -o '/scratch/prj/gwas_sumstats/munged'" --output "BODY14.$(date +%Y%m%d).out.txt"

```

Munge file with existing info in spreadsheet, HC1kG reference
```bash
sbatch --time 01:00:00 --partition cpu --job-name="munge" --ntasks 1 --cpus-per-task 5 --mem 64G --wrap "Rscript ~/project/genetic_correlations/scripts/cleaning/runStandardCleanAndMunge.R -f '/scratch/prj/gwas_sumstats/new/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz' -c 'BODY14' -r '/scratch/prj/gwas_sumstats/variant_lists/hc1kgp3.b38.mix.l2.jz2023.gz' -o '/scratch/prj/gwas_sumstats/munged_hc1kg'" --output "BODY14.$(date +%Y%m%d).out.txt"

```
