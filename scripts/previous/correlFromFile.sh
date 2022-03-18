input1=$1
input2=$2
output=$3
ldscoredir=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/usr/helenaG/ldsc-master
###pythonpath="/users/k1507306/Python-2.7.11/python"
pythonpath="/users/$(whoami)/.conda/envs/myenv/bin/python"

$pythonpath $ldscoredir/ldsc.py  --n-blocks 200 --print-delete-vals --rg ${input1},${input2} --ref-ld-chr $ldscoredir/eur_w_ld_chr/ --w-ld-chr $ldscoredir/eur_w_ld_chr/ --out ${output} 
gawk 'c&&c--;/^p1.*p2.*gcov_int_se$/{c=1}' ${output}.log

mv ${output}*delete /mnt/lustre/groups/ukbiobank/sumstats/deletevalues/

