#!/bin/bash 
f=$1 # single file 
REF=$2 
OUTDIR=$3 

module load plink 

bname=$(basename $f)
orig_dir=$(dirname $f)

new_f=$OUTDIR/$bname 
sed 's/^23/X/g' ${f}.bim > ${new_f}_X_recoded.bim
#echo python /hpc/users/sloofl01/software/snpflip/bin/snpflip  --fasta-genome=$REF --bim-file=${new_f}_X_recoded.bim --output-prefix=$new_f.checkflip
#python /hpc/users/sloofl01/software/snpflip/bin/snpflip  --fasta-genome=$REF --bim-file=${new_f}_X_recoded.bim --output-prefix=$new_f.checkflip
echo python /sc/arion/projects/buxbaj01a/software/snpflip/bin/snpflip  --fasta-genome=$REF --bim-file=${new_f}_X_recoded.bim --output-prefix=$new_f.checkflip
python /sc/arion/projects/buxbaj01a/software/snpflip/bin/snpflip  --fasta-genome=$REF --bim-file=${new_f}_X_recoded.bim --output-prefix=$new_f.checkflip

# Missnp was made during merge and failed 
to_flip=${new_f}.checkflip.reverse
to_remove=${new_f}.checkflip.ambiguous 
#plink --bfile ${f} --exclude $to_remove --make-bed --allow-no-sex --out ${new_f}.excluded_tmp 
plink --bfile ${f} --exclude $to_remove --make-bed --allow-no-sex --out ${new_f}.excluded 
#plink --bfile ${new_f}.excluded_tmp --exclude /sc/arion/projects/rg_choj07/sloofl01/covid19/MERGED/data/seaver_sema4_regen_EUR-merge.missnp --make-bed --allow-no-sex --out ${f}.excluded 
plink --bfile ${new_f}.excluded --flip $to_flip --allow-no-sex --make-bed --out ${new_f}.final
