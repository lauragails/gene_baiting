#!/bin/bash 
plink_ref=$1 
pvals=$2 
N=$3 
annot=$4 
set_file=$5 
OUTDIR=$6 

# Go to the right place 
ml magma_gwas/1.0.8
cd $OUTDIR

# Get file prefixes   
bname_pvals=$(basename $pvals .forMAGMA)
bname_set=$(basename $set_file .txt)

magma --bfile $plink_ref --pval $pvals  N=$N --gene-annot $annot --out $bname_pvals --gene-settings adap-permp 

#third, to perform a gene set analysis
#example
magma --gene-results ${bname_pvals}.genes.raw --set-annot $set_file --out ${bname_pvals}_${bname_set}
