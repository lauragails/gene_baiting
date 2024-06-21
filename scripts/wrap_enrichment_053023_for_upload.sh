#!/bin/bash
# Script can be run all in R 
# Step 1: Get all files for mouse (filepaths below edited) 
bait_file="files/bait_lists_053023.tsv"
fu_file="files/Fu_ASD_published_mouse_FDR_lt_0.05.txt"
sat_file="files/Satterstrom_published_102_FDR_lt_0.1.txt"

outdir="/path/to/your/outdir/"
dir.create(outdir) 

# Background genes is mouse only! Only use mouse inputs and convert ASD list a priori  
remove_file=FALSE 
external_universe="files/Table_S5_Statistical_Background_053023_singlecol.csv" 

# Step 2: Find enrichment stats in 3 Fu lists 
script=scripts/resubmission_enrichment_baiting_NO_neg_control_053023.R
n_permutation=20 
OUTDIR=/path/to/outdir/
mkdir -p $OUTDIR  

for file in $fu_file $sat_file 
    do
	bname=$(basename $file .txt)
	outfile=$OUTDIR/final_baiting_results_${bname}.txt 

         Rscript $script $bait_file $file $n_permutation $remove_file mouse $external_universe $outfile # commented code where string mattered
    done 

# Convert to excel (optional but was helpful!)
script=scripts/dir_to_excel.R
INDIR=/path/to/results
bname=baiting_results_053023

Rscript $script $INDIR .txt $bname

# in R 
source("scripts/helper_funcs_baiting_053023.R")
DIR="/path/to/results"
files=Sys.glob(paste0(DIR,"/*txt")) 
for (f in files){
 df=read.table(f,header=TRUE)
 bname=gsub(".txt","",basename(f))
 outfile=paste0(DIR,"/",bname,".png") 
 png(outfile)
 par(mar=c(8, 10, 3, 8)) 
 basic_pval_plot(df)
 dev.off() 
}

# Checking overlap between Satterstrom and Fu lists
library(UpSetR)
fu_file="files/Fu_ASD_published_mouse_FDR_lt_0.05.txt"
sat_file="files/Satterstrom_published_102_FDR_lt_0.1.txt"
external_universe="files/Table_S5_Statistical_Background_053023_singlecol.csv" 
 
fu=read.table(fu_file,header=TRUE)
sat=read.table(sat_file,header=TRUE) 
u=read.table(external_universe,header=TRUE) 

lists=c()
lists[["fu"]]=unlist(fu) 
lists[["sat"]]=unlist(sat)
lists[["univ"]]=unlist(u) 

upset(fromList(lists),order.by="freq")
