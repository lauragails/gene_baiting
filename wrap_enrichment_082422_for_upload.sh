#!/bin/bash
# Script can be run all in R 
# Step 1: Get all files for mouse (filepaths below edited) 
bait_file="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/files/mouse/mouse_baits_no_filter_copied_082422.txt"
asd_human="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/scratch/Fu_ASD_published_human_FDR_lt_0.05.txt"
asd_mouse="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/files/mouse/Fu_ASD_published_mouse_FDR_lt_0.05.txt"

# For simplicity get ASD list in bash and convert using my_awk 
awk '$14<0.05 {print $1}' /sc/arion/projects/buxbaj01a/sloofl01/ASC/files/gene_lists/Fu_2021_published/Sup_tab_11_from_copypaste.txt > $asd_human 

# Convert to mouse 
# Downloaded 082422 http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt 
key=/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/files/mouse/mouse_human_homologs.txt 
~/scripts/my_awk.pl $asd_human $key | cut -f3 > $asd_mouse 

# Human ASD lists from preprint  
thresh=0.01 
#outdir="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/results/082622/"
outdir="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/results/082422/"
dir.create(outdir) 

# Background genes is mouse only! Only use mouse inputs and convert ASD list a priori  
bait_file="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/files/mouse/mouse_baits_no_filter_copied_082422.txt"
asd_mouse="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/files/mouse/Fu_ASD_published_mouse_FDR_lt_0.05.txt"
remove_file="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/files/mouse/SolTid_filter_mouse.txt"
external_universe="/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/files/mouse/universe_mouse_recopied_082422"

# Step 2: Find enrichment stats in 3 Fu lists 
script=/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/scripts/run_enrichment_baiting_with_neg_control_082422.R
thresh=0.05 
n_permutation=10 
outfile=/sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/results/082622/final_baiting_results.txt

mkdir -p /sc/arion/projects/buxbaj01a/sloofl01/ASC_Baiting/results/082622/
Rscript $script $bait_file $asd_mouse $thresh $n_permutation $remove_file mouse $external_universe $outfile 
