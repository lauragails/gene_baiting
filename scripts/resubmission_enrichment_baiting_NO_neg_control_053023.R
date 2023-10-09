args = commandArgs(trailingOnly=TRUE) 
bait_file=args[1]
list_file=args[2]
n_permutation=args[3]
remove_file=args[4]
species=args[5]
external_universe=args[6]
outfile=args[7]

# Load libraries and make n_permutation numeric 
source("scripts/helper_funcs_baiting_053023.R")
load_libraries()
n_permutation=as.numeric(as.character(n_permutation))

# Read in files 
to_remove=NA
if (remove_file != FALSE){
    to_remove=read_list(remove_file)
}

universe=read_list(external_universe,species=species,to_remove=to_remove)

# Read tables 
baits=read.table(bait_file,sep="\t",header=TRUE)

# Read in ASC list - I did the mouse conversion a priori to make things easier  
list_name=gsub(".txt","",basename(list_file))
my_list=read.table(list_file,header=FALSE)

# Get the genes significant in Fu, and run hypergeometric testing  
result=wrap_hypergeometric(universe,baits,my_list,to_remove=to_remove,species=species,n_permutation=n_permutation)
write.table(result,file=outfile,sep="\t",quote=FALSE,row.names=FALSE)
