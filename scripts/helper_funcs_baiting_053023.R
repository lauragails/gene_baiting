load_libraries=function(){
    library("tools") 
    library("org.Hs.eg.db")
    library(HGNChelper)
    library("biomaRt") 
}
get_first_item=function(x,delim=" "){
    str=unlist(strsplit(x,delim))
    item=str[1]
    return(item)
}

wrap_symbol_synonyms=function(vec,species="mouse"){
# Symbol synonyms (below) gets most common gene name/standardizes gene names for given Gene ID inputs. 
# It does NOT fix excel dates. Likewise, there were a few genes (ie with "orf") that were returned as NA, but had an exact string match in both dataset. This wrapper takes care of that scenario. 

# 1. Change the dates
    vec=unlist(vec)
    print(head(vec))
    tbl=checkGeneSymbols(vec, unmapped.as.na = FALSE, map = NULL, species = species)
    print(dim(tbl)) 
    print(length(unlist(vec)))
    stopifnot(nrow(tbl)==length(unlist(vec))) 
    #vec=noquote(strsplit(noquote(tbl$Suggested.Symbol),","))
    vec=tbl$Suggested.Symbol
    print(tail(vec))
    vec=apply(data.frame(vec),1,get_first_item)
    print(tail(vec))
    vec=unlist(vec)
    # Get first item if there are many 

    #vec=symbol_synonyms(unlist(vec),species) # Commenting for now, causing a lot of trouble! Let's see what happens without it 
    stopifnot(nrow(tbl)==length(unlist(vec))) 
    return(vec) # this WILL return NAs, we want to return vector of the same length  
}

# Gene synonyms (copied from COBS paper, Lin 2021) 
# The code for this comes from the user Duff at \url{https://www.biostars.org/p/14971/}.
symbol_synonyms=function(vec, species="mouse",verbose = T){
  dbCon=c() 
  if (species=="mouse"){
      dbCon=org.Mm.eg.db::org.Mm.eg_dbconn()
  } else if (species=="human"){
      dbCon=org.Hs.eg.db::org.Hs.eg_dbconn()
  } else {
      print("species isn't listed in symbol_synonyms")
  }

  sqlQuery='SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  aliasSymbol=DBI::dbGetQuery(dbCon, sqlQuery)

  sapply(1:length(vec), function(i){
    if(verbose & i %% max(floor(length(vec)/10),1) == 0) cat('*')
  
    res=aliasSymbol[which(aliasSymbol[,2] %in% vec[i]), 5]
    if(length(res) == 0) return(NA)

    #if there are more than one, take the most common one
    if(length(res) > 1){
      len_vec=sapply(res, function(x){
        length(which(aliasSymbol[,2] %in% x))
      })

      #if there is still a tie, take the first alphabetic one
      res_final=res[which(len_vec == max(len_vec))]
      if(length(res_final) > 1){
        res_final=sort(res_final, decreasing = F)[1]
      }
      res=res_final
    }
    res
  })
}

# Basic function to convert human to mouse gene names (I did NOT end up using this) 
# Source: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/ 
convertHumanGeneList <- function(x){
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
}

# User defined functions 
my_hypergeometric=function(universe,FDR_sig,bait_list){
    N=length(unique(universe)) 
    m=length(unique(bait_list))
    n=length(unique(FDR_sig))
    k_genes=unique(intersect(FDR_sig,bait_list)) 
    k=length(k_genes)

    pval=phyper(k-1,m,N-m,n,lower.tail= FALSE)
    result=c()
    result$N=N
    result$n=n
    result$m=m
    result$k=k
    result$pval=pval
    return(result)
}

calc_qval=function(universe,FDR_sig,bait_list,n_permutation,thresh){
    # NOTE: Using R's FDR correction, I wrote this as a sanity check. 
    counter=0
    for (i in 1:n_permutation){
        rand=names(sample(unique(universe),length(bait_list)))
        df=my_hypergeometric(universe,FDR_sig,rand)
	pval_rand=df$pval
	if (pval_rand<thresh){counter=counter+1}
    }
    avg_FP=counter/n_permutation # aka q-value 
    return(avg_FP) 
}

read_list=function(list_file,species="mouse",to_remove=NA,header=FALSE){
    vec=read.table(list_file,header=header,sep="\t") # Laura only works with tabs, and ASC list needs that
    vec=unlist(vec) 
    #vec=wrap_symbol_synonyms(vec,species) # just using mouse input 
    
    # Read genes to remove, if needed 
    remove=c()
    keep=which(is.na(vec)==FALSE)
    if (length(to_remove)>0 && is.na(to_remove[1])==FALSE){
        remove=union(remove,NA)
        remove_idx=which(vec %in% unlist(remove)) # includes NA  
        keep=setdiff(keep,remove_idx) # which is NOT NA 
    }
    vec=unique(vec[keep])
    return(vec)
}

wrap_hypergeometric=function(universe,baits,Fu_list,to_remove=NA,species="mouse",n_permutation=10){
# Do enrichments for baits 
result=c()
for (my_list in 1:ncol(baits)){
    # Process gene list 
    l=baits[,my_list]
    #l=wrap_symbol_synonyms(l,species) 
    l=unique(setdiff(l,NA))
    l_prev=l 

    # Make sure all lists are in the universe. Note to_remove is *only* from baiting  
    l=setdiff(l,to_remove)
    Fu_list=intersect(unlist(Fu_list),unlist(universe)) 

    # Get names 
    bname=colnames(baits)[my_list]

    # Hypergeometric test 
    output=my_hypergeometric(universe,Fu_list,l) # all lists are unique-d in function 
    N=output$N
    m=output$m
    n=output$n
    k=output$k
    pval=output$pval
    qval=calc_qval(universe,Fu_list,l,n_permutation,pval)

    # Finish getting enrichment test params  
    m_prev=length(l_prev) # see how many genes aren't in universe 
    pcnt_loss=round((m-m_prev)*100/m_prev,2) 

    # return result
    rslt=c(bname,list_name,N,n,m,k,m_prev,pcnt_loss,pval,qval,n_permutation)
    result=rbind(result,rslt) 
}

result=data.frame(result)
colnames(result)=c("bait_list","gene_list","N","n","m","k","m_prev","pcnt_loss","pval","qval","n_permutations")

# add fdr etc for ***within** a single test/dataframe   
for (method in p.adjust.methods){
        result[,c(method)]=p.adjust(result$pval,method=method)
}
return(result)
}

basic_pval_plot=function(df,alpha=0.05,ntests=14,my_order=c("Anks1b","Syngap1","Shank2","Shank3","Nckap1","Nbea","Ctnnb1","Lrrc4c","Iqsec2","Arhgef9","Ank3","Scn2a","Scn8a","Hnrnpu")){
    # example data
    # ordered_df=df[order(df$pval,decreasing=TRUE),]

    # Table has "bait" in it. Fix my_order kludgily 
    my_order=gsub(".Bait","",my_order)
    my_order=paste0(my_order,".Bait")
    my_order=rev(unlist(as.character(my_order)))
    rownames(df)=df$bait_list 
    #ordered_df=df[order(df$bait_list==my_order),]
    ordered_df=df[my_order,]
    bonf_thresh=-log10(alpha/ntests) 
    baits=gsub(".Bait","",ordered_df$bait_list)

    # Make prettier title. Can add as many rows as necessary, strings should be unique 
    ln=as.character(unlist(unique(df$gene_list)))
    main_sub="Overlap with Bait Lists"
    if (ln=="Satterstrom_published_102_FDR_lt_0.1"){
	main_sub="Overlap with Satterstrom et al, 2020"
    }
    if (ln=="Fu_ASD_published_mouse_FDR_lt_0.05"){
	main_sub="Overlap with Fu et al, 2022" 
    }
    # create barplot
    #plot.new()
    barplot(-log10(ordered_df$pval), horiz = TRUE, las = 1, xlab = expression(-log10(P[hypergeometric])), ylab="",names.arg = baits, main = main_sub,col="deepskyblue3",cex.names=1.3)

    # add vertical lines at different thresholds
    abline(v = bonf_thresh, col = "firebrick1", lwd = 2)
    abline(v = -log10(alpha), col = "chocolate4", lwd = 2, lty="dashed")
}
