#!/usr/bin/env Rscript

library(data.table)
args<-commandArgs(trailingOnly = T)

# Read in file paths
gene_universe_file<-args[1]
set1_file<-args[2]
set2_file<-args[3]

gene_universe<-scan(file=gene_universe_file,what="character",na.strings = "NA",quiet = T)
set1<-scan(file=set1_file,what="character",na.strings = "NA",quiet = T)
set2<-scan(file=set2_file,what="character",na.strings = "NA",quiet = T)

gene_universe<-na.omit(gene_universe)
set1<-na.omit(set1)
set2<-na.omit(set2)

# Ensure there are no values in the sets that are not in the background set
set1<-intersect(set1,gene_universe)
set2<-intersect(set2,gene_universe)

# Make enrichment calculations
set1_all<-length(set1)
set2_all<-length(set2)
gene_universe_all<-length(gene_universe)

in_both<-length(intersect(set1,set2))


# numerator=proportion of set 1 that is in set 2
numerator<-in_both/set1_all

# denominator=poroportion of all genes that are in set 2
denominator<-set2_all/gene_universe_all


outdf<-data.frame(in_both,set1_all,set2_all, gene_universe_all,
                  enrichment=numerator/denominator)

write.table(outdf, file = "",append = T, quote = F,sep = "\t",row.names = F,col.names = F)

