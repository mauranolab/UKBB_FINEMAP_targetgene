#!/usr/bin/env Rscript

library(data.table)
args<-commandArgs(trailingOnly = T)

## Read in gene set files
all_gencode_v24_file<-args[1] # Background set
targeted_genes_file<-args[2]  # Targeted gene set
pc_genes_file<-args[3]				# Positive control gene set

all_gencode_v24_genes<-scan(file=all_gencode_v24_file,what="character",na.strings = "NA",quiet = T)
targeted_genes<-scan(file=targeted_genes_file,what="character",na.strings = "NA",quiet = T)
pc_genes<-scan(file=pc_genes_file,what="character",na.strings = "NA",quiet = T)


## Ensure there are no values in the sets that are not in the background set
targeted_genes<-intersect(targeted_genes,all_gencode_v24_genes)
pc_genes<-intersect(pc_genes,all_gencode_v24_genes)

## Count genes in each group
targeted_genes_count<-length(targeted_genes)
pc_genes_count<-length(pc_genes)
all_gencode_v24_genes_count<-length(all_gencode_v24_genes)

targeted_pc_genes<-length(intersect(targeted_genes,pc_genes))


## numerator=proportion of targeted genes that are also positive control genes
numerator<-targeted_pc_genes/targeted_genes_count

## denominator=poroportion of all genes that are in positive control genes
denominator<-pc_genes_count/all_gencode_v24_genes_count

## Perform enrichment calculations
enrichment=round(numerator/denominator,digits = 2)

outdf<-data.frame(targeted_pc_genes,
									targeted_genes_count,
									pc_genes_count, 
									all_gencode_v24_genes_count,
                  enrichment)

write.table(outdf, file = "",append = T, quote = F,sep = "\t",row.names = F,col.names = F)

