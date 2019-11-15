#!/bin/bash

#set -e -o pipefail

module load bedops/2.4.37

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
alias closest-features='closest-features --header'

## Parameters
trait_snps=$1           # bed/starch file with -log10 p-val or log10 bf in column 5
trait_name=$2           # String
pc_gene_list=$3         # txt file. One column of gene IDs
pc_gene_list_name=$4    # String
all_genes=$5            # bed/starch file with gene ID in column 4
outfile_txt=$6          # String
sig_metric=$7           # String. Either "-log10_P-value" or "log10_Bayes_factor"
snp_set=$8              # String. Either "All_SNPs", "DHS_SNPs", or "Trait-Specific_DHS_SNPs"

MHC_coords_hg19="/vol/mauranolab/vulpen01/ukbb/src_data/gencodev24_hg19/MHC_coordinates.hg19.bed"
MHC_genes_gencodev24="/vol/mauranolab/vulpen01/ukbb/src_data/positive_control_genes/src/MHC_genes.txt"

source targeting_methods.sh
targeting_method="closest_gene" # see targeting_methods.sh


all_genes_list=$TMPDIR/${trait_name}_${pc_gene_list_name}.all_genes_list.txt # 1 column: Gene name (tab delimited)
targeted_genes=$TMPDIR/${trait_name}_${pc_gene_list_name}.targeted_genes.txt # 2 columns: Gene name, value (tab delimited)
positive_control_genes_noMHC=$TMPDIR/${trait_name}_${pc_gene_list_name}.positive_control_genes.txt



if [ $sig_metric == "-log10_P-value" ];then
    sig_cutoff=7.30103
    min=0
    max=30
fi

if [ $sig_metric == "log10_Bayes_factor" ];then
    sig_cutoff=2
    min=-10
    max=14
fi


# Pair each gene with the minimum distance to a significant SNP
if [ $sig_metric == "-log10_P-value" ];then
    bedops -n $trait_snps $MHC_coords_hg19 | \
    awk -F"\t" -v max=$max 'BEGIN{OFS="\t"}{
        log10p=-log($5)/log(10); 
        if(log10p!="inf"){
            print $1,$2,$3,$4, log10p
        }else{print $1,$2,$3,$4,max}}' | \
    awk -F"\t" -v sig_cutoff=$sig_cutoff '$5>=sig_cutoff' | \
    target_genes $targeting_method - $all_genes | \
    awk -F"\t" 'function abs(v) {return v < 0 ? -v : v} BEGIN{OFS="\t"}{if($NF!="NA"){print $4,abs($NF)}}' | \
    sort -nk2 | awk -F"\t" '!seen[$1]++' > $targeted_genes 
else
    bedops -n $trait_snps $MHC_coords_hg19 | \
    awk -F"\t" -v min=$min 'BEGIN{OFS="\t"}{
        if($5=="-Inf"){$5=min};
        print $1,$2,$3,$4,$5
    }' | \
    awk -F"\t" -v sig_cutoff=$sig_cutoff '$5>=sig_cutoff' | \
    target_genes $targeting_method - $all_genes | \
    awk -F"\t" 'function abs(v) {return v < 0 ? -v : v} BEGIN{OFS="\t"}{if($NF!="NA"){print $4,abs($NF)}}' | \
    sort -nk2 | awk -F"\t" '!seen[$1]++' > $targeted_genes  
fi




# Make input to enrichment script
unstarch $all_genes | cut -f 4 | sort | uniq > $all_genes_list

if [ -f $outfile_txt ];then rm $outfile_txt;fi
echo -e "distance\tcutoff\ttargeting_method\tsnp_set\tannotation\ttrait\tTargeted_PositiveControl\tTargeted\tAll_PositiveControl\tAll\tenrichment" > $outfile_txt


targeted_genes_subset=$TMPDIR/${trait_name}_${pc_gene_list_name}.targeted_genes_subset.txt

sequence=`seq 10000 10000 500000`

for value in ${sequence[@]};do
    cat $targeted_genes | awk -F"\t" -v cutoff=$value '$2<=cutoff{print $1}' | sort | uniq > $targeted_genes_subset
    echo -ne "$value\t${sig_metric}_${sig_cutoff}\t$targeting_method\t$snp_set\t$pc_gene_list_name\t$trait_name\t" >> $outfile_txt
    ./gsea.enrichment.R $all_genes_list $targeted_genes_subset $pc_gene_list >> $outfile_txt
done



