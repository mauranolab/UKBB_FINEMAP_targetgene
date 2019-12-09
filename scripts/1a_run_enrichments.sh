#!/bin/bash


module load bedops/2.4.37

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
alias closest-features='closest-features --header'

## Parameters
trait_snps=$1           # bed/starch file with p-val or log10 bf in column 5
                        # trait_snps file should already be filtered to remove MHC and Missense SNPs
trait_name=$2           # String
sig_metric=$3           # String. Describes column 5 of "trait_snps". Either "-log10_P-value" or "log10_Bayes_factor"
snp_set=$4              # String. Describes "trait_snps". Either "All_SNPs", "DHS_SNPs", or "Trait-Specific_DHS_SNPs"
pc_gene_list=$5         # txt file. One column of gene IDs
pc_gene_list_name=$6    # String.
all_genes=$7            # bed/starch file with gene ID in column 4. Contains all genes in transcript model



mkdir -p "intermediate_files"
outfile="intermediate_files/${trait_name}-${pc_gene_list_name}.${sig_metric}_${snp_set}.txt" # 11 columns (tab delimited)
all_genes_list=$TMPDIR/${trait_name}_${pc_gene_list_name}.all_genes_list.txt # 1 column: Gene name
targeted_genes=$TMPDIR/${trait_name}_${pc_gene_list_name}.targeted_genes.txt # 2 columns: Gene name, value (tab delimited)


function target_genes {
    snps=$1
    genes=$2

    # Returns two columns: 1. gene name and 2. distance to nearest SNP

    closest-features --closest --delim "\t" --dist $genes $snps | \
    awk -F"\t" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$NF}' | \
    sort-bed - | uniq | \
    awk -F"\t" 'function abs(v) {return v < 0 ? -v : v} BEGIN{OFS="\t"}{if($NF!="NA"){print $4,abs($NF)}}' | \
    sort -nk2 | awk -F"\t" '!seen[$1]++'
}


# Pair each gene with the minimum distance to a significant SNP
if [ $sig_metric == "P-value" ];then
    sig_cutoff=7.30103
    max=30

    unstarch $trait_snps | \
    grep -v "#" | \
    awk -F"\t" -v max=$max 'BEGIN{OFS="\t"}{
        log10p=-log($5)/log(10); 
        if(log10p=="inf"){log10p=max}
        print $1,$2,$3,$4, log10p
    }' | \
    awk -F"\t" -v sig_cutoff=$sig_cutoff '$5>=sig_cutoff' | \
    target_genes - $all_genes > $targeted_genes 

    sig_metric="-log10_P-value"

elif [ $sig_metric == "log10_Bayes_factor" ];then
    sig_cutoff=2
    min=-10

    unstarch $trait_snps | \
    grep -v "#" | \
    awk -F"\t" -v min=$min 'BEGIN{OFS="\t"}{
        if($5=="-Inf"){$5=min};
        print $1,$2,$3,$4,$5
    }' | \
    awk -F"\t" -v sig_cutoff=$sig_cutoff '$5>=sig_cutoff' | \
    target_genes - $all_genes > $targeted_genes  
else
    echo -e "\nERROR: Invalid significance metric. Exiting script."
    exit 1
fi



# Make input to enrichment script
unstarch $all_genes | cut -f 4 | sort | uniq > $all_genes_list

echo -e "distance\tcutoff\ttargeting_method\tsnp_set\tgene_set\tTrait\tTargeted_PositiveControl\tTargeted\tAll_PositiveControl\tAll\tEnrichment" > $outfile

targeted_genes_subset=$TMPDIR/${trait_name}_${pc_gene_list_name}.targeted_genes_subset.txt

distance_sequence=`seq 10000 10000 500000`
for distance in ${distance_sequence[@]};do
    cat $targeted_genes | awk -F"\t" -v cutoff=$distance '$2<=cutoff{print $1}' | sort | uniq > $targeted_genes_subset
    echo -ne "$distance\t${sig_metric}_${sig_cutoff}\tclosest_gene\t$snp_set\t$pc_gene_list_name\t$trait_name\t" >> $outfile
    ./1b_enrichment_calculation.R $all_genes_list $targeted_genes_subset $pc_gene_list >> $outfile
done



