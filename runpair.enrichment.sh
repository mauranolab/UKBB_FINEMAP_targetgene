#!/bin/bash

#set -e -o pipefail

snp_file=$1 # starch file
targeting_method=$2
positive_control_gene_list=$3 # txt file
all_genes=$4 # starch file
outfile_txt=$5
sig_metric=$6
snp_set=$7
x=$8
y=$9
z=${10}

MHC_coords_hg19="/vol/mauranolab/vulpen01/ukbb/src_data/gencodev24_hg19/MHC_coordinates.hg19.bed"
MHC_genes_gencodev24="/vol/mauranolab/vulpen01/ukbb/src_data/positive_control_genes/src/MHC_genes.txt"

source targeting_methods.sh


trait_name=`basename $snp_file`
trait_name=${trait_name%.sumstats.hg*.starch}
trait_name=${trait_name%_finemap.hg*.starch}
trait_class=${trait_name%%_*}
trait_name=${trait_name#*_}

annotation_name=`basename $positive_control_gene_list`
annotation_name=${annotation_name%.txt}


all_genes_list=$TMPDIR/${trait_name}_${annotation_name}.all_genes_list.txt # 1 column: Gene name (tab delimited)
#positive_control_genes_list=$TMPDIR/${trait_name}_${annotation_name}.positive_control_genes_list.txt # 1 column: Gene name (tab delimited)
targeted_genes=$TMPDIR/${trait_name}_${annotation_name}.targeted_genes.txt # 2 columns: Gene name, value (tab delimited)
#targeted_genes=${trait_name}_${annotation_name}.targeted_genes.txt # 2 columns: Gene name, value (tab delimited)
positive_control_genes_noMHC=$TMPDIR/${trait_name}_${annotation_name}.positive_control_genes.txt


#comm -23 $positive_control_gene_list $MHC_genes_gencodev24 > $positive_control_genes_noMHC
#grep -vf $MHC_genes_gencodev24 $positive_control_gene_list > $positive_control_genes_noMHC


#echo $targeted_genes >> "/dev/stderr"

#echo $targeted_genes


if [ $sig_metric == "-log10_P-value" ];then
    #sig_cutoff=5e-8
    sig_cutoff=7.30103
    min=0
    max=30
    # increment=2

    # sequence=`seq $min $increment $max`
    # sequence+=(7.30103) # 7.30103 == -log10(5e-8)
    # sequence+=(3.30103) # 3.30103 == -log10(5e-4)
fi


if [ $sig_metric == "log10_Bayes_factor" ];then
    sig_cutoff=2

    min=-10
    max=14
    # increment=0.5

    # sequence=`seq $min $increment $max`
fi



# Pair each gene with the minimum distance to a significant SNP
if [ $sig_metric == "-log10_P-value" ];then

    unstarch $snp_file | grep -v "#" | \
    bedops -n - $MHC_coords_hg19 | \
    awk -F"\t" -v max=$max 'BEGIN{OFS="\t"}{log10p=-log($5)/log(10); if(log10p!="inf"){print $1,$2,$3,$4, log10p}else{print $1,$2,$3,$4,max}}' | \
    awk -F"\t" -v sig_cutoff=$sig_cutoff '$5>=sig_cutoff' | \
    #closest-features --closest --delim "\t" --dist $all_genes - | \
    target_genes $targeting_method - $all_genes $x $y $z | \
    awk -F"\t" 'function abs(v) {return v < 0 ? -v : v} BEGIN{OFS="\t"}{if($NF!="NA"){print $4,abs($NF)}}' | \
    sort -nk2 | awk -F"\t" '!seen[$1]++' > $targeted_genes 

    #cat $targeted_genes >> "/dev/stderr"

    
    # sort -k1,1 | \
    # awk -F"\t" 'BEGIN{OFS="\t";gene_name="blank"}{
    #     if(gene_name!=$1){if(NR>1){print gene_name,val}; gene_name=$1;val=$2}
    #     else if(gene_name==$1 && val<$2){gene_name=$1;val=$2}}
    #     END{print gene_name,val}' > $targeted_genes 

else

    unstarch $snp_file | grep -v "#" | \
    bedops -n - $MHC_coords_hg19 | \
    awk -F"\t" -v min=$min 'BEGIN{OFS="\t"}{if($5=="-Inf"){$5=min};print $1,$2,$3,$4,$5}' | \
    awk -F"\t" -v sig_cutoff=$sig_cutoff '$5>=sig_cutoff' | \
    #closest-features --closest --delim "\t" --dist $all_genes - | \
    target_genes $targeting_method - $all_genes $x $y $z | \
    awk -F"\t" 'function abs(v) {return v < 0 ? -v : v} BEGIN{OFS="\t"}{if($NF!="NA"){print $4,abs($NF)}}' | \
    sort -nk2 | awk -F"\t" '!seen[$1]++' > $targeted_genes  

    #cat $targeted_genes >> "/dev/stderr"

    # sort -k1,1 | \
    # awk -F"\t" 'BEGIN{OFS="\t";gene_name="blank"}{
    #     if(gene_name!=$1){if(NR>1){print gene_name,val}; gene_name=$1;val=$2}
    #     else if(gene_name==$1 && val<$2){gene_name=$1;val=$2}}
    #     END{print gene_name,val}' > $targeted_genes  
fi




# Make input to enrichment script
unstarch $all_genes | cut -f 4 | sort | uniq > $all_genes_list
#unstarch $positive_control_genes | cut -f 4 | sort | uniq > $positive_control_genes_list

if [ -f $outfile_txt ];then rm $outfile_txt;fi
echo -e "distance\tcutoff\ttargeting_method\tsnp_set\tannotation\ttrait_class\ttrait\tTargeted_PositiveControl\tTargeted\tAll_PositiveControl\tAll\tenrichment" > $outfile_txt


targeted_genes_subset=$TMPDIR/${trait_name}_${annotation_name}.targeted_genes_subset.txt

sequence=`seq 10000 10000 1000000`



for value in ${sequence[@]};do
    cat $targeted_genes | awk -F"\t" -v cutoff=$value '$2<=cutoff{print $1}' | sort | uniq > $targeted_genes_subset

    echo -ne "$value\t${sig_metric}_${sig_cutoff}\t$targeting_method\t$snp_set\t$annotation_name\t$trait_class\t$trait_name\t" >> $outfile_txt

    ../gsea.enrichment.R $all_genes_list $targeted_genes_subset $positive_control_gene_list >> $outfile_txt

done



