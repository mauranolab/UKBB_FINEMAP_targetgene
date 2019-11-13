#!/bin/bash

#set -e -o pipefail

#Output format:
#Columns 1-4: a subset of the "genes" input.
#Column 5: the maximum value from column 5 of the "snps" input that targeted the gene.


#BUGBUG cannot yet handle infinite values

function target_genes {
    target_method=$1
    snps=$2
    genes=$3
    x=$4 # This variable is for inputs that are not applicable to all targeting methods
    y=$5 # This variable is for inputs that are not applicable to all targeting methods
    z=$6 # This variable is for inputs that are not applicable to all targeting methods

    ctcf="/vol/isg/encode/CTCF_maurano_cell_reports_2015/ctcf.all.hg19.starch"
    hg19_coords="/vol/isg/annotation/bed/hg19/mapping/chromInfo.bed"
    all_dhs="/vol/mauranolab/mauram01/gwas_enrichmentplots/hg19/ENCODE_REMC_FLER_2018may.flt.hotspots2.fdr0.05.hg19.starch"
    coding_snps="/vol/mauranolab/vulpen01/isg_annotation/vep/hg19/gencodev19_CDS/Gencodev19.CDS_vep.hg19.starch"



    # All functions must either return all distance values, or return the smallest ABSOLUTE distance value
    if [ $target_method == "closest_gene" ];then

        # Returns the distance between any gene and all significant SNPs
        gcat $snps | grep -v "#" | \
        closest-features --closest --delim "\t" --dist $genes - | \
        awk -F"\t" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$NF}' | \
        sort-bed - | uniq
        
    elif [[ $target_method == *"eQTL"* ]];then
        
        eqtl=$x

        gcat $snps | grep -v "#" | \
        bedmap --echo --skip-unmapped $eqtl - | \
        sort -nk5 | awk -F"\t" '!seen[$4]++' | \
        sort-bed - | uniq

    elif [[ $target_method == *"Hi-C"* ]];then

        #padding_range=$y
        #min_chicago_score=$z


        rstring=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c 15 ; echo ''`
        Hi_C=$TMPDIR/Hi-C_data.${rstring}.bed


        gcat $x | grep -v "#" | \
        # awk -F"\t" -v min_chicago_score=${min_chicago_score} 'BEGIN{OFS="\t"}$7>min_chicago_score{print $1,$2,$3,$4,$5,$6}' | \
        # awk -F"\t" -v padding_range=$padding_range ' # Finds midpoint and pads +/- by padding range
        # function abs(v) {return v < 0 ? -v : v}
        # BEGIN{OFS="\t"}{
        #     l_midpoint=int(($2+$3)/2)
        #     r_midpoint=int(($5+$6)/2)
            
        #     print $1,l_midpoint-padding_range,l_midpoint+1+padding_range,$4,r_midpoint-padding_range,r_midpoint+1+padding_range
        # }' | \
        # awk -F"\t" 'BEGIN{OFS="\t"}{if($2<0){$2=0}if($5<0){$5=0};print $0}'| \
        cut -f 1-6 | \
        awk -F"\t" 'BEGIN{OFS="\t"}{
            print $1,$2,$3,$4,$5,$6
            print $4,$5,$6,$1,$2,$3
        }' | \
        sort-bed - | uniq > $Hi_C

        #head $Hi_C >&2

        gcat $snps | grep -v "#" | \
        bedmap --echo --echo-map --skip-unmapped $Hi_C - | \
        awk -F"|" 'BEGIN{OFS="\t"}{split($1,a,"\t")}{print a[4],a[5],a[6],".|"$2}' | sort-bed - | \
        bedmap --echo --echo-map --skip-unmapped $genes -  | \
        awk -F"|" 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{
            split($1,gene,"\t");split($3,sig_snp_array,";")
            for(n in sig_snp_array){
                split(sig_snp_array[n],sig_snp,"\t")
                print gene[1],gene[2],gene[3],gene[4],abs(gene[2]-sig_snp[2])}
            }' | sort-bed - | uniq


    fi


}




