#!/bin/bash

#set -e -o pipefail

# Define variables
ukbb_dir="/vol/mauranolab/vulpen01/ukbb"
positive_control_genes_dir="/vol/mauranolab/mauram01/ukbb/finemap_multitrait/positive_control_gene_sets"
finemap_data="finemap_data_minBF"
association_data="association_data"
sig_metric=$1
snp_set=$2
targeting_method=$3
if [ -z $4 ];then
    # If no argument is passed to x, assign filler variable.
    x="x"
else
    x=$4
fi

if [ -z $5 ];then
    # If no argument is passed to x, assign filler variable.
    y="y"
else
    y=$5
fi

if [ -z $6 ];then
    # If no argument is passed to x, assign filler variable.
    z="z"
else
    z=$6
fi


if [ $sig_metric == "log10_Bayes_factor" ];then
    if [[ $targeting_method == *"Hi-C"* ]];then
        #snp_file_array=($(find -L $finemap_data -not -path "*trash*" -not -name "*genetic_credible*" \( -name "disease_T2*_finemap.hg19.starch" -o -iname "*glucose*_finemap.hg19.starch" \) ))
        snp_file_array=($(find -L $finemap_data -not -path "*trash*" -not -name "*genetic_credible*" -not -iname "*Morris*" -iname "*GLUCOSE*_finemap.hg19.starch" ))
    else
        snp_file_array=($(find -L $finemap_data -not -path "*trash*" -not -name "*genetic_credible*" -not -iname "*Morris*" -iname "*GLUCOSE*_finemap.hg19.starch"))
        #snp_file_array=($(find -L $finemap_data -not -path "*old_*" -not -name "*genetic_credible*" \( -name "disease_T2*_finemap.hg19.starch" -o -iname "*glucose*_finemap.hg19.starch" \) ))
    fi
    intermediate_files="intermediate_files_V3_2019may14/${targeting_method}/bf_gsea_${snp_set}"


elif [ $sig_metric == "-log10_P-value" ]; then
    if [[ $targeting_method == *"Hi-C"* ]];then
        #snp_file_array=($(find -L $association_data -not -path "*trash*" -not -name "*genetic_credible*" \( -name "disease_T2*.sumstats.hg19.starch" -o -iname "*glucose*.sumstats.hg19.starch" \) ))
        snp_file_array=($(find -L $association_data -not -path "*trash*" -not -name "*genetic_credible*" -not -iname "*Morris*" -iname "*GLUCOSE*.sumstats.hg19.starch" ))
    else
        snp_file_array=($(find -L $association_data -not -path "*trash*" -not -name "*genetic_credible*" -not -iname "*Morris*" -iname "*GLUCOSE*.sumstats.hg19.starch"))
        #snp_file_array=($(find -L $finemap_data -not -path "*old_*" -not -name "*genetic_credible*" \( -name "disease_T2*_finemap.hg19.starch" -o -iname "*glucose*_finemap.hg19.starch" \) ))

    fi
    intermediate_files="intermediate_files_V3_2019may14/${targeting_method}/pval_gsea_${snp_set}"

fi

# echo ${snp_file_array[@]}


# if [ -d $intermediate_files ];then
#     read -p "Delete existing directory? (y/n): " rmOld
#     if [ $rmOld == "y" ] ;then rm -r $intermediate_files;fi
# fi 

log_dir=$intermediate_files/logs
if [ -d $log_dir ]; then rm -r $log_dir;fi
mkdir -p $log_dir


for i in ${!snp_file_array[@]};do
    snp_file=${snp_file_array[$i]}
    infile_base=`basename $snp_file`
    trait_class=${infile_base%%_*}
    trait_name=${infile_base#*_}
    trait_name=${trait_name%_finemap.hg19.starch}
    trait_name=${trait_name%.sumstats.hg19.starch}
    echo $trait_name

    outdir=$intermediate_files/$trait_class
    mkdir -p $outdir


    # Redefine set of SNPs depending on "snp_set" variable
    if [ $snp_set == "DHS_SNPs" ];then
        dhs_file="/vol/mauranolab/mauram01/ukbb/finemap_multitrait/dhs/ENCODE_REMC_FLER_T2D_2019aug.flt.hotspots2.fdr0.05.merged.hg19.starch"
        dhs_snps="$outdir/$infile_base"

        if [ ! -f $dhs_snps ];then
            unstarch $snp_file | grep -v "#" | cut -f 1-5 | \
            bedmap --echo --skip-unmapped - $dhs_file | starch - > $dhs_snps
        fi
        snp_file=$dhs_snps

    elif [ $snp_set == "Trait-Specific_DHS_SNPs" ];then

        dhs_file="/vol/mauranolab/mauram01/ukbb/finemap_multitrait/dhs/${trait_class}_${trait_name}.dhs.hg19.starch"
        if [ ! -e $dhs_file ];then
            echo "No $dhs_file"
            continue
        fi

        dhs_snps="$outdir/$infile_base"
        if [ ! -f $dhs_snps ];then
            unstarch $snp_file | grep -v "#" | cut -f 1-5 | \
            bedmap --echo --skip-unmapped - $dhs_file | starch - > $dhs_snps
        fi
        snp_file=$dhs_snps
    fi


    #if [ $targeting_method == "overlap_gene_body" ] ;then
    #    all_genes="/vol/mauranolab/vulpen01/ukbb/src_data/gencodev24_hg19/lincRNA_protein_coding.gene.bed"
    #else
    all_genes="/vol/mauranolab/vulpen01/ukbb/src_data/tss/lincRNA_ProteinCoding.hg19.starch"
    #fi


    # # Make files containing all genes within targeting distance of any pricelab SNP. These will act as the background to the enrichment calculation
    # all_genes="$intermediate_files/../lincRNA_ProteinCoding.${trait_name}.hg19.starch"
    # if [ ! -f $all_genes ];then
    #     echo "Making background files..."
    #     filename=`basename $snp_file`
    #     filename=${filename%_finemap.*}
    #     filename=${filename%.sumstats*}

    #     source targeting_methods.sh


    #     if [ $targeting_method == "overlap_gene_body" ] ;then
    #         gencode_genes="/vol/mauranolab/vulpen01/ukbb/src_data/gencodev24_hg19/lincRNA_protein_coding.gene.bed"
    #     else
    #         gencode_genes="/vol/mauranolab/vulpen01/ukbb/src_data/tss/lincRNA_ProteinCoding.hg19.starch"
    #     fi

    #     unstarch $association_data/${filename}.sumstats.hg19.starch | \
    #     grep -v "#" | \
    #     target_genes $targeting_method - $gencode_genes | \
    #     cut -f 1-4 | uniq | \
    #     starch - > $all_genes
    # fi

    #echo $positive_control_genes_dir

    #positive_control_genes_file_array=($(find -L $positive_control_genes_dir -type f  -name "*.txt"  -not -path "*src/*" -not -path "*gtex*" ))
    #positive_control_genes_file_array=($(find -L $positive_control_genes_dir -not -path "*src/*" -not -path "*tmp/*" -path "*flannick2019*" -name "*.txt" ))

    # if [ $targeting_method == "pad_10k" ];then
    #     positive_control_genes_file_array=($(find -L $positive_control_genes_dir -name "*.txt" -not -name "*gtex*" -not -path "*src/*" -not -path "*drug_genes*" -not -path "*disease_genes*" -not -path "*pancreas*" -not -name "*ibd-disease*" -not -name "*afib-disease*" -not -name "*ra-disease*" -not -name "*Tissue_Non-Sp*" -not -path "*bak.2019summer*" -not -path "*trash*"))
    # elif [[ $targeting_method == *"Hi-C"* ]];then
    #     positive_control_genes_file_array=($(find -L $positive_control_genes_dir -name "*_all_pc.txt" -path "*all_pc_combined*" -not -path "*bak.2019summer*" -not -path "*trash*" -not -path "*src/*"))
    # elif [[ $snp_file == *"T2D"* ]];then
    #     positive_control_genes_file_array=($(find -L $positive_control_genes_dir \( -iname "t2d*.txt" -o -iname "*diabetes*.txt" \) -not -path "*src/*" -not -name "*genes.txt" -not -path "*bak.2019summer*" -not -path "*trash*"))
    #     #positive_control_genes_file_array=($(find -L $positive_control_genes_dir -iname "*diabetes*.txt" -not -path "*src/*" -not -name "*genes.txt")) 
    # else
    #     positive_control_genes_file_array=($(find -L $positive_control_genes_dir -name "*_all_pc.txt"  -not -name "*gtex*" -not -path "*gtex*" -not -path "*encode*" -not -path "*src/*" -not -path "*drug_genes*" -not -path "*disease_genes*" -not -path "*pancreas*" -not -name "*ibd-disease*" -not -name "*afib-disease*" -not -name "*ra-disease*" -not -name "*Tissue_Non-Sp*" -not -path "*bak.2019summer*" -not -path "*trash*"))
    # fi

    #echo ${positive_control_genes_file_array[@]}

    grep $trait_name "trait_genes_key.txt" | \
    awk -F"\t" 'BEGIN{OFS="\n"}{$1=NA;print $0}' > $TMPDIR/${trait_name}.pc_gene_sets.txt


    while read -r pc_set;do


        if [ "$pc_set" == "NA" ] || [ "$pc_set" == "" ];then
            continue
        fi

        pc_set_file=$(find -L $positive_control_genes_dir -name "*${pc_set}.txt" -name "*.txt" -not -name "*string*" -not -path "*src/*" -not -path "*bak.*" -not -path "*trash*")

        if [ "$pc_set_file" == "" ];then
            continue
        fi


        if [[ $targeting_method == *"Hi-C"* ]];then
            hic_dir="../hic_data"
            hic_file_array=($(find -L $hic_dir -iname "*${trait_name}*starch" ))


            for hic in ${hic_file_array[@]};do
                paper=`basename $hic`
                paper=${paper%.PCHi-C.*}
                paper=${paper%.Hi-C.*}
                paper=${paper#*.}
                targeting_method_paper=${targeting_method}_${paper}
                hic_outdir=${outdir}/${trait_name}_${paper}

                mkdir -p $hic_outdir
                outfile_txt=$hic_outdir/${trait_name}_${pc_set}.txt

                echo $outfile_txt

                qsub -S /bin/bash -N ${trait_name}_${pc_set} -o $log_dir "./runpair.enrichment.sh $snp_file $targeting_method_paper $pc_set_file $all_genes $outfile_txt $sig_metric $snp_set $hic"

            done
        elif [[ $targeting_method == *"eQTL"* ]] ;then
            # eqtl_dir="../eQTL_data"
            # eqtl_file_array=($(find -L $eqtl_dir -iname "*${trait_name}*.starch" ))

            # for eqtl in ${eqtl_file_array[@]};do
            #      outfile_txt=$outdir/${trait_name}_${pc_set}.txt
            #     echo $outfile_txt
            #     qsub -S /bin/bash -N ${trait_name}_${pc_set} -o $log_dir "./runpair.enrichment.sh $snp_file $targeting_method $pc_set_file $all_genes $outfile_txt $sig_metric $snp_set $eqtl"
            # done


            grep $trait_name "trait_eqtl_key.txt" | \
            awk -F"\t" 'BEGIN{OFS="\n"}{$1=$1;print $0}' > $TMPDIR/${trait_name}.pc_eqtl_sets.txt

            while read -r eqtl_set;do
                if [ "$eqtl_set" == "" ];then
                    continue
                fi

                eqtl_dir="../eQTL_data"
                eqtl_file_array=($(find -L $eqtl_dir -iname "*${eqtl_set}*.starch" ))
                eqtl_outdir=$outdir/${trait_name}
                mkdir -p $eqtl_outdir

                for eqtl in ${eqtl_file_array[@]};do
                    outfile_txt=$eqtl_outdir/${trait_name}_${pc_set}_${eqtl_set}.txt
                    echo $outfile_txt
                    qsub -S /bin/bash -N ${trait_name}_${pc_set} -o $log_dir "./runpair.enrichment.sh $snp_file ${targeting_method}_${eqtl_set} $pc_set_file $all_genes $outfile_txt $sig_metric $snp_set $eqtl"
                done

            done < $TMPDIR/${trait_name}.pc_eqtl_sets.txt

        else
            outfile_txt=$outdir/${trait_name}_${pc_set}.txt
            echo $outfile_txt
            qsub -S /bin/bash -N ${trait_name}_${pc_set} -o $log_dir "./runpair.enrichment.sh $snp_file $targeting_method $pc_set_file $all_genes $outfile_txt $sig_metric $snp_set $x $y $z"
        fi

        





    done < $TMPDIR/${trait_name}.pc_gene_sets.txt



    # if [[ $targeting_method == *"Hi-C"* ]];then
    #     for pc_set_file in ${positive_control_genes_file_array[@]};do
    #         echo $pc_set_file
    #         pc_set=`basename $pc_set_file`
    #         pc_set=${pc_set%.txt}
            

    #         hic_dir="../hic_data"
    #         hic_file_array=($(find -L $hic_dir -name "*${trait_name}*.starch" ))


    #         for hic in ${hic_file_array[@]};do
    #             paper=`basename $hic`
    #             paper=${paper%.PCHi-C.*}
    #             paper=${paper%.Hi-C.*}
    #             paper=${paper#*.}
    #             targeting_method_paper=${targeting_method}_${paper}
    #             hic_outdir=${outdir}/${trait_name}_${paper}

    #             mkdir -p $hic_outdir
    #             outfile_txt=$hic_outdir/${trait_name}_${pc_set}.txt

    #             echo $outfile_txt

    #             qsub -S /bin/bash -N ${trait_name}_${pc_set} -o $log_dir "./runpair.enrichment.sh $snp_file $targeting_method_paper $pc_set_file $all_genes $outfile_txt $sig_metric $snp_set $hic"

    #         done
    #     done
    # else
    #     for pc_set_file in ${positive_control_genes_file_array[@]};do
    #         echo $pc_set_file
    #         pc_set=`basename $pc_set_file`
    #         pc_set=${pc_set%.txt}
    #         outfile_txt=$outdir/${trait_name}_${pc_set}.txt


    #         echo $outfile_txt

    #         qsub -S /bin/bash -N ${trait_name}_${pc_set} -o $log_dir "./runpair.enrichment.sh $snp_file $targeting_method $pc_set_file $all_genes $outfile_txt $sig_metric $snp_set $x $y $z"
    #     done
    # fi

done


