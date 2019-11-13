#!/bin/bash

set -eu -o pipefail


################################################################################
#                        BAYES FACTOR SNPs                                     #
################################################################################
sig_metric="log10_Bayes_factor"

targeting_method="closest_gene"
./MASTER.1.sh $sig_metric "All_SNPs" $targeting_method
./MASTER.1.sh $sig_metric "DHS_SNPs" $targeting_method
./MASTER.1.sh $sig_metric "Trait-Specific_DHS_SNPs" $targeting_method


#################### Targeting method: Hi-C ########################
hi_c_dir="../hic_data"

min_chicago_score=0
padding_distance=0
targeting_method="Hi-C"

./MASTER.1.sh $sig_metric "All_SNPs" "PC${targeting_method}"
./MASTER.1.sh $sig_metric "DHS_SNPs" "PC${targeting_method}"
./MASTER.1.sh $sig_metric "Trait-Specific_DHS_SNPs" "PC${targeting_method}"

#################### Targeting method: eQTL ########################
hi_c_dir="../eQTL_data"

targeting_method="eQTL"

./MASTER.1.sh $sig_metric "All_SNPs" "${targeting_method}"
./MASTER.1.sh $sig_metric "DHS_SNPs" "${targeting_method}"
./MASTER.1.sh $sig_metric "Trait-Specific_DHS_SNPs" "${targeting_method}"

################################################################################
#                            P-VALUE SNPs                                      #
################################################################################
sig_metric="-log10_P-value"

targeting_method="closest_gene"
./MASTER.1.sh $sig_metric "All_SNPs" $targeting_method


#################### Targeting method: Hi-C ########################
hi_c_dir="../hic_data"

min_chicago_score=5
padding_distance=0
targeting_method="Hi-C"

./MASTER.1.sh $sig_metric "All_SNPs" "PC${targeting_method}"


#################### Targeting method: Hi-C ########################
hi_c_dir="../eQTL_data"

targeting_method="eQTL"

./MASTER.1.sh $sig_metric "All_SNPs" "${targeting_method}"


