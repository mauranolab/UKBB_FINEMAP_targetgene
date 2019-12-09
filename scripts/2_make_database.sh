#!/bin/bash


src_dir="../data/intermediate_files"
files_array=($(find $src_dir -name "*.txt"))
cat ${files_array[@]} | \
awk -F"\t" '(NR==1 || substr($1,1,8)!="distance") && $NF!=0 && $NF!="NA"' > "../data/all_enrichments.db.txt"
