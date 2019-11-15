#!/bin/bash


src_dir="intermediate_files"

files_array=($(find $src_dir -name "*.txt" -not -path "*trash*"))

cat ${files_array[@]} | \
awk -F"\t" '(NR==1 || substr($1,1,8)!="distance") ' | grep -iv "string" > all_enrichments.db.txt
#awk -F"\t" '(NR==1 || substr($1,1,8)!="distance") && $NF!=0 && $NF!="NA"' > all_enrichments.db.txt
