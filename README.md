# Enrichment of positive control genes by distance from SNP

## Running this script
**Navigate to scripts/ and run 1a_run_enrichments.sh with the following arguments (ordered):**
1. Path to SNP bed/starch file with p-val or log10 bf in column 5. File should already be filtered to remove MHC and Missense SNPs.
2. Trait name. Corresponds to column 1 in trait_genes_key.txt. (Must match exactly)
3. Either "-log10_P-value" or "log10_Bayes_factor". Describes column 5 of "trait_snps". 
4. Either "All_SNPs", "DHS_SNPs", or "Trait-Specific_DHS_SNPs". Describes first argument. 
5. Path to txt file of positive control genes. Format: one column of gene IDs
6. Name of positive control set. Corresponds to columns 2&3 in trait_genes_key.txt. (Must match exactly)
7. Path to bed/starch file with gene ID in column 4. Contains all genes in transcript model
8. Directory name for output files.

Example:
./1a_run_enrichments.sh "$snp_dir/glucose_dhs_snps_pvals.bed" \
	"GLUCOSE" \
	"-log10_P-value" \
	"DHS_SNPs" \
	"$positive_control_gene_dir/glucose_expression_genes.txt" \
	"Glucose_expression" \
	"$transcript_models/gencodev24_genes.txt" \
	"../data/intermediate_files/closest_gene/pval_gsea_DHS_SNPs/"

Run this script for each trait/positive control gene set pairing.
Name the output directory to separate files using different SNP sets (e.g. All_SNPs vs DHS_SNPs). This will avoid collision in output files.



**Run ./2_make_database.sh**
Combines all files within "../data/intermediate_files/" (the output of ./1a_run_enrichments.sh).
Produces file "../data/all_enrichments.db.txt"

**Run 3_plot_dist_effect.R**
Reads "../data/all_enrichments.db.txt" to produce plot at "../plots/distance_effect.pdf"

