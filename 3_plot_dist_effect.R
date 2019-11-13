#!/usr/bin/env Rscript

library(ggplot2)
theme_set(theme_classic())
library(ggthemes)
library(data.table)
library(dplyr)

outdir<-"summary_plots"
system(paste0("mkdir -p ", outdir))

process_infile<-function(all_enrichments, pc_genes_key){
	infile_df<-fread(all_enrichments)
	infile_df$snp_set[which(infile_df$snp_set=="All_SNPs")]<-"All SNPs"
	infile_df$snp_set[which(infile_df$snp_set=="DHS_SNPs")]<-"DHS SNPs"
	infile_df$snp_set[which(infile_df$snp_set=="Trait-Specific_DHS_SNPs")]<-"DHS SNPs\n(Trait-Specific)"
	
	infile_df<-infile_df[which(infile_df$distance>=10000),]
	infile_df<-infile_df[which(infile_df$Targeted_PositiveControl>0),]
	
	#write.table(infile_df[which(distance==1000000),1:7],col.names = T, row.names = F, quote = F, sep = "\t")
	
	# Exclude pval DHS and pval trait specific DHS
	infile_df$cutoff[which(infile_df$cutoff=="-log10_P-value_7.30103")]<-"-log10_P-value_5e-8"
	infile_df<-infile_df[which(cutoff=="log10_Bayes_factor_2" | (cutoff=="-log10_P-value_5e-8" & snp_set=="All SNPs")),]
	infile_df<-infile_df[grepl("closest_gene",infile_df$targeting_method),]
	infile_df<-infile_df[!grepl("_all_pc",annotation),]
	infile_df<-infile_df[!grepl("string",annotation),]
	
	#select_enrichments<-fread("trait_genes_key.compare_gene_sets.txt")
	select_enrichments<-fread(pc_genes_key)
	select_enrichments_melt<-data.table::melt(select_enrichments, id.vars="Trait")
	# select_enrichments_melt<-data.table::melt(select_enrichments[,c("Trait",
	# 																				 "Expression (GTEx)",
	# 																				 "Expression (scRNA-seq)",
	# 																				 "Disease + Drug (combined)",
	# 																				 "Mouse diabetes (Flannick 2019)",
	# 																				 "String")], 
	# 																					id.vars="Trait")
	colnames(select_enrichments_melt)<-c("trait","gene_set","annotation")
	
	select_enrichments_melt<-select_enrichments_melt[!which(trait%in%c("AFIB","RA","IBD")),]
	
	select_enrichments_melt$gene_set<-as.character(select_enrichments_melt$gene_set)
	select_enrichments_melt$gene_set<-sub("\\(","\n\\(",select_enrichments_melt$gene_set) 
	infile_df<-merge(infile_df,select_enrichments_melt, by=c("trait","annotation"),
									 all=F)
	#write.table(infile_df,file = "tmp.txt", quote = F, sep = "\t",row.names = F, col.names = T)
	
	
	infile_df$trait<-gsub("HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT","HIGH_LIGHT_SCATTER\nRETICULOCYTE_COUNT",infile_df$trait)
	infile_df$trait<-gsub("MEAN_CORPUSCULAR_HEMOGLOBIN","MEAN_CORPUSCULAR\nHEMOGLOBIN",infile_df$trait)
	infile_df$trait<-gsub("AFIB","Atrial\nFibrillation",infile_df$trait)
	infile_df$trait<-gsub("^RA$","Rheumatoid\nArthritis",infile_df$trait)
	infile_df$trait<-gsub("T2Dbmiadj","Type 2 Diabetes",infile_df$trait)
	infile_df$trait<-gsub("DIASTOLIC","Blood Pressure\n(Diastolic)",infile_df$trait)
	infile_df$trait<-gsub("SYSTOLIC","Blood Pressure\n(Systolic)",infile_df$trait)
	infile_df$trait<-gsub("HEIGHTz","Height",infile_df$trait)
	infile_df$trait<-gsub("LOW_TSH","Hypothyroidism",infile_df$trait)
	infile_df$trait<-gsub("RED_COUNT","RBC Count",infile_df$trait)
	infile_df$trait<-gsub("LDL_CHOLESTEROL","LDL\nCholesterol",infile_df$trait)
	infile_df$trait<-gsub("TRIGLYCERIDES","Triglycerides",infile_df$trait)
	infile_df$trait<-gsub("DBILIRUBIN","Bilirubin (D)",infile_df$trait)
	infile_df$trait<-gsub("TBILIRUBIN","Bilirubin (T)",infile_df$trait)
	infile_df$trait<-gsub("GLUCOSE","Glucose",infile_df$trait)
	infile_df$trait<-gsub("CALCIUM","Calcium",infile_df$trait)
	
	infile_df$cutoff<-sub("log10_Bayes_factor_","log10 BF > ",infile_df$cutoff)
	infile_df$cutoff<-sub("-log10_P-value_","P-value < ",infile_df$cutoff)
	
	infile_df$col_facet<-paste(infile_df$snp_set,infile_df$cutoff, sep = "\n")
	
	facet_order<-c("All SNPs\nP-value < 5e-8",
								 "All SNPs\nlog10 BF > 2",
								 "DHS SNPs\nlog10 BF > 2",
								 "DHS SNPs\n(Trait-Specific)\nlog10 BF > 2"
	)
	infile_df$col_facet<-factor(infile_df$col_facet,levels = facet_order, ordered = T)

	
	return(infile_df)
	
}



################################
# Absolute distance effect

plot_df<-process_infile("all_enrichments.db.txt","trait_genes_key.selected.txt")


gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

color_vec<-gg_color_hue(11)

colors<-rep(NA, times=length(plot_df$trait))
colors[which(plot_df$trait=="Bilirubin (D)")]<-color_vec[1]
colors[which(plot_df$trait=="Bilirubin (T)")]<-color_vec[1]
colors[which(plot_df$trait=="Calcium")]<-color_vec[2]
colors[which(plot_df$trait=="eBMD")]<-color_vec[3]
colors[which(plot_df$trait=="Blood Pressure\n(Systolic)")]<-color_vec[4]
colors[which(plot_df$trait=="Blood Pressure\n(Diastolic)")]<-color_vec[4]
colors[which(plot_df$trait=="Height")]<-color_vec[5]
colors[which(plot_df$trait=="Hypothyroidism")]<-color_vec[6]
colors[which(plot_df$trait=="RBC Count")]<-color_vec[7]
colors[which(plot_df$trait=="Type 2 Diabetes")]<-color_vec[8]
colors[which(plot_df$trait=="LDL\nCholesterol")]<-color_vec[9]
colors[which(plot_df$trait=="Triglycerides")]<-color_vec[10]
colors[which(plot_df$trait=="Glucose")]<-color_vec[11]

names(colors)<-plot_df$trait




p<-ggplot(plot_df, aes(x = distance/1000, y = enrichment, color = trait)) +
	geom_line() +
	geom_hline(yintercept = 1, color="grey", linetype = "dashed") +
	theme(
		text = element_text(size = 12),
		strip.background = element_blank(),
		# plot.margin = margin(c(0,0,0,0), unit = "mm"),
		# legend.margin = margin(c(0,0,0,0), unit = "mm"),
		legend.position = "right",
		legend.title = element_blank(),
		legend.text = element_text(size = 9, margin = margin(c(0,0.25,0,0.25), unit = "in")),
		strip.text.y = element_text(size = 10),
		axis.text.x = element_text(angle = 45, hjust = 1),
		aspect.ratio = 1
	) +
	facet_grid(cols = vars(col_facet), rows = vars(gene_set), scales = "free_y") +
	labs(x="Distance (Kbp)", y = "Fold enrichment") +
	scale_color_manual(name = "Trait", values = colors, aesthetics = c("colour"), guide=guide_legend(ncol=1)) +
	coord_cartesian(xlim = c(0,250))

pdf(paste0(outdir,"/absolute_distance_effect.minBF.pdf"), height = 4, width = 8)
p
dev.off()

print_df<-plot_df
#print_df$snp_set[which(print_df$snp_set=="DHS SNPs\n(Trait-Specific)")]<-"DHS SNPs (Trait-Specific)"
print_df$snp_set<-gsub("\n"," ",print_df$snp_set)
print_df$gene_set<-gsub("\n"," ",print_df$gene_set)
print_df$trait<-gsub("\n"," ",print_df$trait)
print_df$enrichment<-round(print_df$enrichment,digits = 2)
#write.table(print_df[,c("trait","distance","cutoff","snp_set" ,"enrichment","gene_set")], col.names = T, row.names = F, quote = F, file = "",sep = "\t")
write.table(print_df[,c("trait","distance","cutoff","snp_set","Targeted_PositiveControl" ,"enrichment","gene_set")], col.names = T, row.names = F, quote = F, file = "absolute_distance_effect.minBF.txt",sep = "\t")



################################
# Relative distance effect

# Ordering plot
trait_order<-c("Type 2 Diabetes\n(BMI adj.)",
							 "Type 2 Diabetes",
							 "Glucose",
							 "Atrial\nFibrillation",
							 "THYROID ANY SELF REP",
							 "AID ALL",
							 "RA",
							 "Rheumatoid\nArthritis",
							 "IBD",
							 "Blood Pressure\n(Systolic)",
							 "SYSTOLIC",
							 "Blood Pressure\n(Diastolic)",
							 "DIASTOLIC",
							 "LDL CHOLESTEROL",
							 "LDL\nCholesterol",
							 "TRIGLYCERIDES",
							 "Triglycerides",
							 "LYMPHOCYTE COUNT",
							 "EOSINOPHIL COUNT",
							 "MONOCYTE COUNT",
							 "RED COUNT",
							 "RBC Count",
							 "PLATELET COUNT",
							 "MEAN CORPUSCULAR\nHEMOGLOBIN",
							 "eBMD",
							 "HEIGHTz",
							 "Height",
							 "LOW_TSH",
							 "Hypothyroidism\n(Low TSH)",
							 "Hypothyroidism",
							 "Bilirubin (D)",
							 "Bilirubin (T)",
							 "Calcium")

plot_df_selected<-plot_df

plot_df_all<-process_infile("all_enrichments.db.txt","trait_genes_key.txt")

plot_df_unselected<-setdiff(plot_df_all,plot_df_selected)

plot_df_unselected$selected<-"unselected"
plot_df_selected$selected<-"selected"

plot_df<-bind_rows(plot_df_unselected,plot_df_selected)

plot_df$trait<-factor(plot_df$trait,levels = trait_order, ordered = T)

max_p_df <- plot_df %>%
	filter(cutoff == "P-value < 5e-8") %>%
	group_by_(.dots = c("trait","annotation")) %>%
	summarize(max_p = max(enrichment)) %>%
	select(trait,annotation,max_p)

plot_df<-right_join(max_p_df,plot_df)

#head(plot_df)

relative_df <- plot_df %>%
	group_by_(.dots = c("trait","annotation","cutoff","snp_set")) %>%
	mutate(norm_enrichment = enrichment/max_p)

#######
# Bayes factor
#######

#plot_df<-relative_df[which(grepl("Bayes",relative_df$cutoff)),]

write.table(plot_df[which(plot_df$distance==100000),], file = "", sep = "\t", quote = F, row.names = F)

p<-ggplot(relative_df, aes(x = distance/1000, y = norm_enrichment, color = gene_set)) +
	geom_hline(yintercept = 1, color="grey", linetype = "dashed") +
	geom_line(aes(linetype=selected)) +
	scale_colour_manual(values = c("Disease + Drug"="#0072B2","Expression"="#E69F00"),guide=guide_legend(nrow=1)) +
	theme(
		text = element_text(size = 12),
		strip.background = element_blank(),
		#legend.box.margin = margin(c(0,0,0,0), unit = "mm"),
		legend.text = element_text(size = 9, margin = margin(c(0,0.25,0,0.25), unit = "in")),
		legend.position = "top",
		legend.title = element_blank(),
		strip.text.y = element_text(size = 8),
		axis.text.x = element_text(angle = 45, hjust = 1),
		aspect.ratio = 1
	) +
	facet_grid(cols = vars(col_facet), rows = vars(trait), scales = "free_y") +
	labs(x="Distance (Kbp)", y = "Relative enrichment") +
	coord_cartesian(xlim = c(0,250)) +
	scale_y_continuous(expand=c(0,0)) +
	guides(linetype=FALSE)

pdf(paste0(outdir,"/relative_distance_effect.minBF.pdf"), height = 16, width = 6)
p
dev.off()
