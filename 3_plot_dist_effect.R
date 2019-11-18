#!/usr/bin/env Rscript

library(plotly)
theme_set(theme_classic())
library(data.table)
library(dplyr)

outdir<-"/vol/mauranolab/vulpen01/ukbb/github_upload/enrichment_by_distance/summary_plots"
system(paste0("mkdir -p ", outdir))

process_infile<-function(all_enrichments, pc_genes_key){
	infile_df<-fread(all_enrichments)
	infile_df$snp_set[which(infile_df$snp_set=="All_SNPs")]<-"All SNPs"
	infile_df$snp_set[which(infile_df$snp_set=="DHS_SNPs")]<-"DHS SNPs"
	infile_df$snp_set[which(infile_df$snp_set=="Trait-Specific_DHS_SNPs")]<-"DHS SNPs\n(Trait-Specific)"
	
	infile_df<-infile_df[which(infile_df$Targeted_PositiveControl>0),]
	
	
	# Exclude pval DHS and pval trait specific DHS
	infile_df$cutoff[which(infile_df$cutoff=="-log10_P-value_7.30103")]<-"-log10_P-value_5e-8"
	infile_df<-infile_df[grepl("closest_gene",infile_df$targeting_method),]
	infile_df<-infile_df[!grepl("_all_pc",gene_set),]

	
	select_enrichments<-fread(pc_genes_key)
	select_enrichments_melt<-data.table::melt(select_enrichments, id.vars="Trait")
	colnames(select_enrichments_melt)<-c("Trait","annotation","gene_set")
	
	select_enrichments_melt$gene_set<-as.character(select_enrichments_melt$gene_set)
	select_enrichments_melt$gene_set<-sub("\\(","\n\\(",select_enrichments_melt$gene_set) 
	
	
	infile_df<-merge(infile_df,select_enrichments_melt, by=c("Trait","gene_set"),all=F)


	infile_df$Trait<-gsub("T2Dbmiadj","Type 2 Diabetes",infile_df$Trait)
	infile_df$Trait<-gsub("DIASTOLIC","Blood Pressure\n(Diastolic)",infile_df$Trait)
	infile_df$Trait<-gsub("SYSTOLIC","Blood Pressure\n(Systolic)",infile_df$Trait)
	infile_df$Trait<-gsub("HEIGHTz","Height",infile_df$Trait)
	infile_df$Trait<-gsub("LOW_TSH","Hypothyroidism",infile_df$Trait)
	infile_df$Trait<-gsub("RED_COUNT","RBC Count",infile_df$Trait)
	infile_df$Trait<-gsub("LDL_CHOLESTEROL","LDL Cholesterol",infile_df$Trait)
	infile_df$Trait<-gsub("TRIGLYCERIDES","Triglycerides",infile_df$Trait)
	infile_df$Trait<-gsub("DBILIRUBIN","Bilirubin (D)",infile_df$Trait)
	infile_df$Trait<-gsub("GLUCOSE","Glucose",infile_df$Trait)
	infile_df$Trait<-gsub("CALCIUM","Calcium",infile_df$Trait)

	infile_df$cutoff<-sub("log10_Bayes_factor_","log10 BF > ",infile_df$cutoff)
	infile_df$cutoff<-sub("-log10_P-value_","P-value < ",infile_df$cutoff)

	infile_df$col_facet<-paste(infile_df$snp_set,infile_df$cutoff, sep = "\n")

	facet_order<-c("All SNPs\nP-value < 5e-8",
								 "All SNPs\nlog10 BF > 2",
								 "DHS SNPs\nlog10 BF > 2",
								 "DHS SNPs\n(Trait-Specific)\nlog10 BF > 2")
	infile_df$col_facet<-factor(infile_df$col_facet,levels = facet_order, ordered = T)


	return(infile_df)

}



################################
# Absolute distance effect

plot_df<-process_infile("all_enrichments.db.txt","trait_genes_key.selected.txt")

head(plot_df)

gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

color_vec<-gg_color_hue(11)

colors<-rep(NA, times=length(plot_df$Trait))
colors[which(plot_df$Trait=="Bilirubin (D)")]<-color_vec[1]
colors[which(plot_df$Trait=="Calcium")]<-color_vec[2]
colors[which(plot_df$Trait=="eBMD")]<-color_vec[3]
colors[which(plot_df$Trait=="Blood Pressure\n(Systolic)")]<-color_vec[4]
colors[which(plot_df$Trait=="Blood Pressure\n(Diastolic)")]<-color_vec[4]
colors[which(plot_df$Trait=="Height")]<-color_vec[5]
colors[which(plot_df$Trait=="Hypothyroidism")]<-color_vec[6]
colors[which(plot_df$Trait=="RBC Count")]<-color_vec[7]
colors[which(plot_df$Trait=="Type 2 Diabetes")]<-color_vec[8]
colors[which(plot_df$Trait=="LDL Cholesterol")]<-color_vec[9]
colors[which(plot_df$Trait=="Triglycerides")]<-color_vec[10]
colors[which(plot_df$Trait=="Glucose")]<-color_vec[11]

names(colors)<-plot_df$Trait


plot_df$`Distance (Kbp)`<-plot_df$distance/1000

p<-ggplot(plot_df, aes(color = Trait,x = `Distance (Kbp)`, y = Enrichment)) +
	geom_line() +
	geom_hline(yintercept = 1, color="grey", linetype = "dashed") +
	theme(
		text = element_text(size = 12),
		strip.background = element_blank(),
		legend.position = "right",
		legend.title = element_blank(),
		legend.text = element_text(size = 9, margin = margin(c(0,0.25,0,0.25), unit = "in")),
		strip.text.y = element_text(size = 10),
		axis.text.x = element_text(angle = 45, hjust = 1),
		aspect.ratio = 1
	) +
	facet_grid(cols = vars(col_facet), rows = vars(annotation), scales = "free_y") +
	labs(x="Distance (Kbp)", y = "Fold enrichment") +
	scale_color_manual(name = "Trait", values = colors, aesthetics = c("colour"), guide=guide_legend(ncol=1)) +
	coord_cartesian(xlim = c(0,250))

pdf(paste0(outdir,"/distance_effect.pdf"), height = 4, width = 8)
p
dev.off()

p<-p+geom_point(aes(text = paste0("# Targeted positive control genes: ",formatC(Targeted_PositiveControl, format="d", big.mark=","))))

## Clean up html output
rm_traces<-1:nrow(unique(plot_df[,.(col_facet,targeting_method,gene_set)]))
p<-style(p, hoverinfo = "none", traces=rm_traces)

outfile<-paste0(outdir,"/distance_effect.html")
p_interactive <- ggplotly(p, dynamicTicks = T)
htmlwidgets::saveWidget(widget=as_widget(p_interactive), file=outfile, selfcontained = F, title="Enrichment by distance", libdir = "html_dependencies")

