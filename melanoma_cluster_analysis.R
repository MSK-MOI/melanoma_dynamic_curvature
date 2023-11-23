# analyze melanoma immunotherapy gene modules

# analysis code written by Kevin Murgas
# Stony Brook University, Dept. of Biomedical Informatics
# advisor: Allen Tannenbaum PhD
# in collaboration with: Rena Elkin PhD, Nadeem Riaz PhD, Emil Saucan PhD, Joseph Deasy PhD

library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(survival)
library(survminer)

setwd('/Users/kevinmurgas/Documents/Tannenbaum Lab/OllivierRicci/melanoma_immunotherapy')

### load data ###

# load in raw expression data
# downloaded from NCBI GEO accession: GSE91061
expr = read_csv("Data/GSE91061_rld_qnorm_HGNC.csv")
colnames(expr)[1] = "Gene_HGNC"
sample_id = colnames(expr)[-1]
pat_id = str_extract(sample_id, "Pt\\d+")
condition = str_extract(sample_id, "(?<=_).+(?=_)")

# load in patient metadata
# downloaded from supplemental files of:
# Riaz et al. Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab. Cell 2017. PMID: 29033130
# paper link: https://www.cell.com/cell/fulltext/S0092-8674(17)31122-4
metadata = read_xlsx("Data/mmc2.xlsx", skip=2) %>%
  mutate(Response2 = if_else(Response %in% c("CR","PR"), "CR/PR", Response))


# load in gene module cluster list
results_dir = 'Results/melanoma_qnorm_pearsonshiftinverse_PD1Bnbrs'
genes = read_csv(file.path(results_dir, "genes.csv"))
clusters = read_csv( file.path(results_dir,"louvain_weighted_clusters.csv"))
# note: cluster starts from 0, module starts from 1 (module = cluster+1)

# match gene names to expression data
expr_match = expr[match(clusters$Gene, expr$Gene_HGNC),-1] %>% as.data.frame
rownames(expr_match) = clusters$Gene


### ANALYZE GENE MODULES ###

# compute difference in expression in matched patients
# matched_pats[on] - match_pats[pre]
match_pats = intersect(pat_id[condition == 'Pre'], pat_id[condition == 'On'])
pre_ind = 0*(1:length(match_pats)); on_ind=pre_ind
for (i in 1:length(match_pats)) {
  p = match_pats[i]
  pre_ind[i] = which((pat_id == p) & (condition == 'Pre'))
  on_ind[i] = which((pat_id == p) & (condition == 'On'))
}
expr_diff = expr_match[,on_ind] - expr_match[,pre_ind] # absolute expression difference
colnames(expr_diff) = match_pats

# normalize gene differences across patients by scaling (divide by SD), but no centering (shift mean) so as to preserve positive and negative changes
expr_diff_scaled = t(scale(t(expr_diff), center=F)) # scaled but not centered


# compute module scores for each sample as average absolute change in expression of genes within each module
nclust = length(unique(clusters$module))
module_score = 0*expr_diff[1:nclust,]; rownames(module_score) = paste0("module",1:nclust)
module_score_scaled = module_score
for (i in unique(clusters$module)) {
  module_genes = clusters$Gene[which(clusters$module == i)]
  print( paste0("module: ",i, " ngenes=",length(module_genes)))
  print( paste(module_genes, collapse=", "))
  
  # compute module score as mean of z-score
  inds = match(module_genes, rownames(expr_diff))
  module_score[i,] = colMeans( expr_diff[inds,] ) # scores for each sample: average change of module genes
  module_score_scaled[i,] = colMeans( expr_diff_scaled[inds,] )
}



### plots ###

## module genes heatmap

# Figure 2A: heatmap of scaled gene expression change
# select matched patients with "Response" reported (therapeutic response: CR=complete remission, PR=partial remis., SD=stable disease, PD=progressive disease, NE=not evaluated)
pat_response = metadata$Response[match(match_pats, metadata$Patient)]
keep_pats = match_pats[which(pat_response != "NE")]

# split by clinical response
# column annotations: response, mutational subtype
top_anno = HeatmapAnnotation(which = 'column',
                             #clinical_response = metadata$Response[match(keep_pats,metadata$Patient)] %>% factor(levels=c("CR","PR","SD","PD","NE")),
                             clinical_response = metadata$Response2[match(keep_pats,metadata$Patient)] %>% factor(levels=c("CR/PR","SD","PD")),
                             mut_subtype = metadata$`Mutational\r\n Subtype`[match(keep_pats,metadata$Patient)] %>% factor(levels = c("BRAF","NF1","RAS","TripleWt")),
                             col = list(
                               #clinical_response = setNames(RColorBrewer::brewer.pal(name = "Set2", n = 5), c("CR","PR","SD","PD","NE")),
                               clinical_response = setNames(c("turquoise", "goldenrod", "red"), c("CR/PR","SD","PD")),
                               mut_subtype = setNames(RColorBrewer::brewer.pal(name = "Accent", n = 5), c("BRAF","NF1","RAS","TripleWt","NA"))
                               )
                             )
# row annotation for module
left_anno = HeatmapAnnotation(which = 'row',
                              module = clusters$module,
                              col = list(module=setNames(RColorBrewer::brewer.pal(name = "Set1", n = nclust), unique(clusters$module))))

pdf(width=8, height=10, file = file.path(results_dir, "module_genes_heatmap.pdf") )
Heatmap(as.matrix(expr_diff_scaled[,keep_pats]), name = "expr. change (scaled)",
        top_annotation = top_anno, left_annotation=left_anno,
        #column_split = metadata$Response2[match(keep_pats,metadata$Patient)] %>% factor(levels=c("CR/PR","SD","PD","NE")),
        #cluster_column_slices = F, show_column_dend = F,
        show_column_names = T, column_names_gp = gpar(fontsize = 8),
        column_title = "Patient ID", column_title_side = "bottom",
        row_split = clusters$module %>% as.numeric, cluster_row_slices=F,
        clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
        show_row_names=F)
dev.off()
as.matrix(expr_diff_scaled[,keep_pats]) %>% rowMeans() %>% hist


## Figure 3A: heatmap of module scores
pdf(width=10, height=3, file = file.path(results_dir, "module_score_heatmap.pdf") )
Heatmap(
  module_score_scaled[,keep_pats], name="module_score",
  top_annotation = top_anno,
  cluster_rows = F,
  column_split = metadata$Response2[match(keep_pats,metadata$Patient)] %>% factor(levels=c("CR/PR","SD","PD","NE")),
  cluster_column_slices = F, cluster_columns=T, show_column_dend=F,
  column_title = "Patient ID", column_title_side="bottom", show_column_names=T, column_names_gp=gpar(fontsize = 10))
dev.off()



##### PATHWAY ANALYSIS #####
# Figure 2B: enrichment dotplots of each gene module
# using clusterProfiler::enrichGO

# convert HGNC gene symbols to Entrez gene ID
library(org.Hs.eg.db)
genes_entrez <-  mapIds(org.Hs.eg.db,
                        keys=clusters$Gene,
                        keytype="SYMBOL",
                        column="ENTREZID",
                        multiVals="first")

## with clusterProfiler and enrichplot
library(clusterProfiler)
library(enrichplot)

# compute GO enrichment of each module
ego_list=list()
for (i in 1:nclust) {
  cluster_genes <- genes_entrez[clusters$module == i]
  print(paste("Cluster",i,": ngenes =",length(cluster_genes)))
  ego <- enrichGO(gene          = cluster_genes,
                  #universe      = genes_entrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  ego_list[[i]] = ego
}
# create enrichment dotplots of each module
plot_list = list()
for (i in 1:nclust) {
  cluster_genes <- genes_entrez[clusters$module == i]
  plot_list[[i]] = dotplot(ego_list[[i]], font.size=8,
                           showCategory=8) + ggtitle(sprintf("module%i (n=%i)",i,length(cluster_genes))) + 
    scale_color_gradient(low="red", high="blue", limits=c(0, 0.001),
                         labels = function(x) sprintf("%.1e", x)) +
    theme(axis.text.y = element_text(size=6))
}
# combine into one plot
p_enrich = ggpubr::ggarrange(plotlist=plot_list, nrow=3, ncol=2,
                             common.legend = TRUE, legend="right")
pdf(width=7, height=10, file = file.path(results_dir, "module_enrichment.pdf") )
p_enrich
dev.off()

# save enrichment results (in separate .csvs)
for (i in 1:nclust) {
  ego_list[[i]] %>% as.data.frame() %>%
    write_csv( file.path(results_dir,paste0("module",i,"_enrichmentGO.csv")))
}



##### MODULE SURVIVAL ANALYSIS #####

# multivariable survival analysis
surv_df = data.frame(pat = match_pats,
                     status = metadata$`Dead/Alive\r\n(Dead = True)`[match(match_pats,metadata$Patient)],
                     time = metadata$`Time to Death\r\n(weeks)`[match(match_pats,metadata$Patient)])

# combine with module scores
surv_modules = t(module_score_scaled) %>% as.data.frame %>% rownames_to_column("pat") %>%
  inner_join(surv_df, by="pat")

# Figure 3B: forest plot of all modules scores
fit_cph = coxph(Surv(time,status) ~ ., data=surv_modules[,-1])
summary(fit_cph)
pdf(width=5, height=6, file = file.path(results_dir, "module_score_forestplot.pdf"), onefile=FALSE  )
ggforest(fit_cph, data=surv_modules,
         cpositions=c(0.02, 0.12, 0.32))
dev.off()


# Figure 3C: Kaplan-Meier of module 6 (split by median)
pdf(width=6, height=6, file = file.path(results_dir, "module6_kmplot.pdf"), onefile=FALSE )
fit6 = survfit(Surv(time,status) ~ module6_grp, data=surv_modules %>%
                 mutate(module6_grp = if_else(module6 > median(module6), "high", "low")))
ggsurvplot(fit6, pval=T, risk.table=T, surv.median.line='hv')
dev.off()


### FIGURE 4: Module 6 and Treatment Response

# Figure 4A: waterfall plot of each module colored by progressive/stable response
plot_df = module_score_scaled %>% t %>% t %>% reshape2::melt() %>%
  dplyr::rename(module=Var1, pat=Var2) %>%
  left_join(data.frame(pat = match_pats,
                       response = metadata$Response[match(match_pats,metadata$Patient)]),
            by="pat") %>%
  mutate(response = factor(response, levels = c("CR","PR","SD","PD","NE"))) %>%
  mutate(response2 = if_else(response %in% c("CR","PR"), "CR/PR",
                             if_else(response=="SD","SD",
                                     if_else(response=="PD","PD",NA))) %>%
           factor(levels=c("CR/PR","SD","PD"))) %>%
  filter(!is.na(response2))
# create plot for each module
plot_list=vector('list',nclust)
for (i in 1:nclust) {
  module_i = sprintf("module%i",i)
  plot_df_i = plot_df %>% dplyr::filter(module == module_i) %>%
    arrange(desc(response2), value) %>% mutate(order = 1:length(unique(plot_df$pat)))
  plot_list[[i]] = ggplot(plot_df_i) + geom_bar(aes(x=reorder(pat, -order), y=value, fill=response2),
                                                stat="identity", width=0.7, position = position_dodge(width=0.4)) +
    theme_classic() + xlab("Patient") + ylab(sprintf('Module %i score',i)) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_fill_manual(values=c("turquoise", "goldenrod", "red"))
  print(paste0("module: ",module_i))
  
  # print mean for each group and non-parametric Kruskal-Wallis ANOVA results
  print(plot_df_i %>% group_by(response2) %>% summarize(m=mean(value)))
  print(kruskal.test(value ~ response2, data = plot_df_i))
}
pdf(width=10, height=3, file = file.path(results_dir, "module6_waterfall.pdf") )
plot_list[[6]]
dev.off()


# test for normality
# fails in all modules but 1, so we will use non-parametric test: Kruskal-Wallis
{
  library(nortest)
  library(car)
  hist_list=vector('list',nclust)
  qq_list=vector('list',nclust)
  for (i in 1:nclust) {
    x = as.numeric(module_score_scaled[i,])
    print(paste0("module",i))
    print(lillie.test( x ))
    hist(x)
    readline(prompt="Press [enter] to continue")
    qqPlot(x)
    readline(prompt="Press [enter] to continue")
  }
  hist_plots = ggarrange(plotlist=hist_list, nrow=1, ncol=6)
  qq_plots = ggarrange(plotlist=qq_list, nrow=1, ncol=6)
  p0 = ggarrange(hist_plots, qq_plots, nrow=2, ncol=1)
}



### Figure 4B/C/D: heatmap of module 6 genes expression change (delta)
# compare CR/PR to PD response groups (and SD?)
pat_response = metadata$Response[match(match_pats, metadata$Patient)]
m6_genes = clusters$Gene[which(clusters$module == 6)]

# this is to determine gene row order based on avg difference in responders (CR/PR) vs non-responders (PD)
expr_diff_m6_CRPR = expr_diff[m6_genes,pat_response %in% c("CR","PR")]
expr_diff_m6_PD = expr_diff[m6_genes,pat_response %in% c("PD")]

avg_delta_CRPR = rowMeans(expr_diff_m6_CRPR)
avg_delta_PD = rowMeans(expr_diff_m6_PD)

gene_order = order(avg_delta_CRPR - avg_delta_PD, decreasing=T)


# keep row order the same for all heatmaps: difference in average expression change responders - non-responders
# left annotation for average change overall in that group
# top annotation for module 6 score (and melanoma subtype, etc?)
for (i_grp in 1:3) {
  # subset patients for each treatment response group
  grp_pats = switch(i_grp,
                    match_pats[pat_response %in% c("CR","PR")],
                    match_pats[pat_response %in% c("SD")],
                    match_pats[pat_response %in% c("PD")])
  expr_diff_m6_grp = expr_diff[m6_genes, grp_pats] %>% as.matrix
  avg_delta_grp = rowMeans(expr_diff_m6_grp)
  
  # column annotations: response, mutational subtype
  top_anno = HeatmapAnnotation(which = 'column',
                               clinical_response = metadata$Response2[match(grp_pats,metadata$Patient)] %>% factor(levels=c("CR/PR","SD","PD")),
                               mut_subtype = metadata$`Mutational\r\n Subtype`[match(grp_pats,metadata$Patient)] %>% factor(levels = c("BRAF","NF1","RAS","TripleWt")),
                               module6_score = anno_barplot(t(module_score_scaled[6,grp_pats]),
                                                            ylim = range(module_score_scaled[6,]),
                                                            gp = gpar(fill = circlize::colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))(t(module_score_scaled[6,grp_pats])),
                                                                      border = NA, lty = "blank")),
                               col = list(
                                 clinical_response = setNames(c("turquoise", "goldenrod", "red"), c("CR/PR","SD","PD")),
                                 mut_subtype = setNames(RColorBrewer::brewer.pal(name = "Accent", n = 5), c("BRAF","NF1","RAS","TripleWt","NA"))
                               ))
  left_anno = HeatmapAnnotation(which='row',
                                avg_diff = avg_delta_grp,
                                col = list(avg_diff = circlize::colorRamp2(c(-1,0,1),c("red","yellow","forestgreen"))))
  
  hm = Heatmap(expr_diff_m6_grp, name="expr. diff.",
          col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
          left_annotation = left_anno,
          row_order = gene_order, row_names_gp = gpar(fontsize=8),
          top_annotation = top_anno,
          column_order = order(-t(module_score_scaled[6,grp_pats])),
          height = length(m6_genes)*unit(3, "mm"),
          width = length(grp_pats)*unit(4, "mm"))
  
  pdf(width=8, height=10, file = file.path(results_dir, paste0("module6_genes_",i_grp,".pdf")) )
  print(hm)
  dev.off()
}


# differential comparison of responders (CR/PR) vs non-responders (PD)
module6_gene_df = data.frame(gene = m6_genes,
                             avg_CRPR = avg_delta_CRPR,
                             avg_PD = avg_delta_PD) %>%
  mutate(avg_difference = avg_CRPR - avg_PD)
for  (ig in 1:length(m6_genes)){
  p = wilcox.test(as.numeric(expr_diff_m6_CRPR[ig,]),
                                     as.numeric(expr_diff_m6_PD[ig,]))$p.value
  module6_gene_df$pval[ig] = p
}
module6_gene_df$FDR = p.adjust(module6_gene_df$pval, method="BH")

sig_ind = which(module6_gene_df$FDR < 0.05)
module6_gene_df[sig_ind,]
