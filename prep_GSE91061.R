# Data preparation of GSE91061
# Melanoma Immunotherapy Response

# here, we convert Entrez gene IDs to HGNC gene symbols
# we also apply a basic filter to remove lowly expressed genes (<10 counts in >90% of samples)
# then apply quantile normalization

library(tidyverse)
library(readxl)
library(DESeq2)

setwd('melanoma_immunotherapy/Data')

### load sample/patient metadata
meta_file = 'mmc2.xlsx'
metadata = read_xlsx(meta_file, skip=2)

### load data (raw counts)
cts = read_csv('GSE91061_BMS038109Sample.hg19KnownGene.raw.csv')
rld = read_csv('GSE91061_BMS038109Sample.hg19KnownGene.rld.csv')
sample_names = colnames(cts)[-1]
genes_entrez = cts %>% pull(1) %>% as.character

### convert Entrez ID to HGNC symbol
library(org.Hs.eg.db)
genes_hgnc <-  mapIds(org.Hs.eg.db,
                      keys=genes_entrez,
                      column="SYMBOL",
                      keytype="ENTREZID",
                      multiVals="first")
c(sum(!is.na(genes_hgnc)), sum(is.na(genes_hgnc))) # mapped 22068 genes, 119 NA


### convert data to matrix
# counts
cts_mat = cts %>% dplyr::select(-`...1`) %>% as.matrix
rownames(cts_mat) <- genes_hgnc
cts_mat = cts_mat[-is.na(genes_hgnc),] # remove genes mapping to NA
# rlog data
rld_mat = rld %>% dplyr::select(-`...1`) %>% as.matrix
rownames(rld_mat) <- genes_hgnc
rld_mat = rld_mat[-is.na(genes_hgnc),] # remove genes mapping to NA


### pre-filter low-expressed genes
nrow(cts_mat) # 22186 genes
keep = rowMeans(cts_mat >= 10) >= 0.1 # at least 10% samples with a count of 10 or higher
cts_mat = cts_mat[keep,]
rld_mat = rld_mat[keep,]
nrow(cts_mat) # 18413 genes

# quantile normalize rlog data
library(preprocessCore)
rld_qnorm = normalize.quantiles(rld_mat, keep.names = T)


### save data with HGNC symbols
# write_csv(cts_mat %>% as.data.frame %>% rownames_to_column("Gene_HGNC"),
#           "GSE91061_cts_HGNC.csv")
# write_csv(rld_mat %>% as.data.frame %>% rownames_to_column("Gene_HGNC"),
#           "GSE91061_rld_HGNC.csv")
write_csv(rld_qnorm %>% as.data.frame %>% rownames_to_column("Gene_HGNC"),
          "GSE91061_rld_qnorm_HGNC.csv")
