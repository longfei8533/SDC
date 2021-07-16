## Reference DOI: 10.1038/ncomms9971

library(dplyr)
library(estimate)
library(data.table)
library(org.Hs.eg.db)
source("function/estimate_.R")


setwd("/media/CADD/longfei/project/sex_web/")
tpm_dat <- fread("rawdata/toil/TcgaTargetGtex_rsem_gene_tpm.gz",data.table=FALSE) 
  
#phenotype_dat <- fread("rawdata/toil/TcgaTargetGTEX_phenotype.txt.gz")

ensg <- tpm_dat$sample %>% substr(.,1,15)

gene_symbol <- AnnotationDbi::mapIds(x <- org.Hs.eg.db,
                                     keys = ensg,
                                     column = "SYMBOL", 
                                     keytype = "ENSEMBL",
                                     multiVals="first"
                                     )
tpm_matrix <- tpm_dat[!is.na(gene_symbol) & !duplicated(gene_symbol),-1]

gene_symbol <- na.omit(gene_symbol) %>% unique()
rownames(tpm_matrix) <- gene_symbol

filter_mat <- filterCommonGenes_(tpm_matrix)

estimateRes<- estimateScore_(filter_mat, 
                            platform="affymetrix")

estimateRes <- t(estimateRes)

saveRDS(estimateRes,file = "result/purity.rds")

