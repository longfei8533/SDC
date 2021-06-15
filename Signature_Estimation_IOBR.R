# IOBR integrates 255 published signature gene sets, involving tumor microenvironment, 
# tumor metabolism, m6A, exosomes, microsatellite instability, and tertiary lymphoid structure. 
# Running the function signature_collection_citation to attain the source papers. 
# The function signature_collection returns the detail signature genes of all given signatures.

# Obtain the included signatures first.The signature collection is mainly classified into 3 categories: 
#   TME-associated, tumor-metabolism, and tumor-intrinsic signatures.

library(IOBR)
library(data.table)
library(dplyr)
library(org.Hs.eg.db)

setwd("/media/CADD/longfei/project/sex_web/")

tpm_dat <- fread("rawdata/toil/TcgaTargetGtex_rsem_gene_tpm.gz") %>% 
  tibble::column_to_rownames(.,"sample")   ## log2(tpm+0.001)


ensg <- rownames(tpm_dat) %>% substr(.,1,15)

gene_symbol <- AnnotationDbi::mapIds(x <- org.Hs.eg.db,
                                     keys = ensg,
                                     column = "SYMBOL", 
                                     keytype = "ENSEMBL",
                                     multiVals="first"
)

tpm_dat <- tpm_dat[!is.na(gene_symbol) & !duplicated(gene_symbol),] %>% 
  as.data.frame()
gene_symbol <- na.omit(gene_symbol) %>% unique()
rownames(tpm_dat) <- gene_symbol


## sample annotation and filter

anno_dat <- readRDS("result/imm/imm_pheno.rds")

tcga_tpm_dat <- tpm_dat[,match(anno_dat$sample,colnames(tpm_dat))]

sig_res<-calculate_sig_score(pdata           = NULL,
                             eset            = tcga_tpm_dat,
                             signature       = signature_collection,
                             method          = "integration",
                             mini_gene_count = 2,
                             adjust_eset = TRUE)
saveRDS(sig_res,"result/imm/signature_collection.rds")
