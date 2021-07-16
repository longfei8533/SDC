library(IOBR)
library(EPIC)
library(estimate) 
library(MCPcounter)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(data.table)
library(org.Hs.eg.db)


## read and conver log2(tpm+0.001) to tpm
log2_tpm_dat <- fread("rawdata/toil/TcgaTargetGtex_rsem_gene_tpm.gz") %>% 
  tibble::column_to_rownames(.,"sample")
  
tpm_dat <- 2^(log2_tpm_dat)-0.001

rm(log2_tpm_dat)

## gene id conversion
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

anno_dat <- fread("rawdata/pancan/pancan_Curated clinical data") %>% 
  dplyr::filter(as.numeric(substr(sample,14,15)) < 10) %>%   # remove no-tumor sample
  dplyr::filter(sample %in% colnames(tpm_dat))


sex_sum <- dplyr::group_by(anno_dat,`cancer type abbreviation`) %>% 
  dplyr::summarise(sum_female = sum(gender == "FEMALE"),
                   sum_male = sum(gender == "MALE")) %>% 
  dplyr::filter(sum_female > 10 & sum_male > 10)

anno_dat <- anno_dat[anno_dat$`cancer type abbreviation` %in% sex_sum$`cancer type abbreviation`,]

tcga_tpm_dat <- tpm_dat[,match(anno_dat$sample,colnames(tpm_dat))]

saveRDS(anno_dat,file = "result/imm/imm_pheno.rds")

## decode TME

xcell <-  deconvo_xcell(tcga_tpm_dat, arrays = FALSE)
saveRDS(xcell,"result/imm/xcell.rds")


mcp <-  deconvo_mcpcounter(tcga_tpm_dat)
saveRDS(mcp,"result/imm/mcp.rds")

epic <-  deconvo_epic(tcga_tpm_dat, tumor = TRUE)
saveRDS(epic,"result/imm/epic.rds")

cibersort <-  deconvo_cibersort(tcga_tpm_dat, absolute = FALSE, arrays =FALSE,perm = 200)
saveRDS(cibersort,"result/imm/cibersort.rds")

cibersort_abs <-  deconvo_cibersort(tcga_tpm_dat, absolute = TRUE, 
                                  abs_method = "sig.score",arrays = FALSE,perm = 200) 
saveRDS(cibersort_abs,"result/imm/cibersort_abs.rds")

ips <-  deconvo_ips(tcga_tpm_dat, plot = FALSE)
saveRDS(ips,"result/imm/ips.rds")

quantiseq <-  deconvo_quantiseq(tcga_tpm_dat, tumor=TRUE, arrays=FALSE, scale_mrna=TRUE) 
saveRDS(quantiseq,"result/imm/quantiseq.rds")

estimate <- deconvo_tme(eset = tcga_tpm_dat, method = "estimate")
saveRDS(estimate,"result/imm/estimate.rds")

timer_available_cancers <- c('kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg',
                             'lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct',
                             'ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca',
                             'uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol')
timer_tpm <- tcga_tpm_dat[,anno_dat$`cancer type abbreviation` %in% toupper(timer_available_cancers)]
timer_anno <- anno_dat[anno_dat$`cancer type abbreviation` %in% toupper(timer_available_cancers),]
timer <-  deconvo_timer(timer_tpm, indications = timer_anno$`cancer type abbreviation`)
saveRDS(timer,"result/imm/timer.rds")


tme_combine <- cibersort %>%
  #inner_join(.,cibersort_abs,by= "ID") %>%
  inner_join(.,mcp,by          = "ID") %>%
  inner_join(.,xcell,by        = "ID") %>%
  inner_join(.,epic,by         = "ID") %>%
  inner_join(.,timer,by        = "ID") %>%
  inner_join(.,quantiseq,by    = "ID") %>%
  inner_join(.,ips,by          = "ID") %>% 
  inner_join(.,estimate,by       = "ID")

print("Successful!")
print(Sys.time())

## save

saveRDS(tme_combine,"result/imm/tem_combine.rds")





