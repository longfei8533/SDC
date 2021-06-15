library(dplyr)
library(data.table)
library(org.Hs.eg.db)

setwd("/media/CADD/longfei/project/sex_web/")

## read and conver log2(tpm+0.001) to tpm
tpm_dat <- fread("rawdata/toil/TcgaTargetGtex_rsem_gene_tpm.gz") %>% # log2(tpm+0.001)
  tibble::column_to_rownames(.,"sample")

## gene id conversion
ensg <- rownames(tpm_dat) %>% substr(.,1,15)

gene_entrezid <- AnnotationDbi::mapIds(x <- org.Hs.eg.db,
                                     keys = ensg,
                                     column = "ENTREZID", 
                                     keytype = "ENSEMBL",
                                     multiVals="first"
)


tpm_dat <- tpm_dat[!is.na(gene_entrezid) & !duplicated(gene_entrezid),] %>% 
  as.data.frame()
gene_entrezid <- na.omit(gene_entrezid) %>% unique()
rownames(tpm_dat) <- gene_entrezid

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

data.table::fwrite(tcga_tpm_dat,"tempdata/tcga_tpm_entrezid.gz",sep="\t",quote = FALSE,
                   row.names = TRUE,compress = "gzip")







