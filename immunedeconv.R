# an R package for unified access to computational methods for estimating 
# immune cell fractions from bulk RNA sequencing data. 
# https://github.com/icbi-lab/immunedeconv
# https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html
# The method can be one of
est_method <- c("quantiseq", "timer", "cibersort", "cibersort_abs", "mcp_counter", "xcell","epic")

# The input data is a gene × sample gene expression matrix. In general values should be
# 1. TPM-normalized
# 2. not log-transformed.

library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(immunedeconv) 
library(tibble)

## set path
setwd("/media/CADD/longfei/project/sex_web/")
set_cibersort_binary("/media/CADD/longfei/project/sex_web/code/CIBERSORT.R")
set_cibersort_mat("/media/CADD/longfei/project/sex_web/code/LM22.txt")

## read and conver log2(tpm+0.001) to tpm
log2_tpm_dat <- fread("rawdata/toil/TcgaTargetGtex_rsem_gene_tpm.gz",
  drop = "sample"
  )
tpm_dat <- 2^(log2_tpm_dat)-0.001
rm(log2_tpm_dat)

ensg <- fread("rawdata/toil/TcgaTargetGtex_rsem_gene_tpm.gz",data.table=FALSE,
  select = "sample"
  ) 
ensg <- ensg$sample %>% substr(.,1,15)

#############################################################

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

tpm_matrix <- tpm_dat[,1:10]

# imm_res_list <- list()
# for(i in est_method ){
#   if(i %in% c("timer","cibersort", "cibersort_abs")){
#     print(paste(i,"Pass!"))
#     next()
#   }else{
#     res <- immunedeconv::deconvolute(tpm_matrix, i) 
#     print(paste(i,Sys.time()))
#     }
#   res <- list(res)
#   names(res) <- i
#   imm_res_list <- c(imm_res_list,res)
# }
  
imm_res_list <- immunedeconv::deconvolute(tpm_matrix, "cibersort")
saveRDS(imm_res_list,file = "result/imm/imm_res_list.rds")



