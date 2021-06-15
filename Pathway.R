library(data.table)
library(tidyverse)

setwd("/media/CADD/longfei/project/sex_web/")

ssgsea_dat <- fread("rawdata/Pathway/PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level.txt.gz") %>% as.data.frame()
pw_s <- ssgsea_dat[,1]
ssgsea_dat <- ssgsea_dat[,-1]

zscore_dat <- fread("rawdata/Pathway/PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level_Z.txt.gz") %>% as.data.frame()
pw_z <- zscore_dat[,1]
zscore_dat <- zscore_dat[,-1]

all(pw_s == pw_z)
all(colnames(ssgsea_dat) == colnames(zscore_dat))

## Phenotype data
anno_dat <- fread("rawdata/pancan/pancan_Curated clinical data") %>% 
  dplyr::filter(as.numeric(substr(sample,14,15)) < 10) %>%   # remove no-tumor sample
  dplyr::filter(sample %in% colnames(zscore_dat))


sex_sum <- dplyr::group_by(anno_dat,`cancer type abbreviation`) %>% 
  dplyr::summarise(sum_female = sum(gender == "FEMALE"),
                   sum_male = sum(gender == "MALE")) %>% 
  dplyr::filter(sum_female > 10 & sum_male > 10)

anno_dat <- anno_dat[anno_dat$`cancer type abbreviation` %in% sex_sum$`cancer type abbreviation`,]


tcga_tpm_dat <- tpm_dat[,match(anno_dat$sample,colnames(tpm_dat))]

ssgsea_dat <- ssgsea_dat[,match(anno_dat$sample,colnames(ssgsea_dat))]
zscore_dat <- zscore_dat[,match(anno_dat$sample,colnames(zscore_dat))]

## diff

diff_res_list <- list()
for( i in unique(anno_dat$`cancer type abbreviation`)){
  print(i)
  t_dat <- zscore_dat[,anno_dat$`cancer type abbreviation` == i]
  a_d <- anno_dat[anno_dat$`cancer type abbreviation` == i,]
  res <- apply(t_dat,1,function(x){
    dat <- data.frame(var = x , gender = a_d$gender)
    tt <- wilcox.test(var ~ gender ,data = dat)
    pv <- tt$p.value
    ef <- group_by(dat,gender) %>% 
      summarise(m = mean(var,na.rm=TRUE))
    ef <- ef$m[ef$gender == "FEMALE"] - ef$m[ef$gender == "MALE"]
    if(is.na(pv)){pv <- 1}
    if(is.na(ef)){ef <- 0}
    c(ef,pv)
  }) %>% t() %>% as.data.frame()
  rownames(res) <- pw_s
  colnames(res) <- c("effect(F-M)","pvalue")
  res$p.adjust <- p.adjust(res$pvalue,method = "BH")
  diff_res_list[[i]] <- res
}

pathway_diff_res <- list("sample_anno" = anno_dat , "wilcox.test_res" = diff_res_list)

saveRDS(pathway_diff_res,file = "result/diff/pathway_diff_res.rds")

