library(data.table)
library(tidyverse)

## CNV data
cnv_dat <- fread("rawdata/pancan/TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz") %>% 
  tibble::column_to_rownames(.,"Sample")


## Phenotype data
anno_dat <- fread("rawdata/pancan/pancan_Curated clinical data") %>% 
  dplyr::filter(as.numeric(substr(sample,14,15)) < 10) %>%   # remove no-tumor sample
  dplyr::filter(sample %in% colnames(cnv_dat))


sex_sum <- dplyr::group_by(anno_dat,`cancer type abbreviation`) %>% 
  dplyr::summarise(sum_female = sum(gender == "FEMALE"),
                   sum_male = sum(gender == "MALE")) %>% 
  dplyr::filter(sum_female > 10 & sum_male > 10)

anno_dat <- anno_dat[anno_dat$`cancer type abbreviation` %in% sex_sum$`cancer type abbreviation`,]

cnv_dat <- cnv_dat[,match(anno_dat$sample,colnames(cnv_dat))]

##


function(data,tissue){
  
}

# 
diff_res_list <- list()
for( i in unique(anno_dat$`cancer type abbreviation`)){
  print(i)
  t_dat <- cnv_dat[,anno_dat$`cancer type abbreviation` == i]
  a_d <- anno_dat[anno_dat$`cancer type abbreviation` == i,]
  res <- apply(t_dat,1,function(x){
    dat <- data.frame(var = x , gender = a_d$gender)
    tt <- t.test(var ~ gender ,data = dat)
    pv <- tt$p.value
    ef <- tt$estimate["mean in group FEMALE"] - tt$estimate["mean in group MALE"]
    if(is.na(pv)){pv <- 1}
    if(is.na(ef)){ef <- 0}
    c(ef,pv)
  }) %>% t()
  colnames(res) <- c("effect(F-M)","pvalue")
  res_list <- list(res)
  names(res_list) <- i
  diff_res_list <- c(diff_res_list,res_list)
}

cnv_diff_res <- list("sample_anno" = anno_dat , "ttest_res" = diff_res_list)

saveRDS(cnv_diff_res,file = "result/diff/cnv_diff_res.rds")




