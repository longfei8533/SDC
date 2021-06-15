library(data.table)
library(tidyverse)
library(DESeq2)
setwd("/media/CADD/longfei/project/sex_web/")

## Phenotype data
phenotype_dat <- fread("rawdata/toil/TcgaTargetGTEX_phenotype.txt.gz") %>% 
  dplyr::filter(`_gender` %in% c("Female","Male"),`_study` %in% c("TCGA","GTEX")) %>% 
  as.data.frame()

phenotype_dat1 <- dplyr::filter(phenotype_dat, `_study` == "TCGA") %>% # Remove TCGA normal tissues
  dplyr::filter( as.numeric(substr(sample,14,15)) < 10)
phenotype_dat <- dplyr::filter(phenotype_dat , `_study` == "GTEX") %>% 
  rbind(phenotype_dat1)
 
purity_data <- readRDS("result/purity.rds")
phenotype_dat <- merge(phenotype_dat,purity_data, by.x= "sample",by.y= "row.names")

## Read and convert log2(expected_count+1) to count
count_dat <- fread("rawdata/toil/TcgaTargetGtex_gene_expected_count.gz",
                   drop = "sample",
                   nThread = 15) 
count_dat <- 2^count_dat - 1

gene <- fread("rawdata/toil/TcgaTargetGtex_gene_expected_count.gz",
              select = "sample")
gene <- gene$sample %>% substr(.,1,15)


#############################################################



## Function
dds_fun <- function(type,tissue,purity= FALSE){
  if(tissue == "All"){
    sample_chose <- dplyr::filter(phenotype_dat, `_study` == type)
  }else{
    sample_chose <- dplyr::filter(phenotype_dat, `_study` == type , `_primary_site` == tissue  )
  }

  sam <- sample_chose[sample_chose$sample %in% colnames(count_dat),]
  sam$sex <-  factor(sam$`_gender`)
  countData <- count_dat[,match(sam$sample,colnames(count_dat))]
  if(purity){
    dds <- DESeqDataSetFromMatrix(countData, 
                                  colData = sam, design = ~TumorPurity + sex)
  }else{
    dds <- DESeqDataSetFromMatrix(countData, 
                                colData = sam, design = ~ sex)}
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("sex","Female","Male")) 
  res.data <- data.frame(res)
  rownames(res.data) <- gene
  return(res.data)
}

##  Phenotypic statistics and filter
phenotype_stat1 <- dplyr::group_by(phenotype_dat,`_study`,`_primary_site`) %>% 
  dplyr::summarise(n_f = sum(`_gender` == "Female"),n_m = sum(`_gender` == "Male")) %>% 
  dplyr::filter(n_f > 10 & n_m >10) %>% 
  as.data.frame()
phenotype_stat2 <- dplyr::group_by(phenotype_stat1,`_study`) %>% 
  dplyr::summarise(`_primary_site` = "All", 
                   n_f = sum(n_f) ,
                   n_m = sum(n_m) 
                  ) %>% 
  as.data.frame()
phenotype_stat <- rbind(phenotype_stat2,phenotype_stat1)
phenotype_stat$n_sample <- phenotype_stat$n_f + phenotype_stat$n_m

saveRDS(phenotype_stat,file = "result/diff/phenotype_stat.rds")

## Run function
diff_res_list <- list()

for(p in c(FALSE,TRUE)){
  for(i in 1:nrow(phenotype_stat)){
    print(i)
    type <- phenotype_stat$`_study`[i]
    tissue = phenotype_stat$`_primary_site`[i]
  
    res <- dds_fun(type = type,
          tissue = tissue ,
          purity = p)
    res_list <- list(res)
    names(res_list) <- paste(type,tissue,p,sep = "_")
    diff_res_list <- c(diff_res_list,res_list)
  }
}

## Save resultes

saveRDS(diff_res_list,file = "result/diff/diff_res_list.rds")



