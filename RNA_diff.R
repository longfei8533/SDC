library(data.table)
library(tidyverse)
library(DESeq2)

phenotype_dat <- fread("rawdata/toil/TcgaTargetGTEX_phenotype.txt.gz") %>% 
  dplyr::filter(`_gender` %in% c("Female","Male"),`_study` %in% c("TCGA","GTEX")) %>% 
  as.data.frame()

tcga_gtex_paired <- readRDS("rawdata/toil/tcga_gtex_paired.rds") %>% 
  dplyr::filter(!(substr(sample,1,4) == "TCGA" & type2 == "normal"))

phenotype_dat <- merge(tcga_gtex_paired,phenotype_dat,by = "sample")
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


##  Phenotypic statistics and filter
phenotype_stat1 <- dplyr::group_by(phenotype_dat,`_study`,tissue) %>% 
  dplyr::summarise(n_f = sum(`_gender` == "Female"),n_m = sum(`_gender` == "Male")) %>% 
  dplyr::filter(n_f > 10 & n_m >10) %>% 
  as.data.frame()
phenotype_stat2 <- dplyr::group_by(phenotype_stat1,`_study`) %>% 
  dplyr::summarise(`tissue` = "All", 
                   n_f = sum(n_f) ,
                   n_m = sum(n_m) 
  ) %>% as.data.frame()
phenotype_stat <- rbind(phenotype_stat2,phenotype_stat1)
phenotype_stat$n_sample <- phenotype_stat$n_f + phenotype_stat$n_m

saveRDS(phenotype_stat,file = "result/diff/phenotype_stat.rds")

##

dds_fun <- function(count_dat,phenotype_dat,study,tissue,purity){
  if(tissue == "All"){
    sample_chose <- dplyr::filter(phenotype_dat, .data$`_study` == study)
  }else{
    sample_chose <- dplyr::filter(phenotype_dat, .data$`_study` == study , .data$tissue == {{tissue}} )
  }
  
  sam <- sample_chose[sample_chose$sample %in% colnames(count_dat),]
  sam$sex <- factor(sam$`_gender`)
  countData <- count_dat[,sam$sample,with=FALSE] %>% 
    apply(2, as.integer)
  
  
  if(purity){
  dds <- DESeqDataSetFromMatrix(countData, 
                                  colData = sam, design = ~TumorPurity + sex)
  }else{
  dds <- DESeqDataSetFromMatrix(countData, 
                                  colData = sam, design = ~ sex)
  }
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("sex","Female","Male")) %>% 
    as.data.frame()
  rownames(res) <- gene
  return(res)
}

diff_res_list <- list()
for(purity in c(FALSE,TRUE)){
  for(i in 1:nrow(phenotype_stat)){
    study <- phenotype_stat$`_study`[i]
    tissue <-  phenotype_stat$tissue[i]
    res <- dds_fun(
      count_dat,phenotype_dat,study = study,
      tissue = tissue, purity = purity )
    
    list_name <- paste(study,tissue,purity,sep = "_")
    diff_res_list[[list_name]] <- res
    print(list_name)
  }
}

saveRDS(diff_res_list,file = "result/diff/diff_res_list.rds")
