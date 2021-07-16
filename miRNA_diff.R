library(TCGAbiolinks)
library(tidyverse)
library(data.table)
library(DESeq2)
purity_dat <- readRDS("../result/purity.rds")
anno_dat <- fread("../rawdata/Survival_SupplementalTable_S1_20171025_xena_sp") %>% 
  dplyr::filter(as.numeric(substr(sample,14,15)) < 10) %>% 
  dplyr::filter(gender %in% c("FEMALE","MALE")) %>% 
  as.data.frame()

purity_data <- readRDS("../result/purity.rds")
anno_dat <- base::merge(anno_dat,purity_data, by.x= "sample",by.y= 0)

## data_download
for (i in unique(anno_dat$`cancer type abbreviation`)[1:33]) {
  print(i)
  CancerProject <- paste0("TCGA-", i)
  DataDirectory <- paste0("../rawdata/PANCAN/miRNA/", CancerProject)
  FileNameData <- paste0(DataDirectory, "_", "miRNA_gene_quantification", ".rda")
  if(file.exists(FileNameData)){next}

  query <- GDCquery(
    project = CancerProject,
    data.category = "Gene expression",
    data.type = "miRNA gene quantification",
    data.format = "txt",
    platform = "Illumina HiSeq",
    file.type = "13.mirna.quantification.txt",
    legacy = TRUE
  )

  GDCdownload(
    query = query,
    directory = DataDirectory
  )
  dataAssy <- GDCprepare(
    query = query,
    save = TRUE,
    save.filename = FileNameData,
    summarizedExperiment = TRUE,
    directory = DataDirectory
  )
  
}

## 
dds_fun <- function(anno,dataAssy,purity){
 
  anno$gender <- as.factor(anno$gender)
  dataAssy <- dataAssy[,match(anno$sample,colnames(dataAssy))]
  
  if(purity){
    dds <- DESeqDataSetFromMatrix(countData = dataAssy, 
                                  colData = anno, design = ~TumorPurity + gender)
  }else{
    dds <- DESeqDataSetFromMatrix(countData = dataAssy, 
                                  colData = anno, design = ~ gender)}
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("gender","FEMALE","MALE")) 
  res.data <- data.frame(res)
  rownames(res.data) <- mirname
  return(res.data)
}

diff_res_list <- list()

for(p in c(FALSE,TRUE)){
  for(i in unique(anno_dat$`cancer type abbreviation`)){
    print(i)
    load(paste0("../rawdata/PANCAN/miRNA/TCGA-",i,"_miRNA_gene_quantification.rda"))
    dataAssy <- as.data.frame(data)
    mirname <- dataAssy$miRNA_ID
    dataAssy <- dplyr::select(dataAssy, starts_with("read_count"))
    colnames(dataAssy) <- substr(colnames(dataAssy), 12, 26)
    anno <- anno_dat[anno_dat$sample %in%  colnames(dataAssy),]
    if(nrow(anno)<10){next}
    phenotype_stat <- dplyr::group_by(anno,`cancer type abbreviation`) %>%
      dplyr::summarise(n_f = sum(`gender` == "FEMALE"), n_m = sum(`gender` == "MALE")) %>%
      as.data.frame()
    if(phenotype_stat$n_f < 10 | phenotype_stat$n_m < 10){next}
    res <- dds_fun(anno,dataAssy,purity = p)
    res_list <- list(res)
    names(res_list) <- paste(i,p,sep = "_")
    diff_res_list <- c(diff_res_list,res_list)
  }
}

## Save resultes

saveRDS(diff_res_list,file = "../result/diff/miRNA_diff_res_list.rds")
