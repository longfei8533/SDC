library(tidyverse)
library(ChAMP)
library(data.table)


pd <- fread("rawdata/pancan/pancan_Curated clinical data") %>%
  dplyr::filter(gender %in% c("FEMALE", "MALE"))

## download file

for ( i in unique(pd$`cancer type abbreviation`)){
  print(i)
  url <- paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.",
               i,
               ".sampleMap%2FHumanMethylation450.gz"
               )
  destfile <- paste0("rawdata/DNA_Methyation/TCGA.",i,".sampleMap%2FHumanMethylation450.gz")
  if(!file.exists(destfile)){
  download.file(url = url, destfile = destfile)}
}

## function
diff_fun <- function(dat, pd) {

  myLoad <- champ.filter(beta = dat, pd = pd)
  myLoad <- champ.impute(
    beta = myLoad$beta,
    pd = myLoad$pd,
    method = "Combine",
    k = 10,
    ProbeCutoff = 0.25,
    SampleCutoff = 0.2
  )

  myNorm <- champ.norm(myLoad$beta,cores = 4)

  myDMP <- champ.DMP(
    beta = myNorm,
    pheno = myLoad$pd$Sample_Group,
    compare.group = c("FEMALE", "MALE"),
  )

  myDMR <- champ.DMR(
    beta = myNorm,
    pheno = myLoad$pd$Sample_Group,
    compare.group = c("FEMALE", "MALE"),
    cores = 4
  )

  myBlock <- champ.Block(
    beta = myNorm,
    pheno = myLoad$pd$Sample_Group,
    cores = 4
  )
  list(DMP = myDMP, DMR = myDMR, myBlock = myBlock)
}


diff_res_list <- list()
#

for(i in unique(pd$`cancer type abbreviation`)){
  print(i)
  datname <- paste0("rawdata/DNA_Methyation/TCGA.",i,".sampleMap%2FHumanMethylation450.gz")
  dat <- fread(datname) %>%
    column_to_rownames(var = "sample")

  pd1 <- pd[pd$sample %in% colnames(dat), ] %>%
    as.data.frame() %>%
    dplyr::select(sample, gender)
  names(pd1) <- c("Sample_Name", "Sample_Group")
  if(sum(pd1$Sample_Group == "FEMALE") < 20 | sum(pd1$Sample_Group == "MALE") < 20 ){next}
  dat1 <- dat[, match(pd1$Sample_Name, colnames(dat))] %>% as.matrix()
  
  res <- tryCatch(diff_fun(dat= dat1 , pd = pd1),
                  error  = function(e){
                    NULL
                  })
  if(is.null(res)){next}
  
  diff_res_list[[i]] <- res
  
}



xx <- list()

for( i in names(diff_res_list)){
  x <- list(DMP = diff_res_list[[i]][["DMP"]][[1]],
       DMR = diff_res_list[[i]][["DMR"]][[1]],
       myBlock  = diff_res_list[[i]][["myBlock"]][[1]])
  xx[[i]] <- x
}

saveRDS(xx,file = "result/diff/DNAmeth_diff_res_list.rds")



