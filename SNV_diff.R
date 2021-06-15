library(maftools)
library(data.table)
library(tidyverse)
dat <- fread("rawdata/pancan/mc3.v0.2.8.PUBLIC.maf.gz")
dat$Tumor_Sample_Barcode <- substr(dat$Tumor_Sample_Barcode,1,15)

anno_dat <- fread("rawdata/pancan/pancan_Curated clinical data") %>% 
  dplyr::filter(as.numeric(substr(sample,14,15)) < 10) %>%   # remove no-tumor sample
  dplyr::filter(sample %in% dat$Tumor_Sample_Barcode)


sex_sum <- dplyr::group_by(anno_dat,`cancer type abbreviation`) %>% 
  dplyr::summarise(sum_female = sum(gender == "FEMALE"),
                   sum_male = sum(gender == "MALE")) %>% 
  dplyr::filter(sum_female > 10 & sum_male > 10)

anno_dat <- anno_dat[anno_dat$`cancer type abbreviation` %in% sex_sum$`cancer type abbreviation`,]


# 1) Comparing two cohorts (Female Male)

compare_fun <- function(cancer_type,maf_dat,anno_dat){
  print(cancer_type)
  a <- anno_dat[anno_dat$`cancer type abbreviation` == cancer_type,]
  f_barcode <- a$sample[a$gender == "FEMALE"]
  m_barcode <- a$sample[a$gender == "MALE"]
  m_f <- read.maf(maf_dat[maf_dat$Tumor_Sample_Barcode %in% f_barcode,],verbose = FALSE)
  m_m <- read.maf(maf_dat[maf_dat$Tumor_Sample_Barcode %in% m_barcode,],verbose = FALSE)
  pt.vs.rt <- mafCompare(m1 = m_f, 
                      m2 = m_m, 
                      m1Name = "FEMALE", 
                      m2Name = "MALE",
                      minMut = 5)
  
  result <- dplyr::filter(pt.vs.rt$results,.data$pval < 0.05)
  ## plot
  
  if(nrow(result)>1){
    if(nrow(result)>10){result <- result[1:10,]}
    pt10 <- list(results = result , SampleSummary = pt.vs.rt$SampleSummary)
    if("TP53" %in% result$Hugo_Symbol){
      pt0.05gene <- result$Hugo_Symbol
    }else{
      pt0.05gene <- c(result$Hugo_Symbol ,"TP53")
    } # "TP53" was used to plot correctly, because of a bug of maftools

    png(filename = paste0("result/diff/snv/",cancer_type,"_forestPlot.png"),width =3200 ,height = 2000,res = 300)
    forestPlot(mafCompareRes = pt10)
    dev.off()
    # 
    #png(filename = paste0("result/diff/snv/",cancer_type,"_coOncoplot.png"),width =3200 ,height = 2000,res= 300)
    # coOncoplot(m1 = m_f,
    #            m2 = m_m,
    #            genes = pt0.05gene,
    #            m1Name = "FEMALE",
    #            m2Name = "MALE")
    #dev.off()
    # 
    png(filename = paste0("result/diff/snv/",cancer_type,"_coBarplot.png"),width =3200 ,height = 2000,res = 300)
    coBarplot(m1 = m_f,
               m2 = m_m,
               genes = pt0.05gene,
               m1Name = "FEMALE",
               m2Name = "MALE")
    dev.off()
  }
  ## muttaion burden
  f_b <- getSampleSummary(m_f) %>% 
    select(.data$Tumor_Sample_Barcode,.data$total) %>% 
    mutate(gender = "FEMALE")
  m_b <- getSampleSummary(m_m) %>% 
    select(.data$Tumor_Sample_Barcode,.data$total) %>% 
    mutate(gender = "MALE")
  b_res <- rbind(f_b,m_b)
  
  list(pt.vs.rt = pt.vs.rt, burden = b_res)
}


snv_res <- list()
for(i in unique(anno_dat$`cancer type abbreviation`) ){
  rc <- compare_fun(i,dat,anno_dat = anno_dat)
  rcl <- list(rc)
  snv_res <- c(snv_res,rcl)
}

names(snv_res) <- unique(anno_dat$`cancer type abbreviation`)

saveRDS(snv_res,file = "result/diff/snv/snv_res.rds")


