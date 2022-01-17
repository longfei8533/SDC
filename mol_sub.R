library(tidyverse)
library(data.table)
tcga_clinical_dat <- read.table("rawdata/pancan/pancan_Curated clinical data",
                                sep = "\t",header = T)
subtypes <- fread("rawdata/pancan/TCGASubtype.20170308.tsv",
                                sep = "\t",header = T) %>% as.data.frame()
subtypes <- subtypes[-grep("Normal",subtypes$Subtype_Selected),]
subtypes$pan.samplesID <- substr(subtypes$pan.samplesID,1,12)
subtypes$sex <- tcga_clinical_dat$gender[match(subtypes$pan.samplesID,
                                               tcga_clinical_dat$X_PATIENT)]


subtypes$Subtype_Selected <- lapply(subtypes$Subtype_Selected,function(x){
  strsplit(x,split = "[.]")[[1]][2]
}) %>% unlist()
subtypes <- dplyr::filter(subtypes,sex %in% c("FEMALE","MALE"),
                             !Subtype_Selected %in% c("NA","-")
                             ) 

filter_cancer <- group_by(subtypes,cancer.type) %>% 
  summarise(nf = sum(sex == "FEMALE"),nm = sum(sex == "MALE")) %>% 
  as.data.frame()

filter_cancer <- filter_cancer$cancer.type[(filter_cancer$nf > 10) & (filter_cancer$nm > 10)]

subtypes <- subtypes[subtypes$cancer.type %in% filter_cancer,]

group_by(subtypes,cancer.type) %>% 
  summarise(type = names(table(Subtype_Selected)),
    n = table(Subtype_Selected)
            ) %>% View()
saveRDS(subtypes,file = "result/mol_sub_dat.rds")

