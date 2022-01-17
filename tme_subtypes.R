library(tidyverse)
library(data.table)
# Please pre-convert XLXS files to TAB split files. (Summary of TME subtype for each analyzed tumor (Table S6))
subtypes <- fread("rawdata/pancan/pan-cancer microenvironment subtypes.txt",
                  sep = "\t",header = T) %>% as.data.frame()
table(subtypes$MFP)
table(subtypes$Gender)


filter_cancer <- group_by(subtypes,TCGA_project) %>% 
  summarise(nf = sum(Gender == "F"),nm = sum(Gender == "M")) %>% 
  as.data.frame()
filter_cancer <- filter_cancer$TCGA_project[(filter_cancer$nf > 10) & (filter_cancer$nm > 10)]


subtypes <- subtypes[subtypes$TCGA_project %in% filter_cancer,]

saveRDS(subtypes,file = "result/tme_sub_dat.rds")
