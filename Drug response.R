library(data.table)
library(tidyverse)
setwd("/media/CADD/longfei/project/sex_web")

# ccle
ccle_drug_dat <- fread("rawdata/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv") %>% 
  dplyr::select(cell_line = `Primary Cell Line Name`, 
                compound = Compound,
                target = Target,
                IC50 = `IC50 (uM)`
                )  %>% 
  dplyr::mutate(IC50 = log(IC50))

ccle_anno_dat <- fread("rawdata/CCLE/Cell_lines_annotations_20181226.txt") %>% 
  dplyr::mutate(dataset = "CCLE") %>% 
  dplyr::select(dataset,
                cell_line = Name,
                tcga_code,
                gender = Gender
                )
ccle_dat <- dplyr::inner_join(ccle_anno_dat,ccle_drug_dat)
ccle_dat$gender <- mapvalues(ccle_dat$gender,c("female","male",""),
                             c("FEMALE","MALE","Unknown"))

  

# gdsc 

gdsc1_drug_dat <- fread("rawdata/GDSC/GDSC1_fitted_dose_response_25Feb20.csv")
gdsc2_drug_dat <- fread("rawdata/GDSC/GDSC2_fitted_dose_response_25Feb20.csv")
gdsc_drug_dat <- rbind(gdsc1_drug_dat,gdsc2_drug_dat) %>% 
  dplyr::select(dataset = DATASET, 
                cell_line = CELL_LINE_NAME,
                compound = DRUG_NAME,
                target = PUTATIVE_TARGET,
                IC50 = LN_IC50
                )

gdsc_anno_dat <- fread("rawdata/GDSC/Cell_Lines_Details.csv")
gdsc_gender_dat <- fread("rawdata/GDSC/CosmicSample.tsv.gz") %>% 
  #dplyr::filter(sample_type == "cell-line") #%>% 
  dplyr::select(sample_id, gender)

gdsc_anno_dat <- dplyr::left_join(gdsc_anno_dat,gdsc_gender_dat,c("COSMIC identifier" = "sample_id")) %>% 
  dplyr::select(cell_line = `Sample Name`,
                tcga_code = `Cancer Type\n(matching TCGA label)`,
                gender
                )
gdsc_dat <- dplyr::inner_join(gdsc_anno_dat,gdsc_drug_dat) %>% 
  dplyr::select({names(ccle_dat)})

gdsc_dat$gender <- mapvalues(gdsc_dat$gender,c("f","m","u",NA),
                             c("FEMALE","MALE","Unknown","Unknown"))




drug_res_dat <- rbind(ccle_dat,gdsc_dat) %>% 
  dplyr::filter(!(tcga_code %in% c("","UNABLE TO CLASSIFY")) , !is.na(tcga_code))
       
saveRDS(drug_res_dat,"result/drug/drug_response.rds")


