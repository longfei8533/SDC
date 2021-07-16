library(data.table)
library(tidyverse)

cell_anno_dat <- fread("rawdata/GDSC/CosmicSample.tsv.gz") %>% 
  dplyr::filter(sample_type == "cell-line")
GDSC_cell_dat <- fread("rawdata/GDSC/Cell_Lines_Details.csv")

## 统计多少细胞系性别已知


match_dat <- cell_anno_dat[cell_anno_dat$sample_id %in% GDSC_cell_dat$`COSMIC identifier`,] %>% 
  dplyr::filter(gender %in% c("f","m"))

table(match_dat$gender)

# f   m   u 
# 264 338 112 

# 因为存在大面积性别缺失，所以通过 RNA-array 确定性别

library(affy)


## 注释基因名称和染色体位置
# library(devtools)
# remote::install_github("jmzeng1314/AnnoProbe")
library(AnnoProbe)
ids=idmap('GPL13667',type = 'soft')
ids_chr <- AnnoProbe::annoGene(ids$symbol,ID_type = "SYMBOL")
ids_sel <- dplyr::filter(ids_chr,chr %in% c("chrY")) %>% 
  dplyr::inner_join(ids,by=c("SYMBOL"="symbol"))


sample_anno <- fread("rawdata/GDSC/E-MTAB-3610.sdrf.txt")
dir_cels="rawdata/GDSC/cell_line_expression/"
affy_data = ReadAffy(celfile.path=dir_cels)
eset = rma(affy_data)

eset_dat <- as.data.frame(eset) 
eset_dat_Y <- eset_dat[,colnames(eset_dat) %in% paste0("X",ids_sel$ID)]

match_sample <- sample_anno[match(rownames(eset_dat_Y),sample_anno$`Array Data File`),]




sex <- ifelse(match_sample$`Characteristics[cell line]` %in% match_dat$sample_name[match_dat$gender == "f"],
              "F",ifelse(match_sample$`Characteristics[cell line]` %in% match_dat$sample_name[match_dat$gender == "m"],"M","U"))

## PCA 分析

library(factoextra)
pca_dat <- prcomp(x=eset_dat_Y,rank.=10,scale. = FALSE)
sc1 <- fviz_eig(pca_dat, addlabels = T,
                linecolor="#E64B35FF",
                barfill="#4DBBD5FF",
                barcolor="#4DBBD5FF",
                ggtheme=ggplot2::theme(
                  axis.title = element_text(size=14,face="bold"),
                  axis.text = element_text(face="bold",size=14,color = "black"),
                  axis.line = element_line(size=0.8,color="black"),
                  axis.ticks= element_line(size=0.8,colour = "black"),
                  panel.grid =element_blank(),
                  panel.background = element_blank(),
                  title = element_blank(),
                ),
                ylim=c(0,100)
)
  
pca1 <- fviz_pca_ind(pca_dat,geom=c("point"),
                     habillage = sex,
                     palette = c("red","blue","gray"),
                     #addEllipses = TRUE,
                     ggtheme=ggplot2::theme(
                       axis.title = element_text(size=12,face="bold"),
                       axis.text = element_text(face="bold",size=12,color = "black"),
                       axis.line = element_line(size=0.8,color="black"),
                       axis.ticks= element_line(size=0.8,colour = "black"),
                       panel.grid =element_blank(),
                       panel.background = element_blank(),
                       title = element_blank()
                     )
)

png("GDSC_pca.png",width = 20,height = 15,units = "cm",res = 300)
pca1
dev.off()