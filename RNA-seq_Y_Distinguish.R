library(data.table)
library(tidyverse)
library(DESeq2)
setwd("/media/CADD/longfei/project/sex_web/")

phenotype_dat <- fread("rawdata/toil/TcgaTargetGTEX_phenotype.txt.gz") %>% 
  dplyr::filter(`_gender` %in% c("Female","Male"),`_study` %in% c("GTEX","TCGA")) %>% 
  as.data.frame()

count_dat <- fread("rawdata/toil/TcgaTargetGtex_rsem_gene_tpm.gz",
                   drop = "sample",
                   nThread = 15) 


gene <- fread("rawdata/toil/TcgaTargetGtex_gene_expected_count.gz",
              select = "sample")
gene <- gene$sample %>% substr(.,1,15)
library(AnnoProbe)
ids_chr <- AnnoProbe::annoGene(gene,ID_type = "ENSEMBL")
ids_sel <- dplyr::filter(ids_chr,chr %in% c("chrY"))


sam <- phenotype_dat[phenotype_dat$sample %in% colnames(count_dat),]

x <- match(sam$sample,colnames(count_dat))

countData <- count_dat[,..x]
countData_Y <- countData[gene %in% ids_sel$ENSEMBL,]


library(factoextra)
pca_dat <- prcomp(x=t(countData_Y),rank.=10,scale. = F)
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
                     habillage = sam$`_gender`,
                     palette = c("red","blue"),
                     alpha.ind=0.5,
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




png(filename = "pca_GTEX_TCGA.png",width = 20,height = 14,units = "cm",res = 300)
pca1
dev.off()
