# miRNA cor with mRNA
library(dplyr)
library(ppcor)
library(purrr)

# miRNA
load("rawdata/miRNA/TCGA-ACC_miRNA_gene_quantification.rda")

mi_dat <- as.data.frame(data)
mirname <- mi_dat$miRNA_ID
mi_dat <- dplyr::select(mi_dat, starts_with("read_count"))
colnames(mi_dat) <- substr(colnames(mi_dat), 12, 26)
mi_dat <- mi_dat[,as.numeric(substr(colnames(mi_dat),14,15)) < 10]

# mRNA

mRNA_dat <- fread("rawdata/toil/TcgaTargetGtex_gene_expected_count.gz",
                   drop = "sample",
                   nThread = 15) 
mRNA_dat <- 2^mRNA_dat - 1 
mRNA_dat <- as.data.frame(mRNA_dat)
gene <- fread("rawdata/toil/TcgaTargetGtex_gene_expected_count.gz",
              select = "sample")
gene <- gene$sample %>% substr(.,1,15)

# purity

purity_data <- readRDS("result/purity.rds")

# com

com <- base::Reduce(base::intersect,
                    list(colnames(mi_dat),colnames(mRNA_dat),rownames(purity_data)))

#
mRNA <- mRNA_dat[,com]
mi <- mi_dat[,com]

pu <- purity_data[com,"TumorPurity"]

tm <- t(mRNA) %>% as.data.frame()
ti <- t(mi)
x <- ti[,1]

library(parallel)

cl <- makeCluster(20)
xx <- parLapply(cl,tm,pcor.test,y=x,z=pu)
stopCluster(cl)




