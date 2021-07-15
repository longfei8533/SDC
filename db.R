# MRNA CPM and miRNA CPM data are written to the database for easy retrieval, 
# rather than being read directly into memory. 

library(DBI)
library(dbplyr)
library(data.table)
library(dplyr)
library(ggstatsplot)
library(edgeR)

cli_dat <- fread("rawdata/pancan/pancan_Curated clinical data")
# -------------------------mRNA----------------------------

mRNA_count <- fread("rawdata/toil/TcgaTargetGtex_gene_expected_count.gz")

con <- dbConnect(RSQLite::SQLite(), "result/RNA_cpm.db")
for(i in unique(cli_dat$`cancer type abbreviation`)){
  ab <- cli_dat$sample[cli_dat$`cancer type abbreviation` == i]
  ab_loc <- which(colnames(mRNA_count) %in% ab)
  ab_dat <- mRNA_count[,..ab_loc]
  ab_dat <- 2^ab_dat - 1
  ab_cpm <- edgeR::cpm(ab_dat,normalized.lib.sizes = TRUE,log = TRUE)
  m_dat <- data.frame(RNA  = mRNA_count$sample,ab_cpm,check.names = F)
  dbWriteTable(con,i, m_dat)
}
dbDisconnect(con)
# ------------------------miRAN----------------
con <- dbConnect(RSQLite::SQLite(), "result/miRNA_cpm.db")
for(i in unique(cli_dat$`cancer type abbreviation`)){
  file_name <- paste0("TCGA-",i,"_miRNA_gene_quantification.rda")
  load(paste0("rawdata/pancan/miRNA/",file_name))
  dat <- as.data.frame(data)
  mirname <- dat$miRNA_ID
  dat <- dplyr::select(dat, starts_with("read_count"))
  colnames(dat) <- substr(colnames(dat), 12, 26)
  dat <- dat[,!duplicated(colnames(dat))]
  dat_cpm <- edgeR::cpm(dat,normalized.lib.sizes = TRUE,log = TRUE)
  mi_dat <- data.frame(mir = mirname, dat_cpm,check.names = F)
  dbWriteTable(con,i,mi_dat)
}
dbDisconnect(con)

# -------------------------- test -------------------------

mcon <- dbConnect(RSQLite::SQLite(), "result/RNA_cpm.db")
micon <- dbConnect(RSQLite::SQLite(), "result/miRNA_cpm.db")


dqm <-  paste0("SELECT * FROM `STAD` WHERE (`RNA` LIKE '%","ENSG00000167323","%')")
dqmi <-  paste0("SELECT * FROM `STAD` WHERE (`mir` = '","hsa-mir-372","')")

m_res <- dbSendQuery(mcon,dqm)
mi_res <- dbSendQuery(micon,dqmi)

m <- dbFetch(m_res) %>% t() %>%
  .[-1, ] %>%
  as.data.frame()
mi <- dbFetch(mi_res) %>%  t() %>%
  .[-1, ] %>%
  as.data.frame()

dbClearResult(m_res)
dbClearResult(mi_res)

dat <- merge(m, mi, by = "row.names", all = FALSE)
colnames(dat) <- c("Sample","mRNA", "miRNA")
dat$mRNA <- as.character(dat$mRNA) %>% as.numeric()
dat$miRNA <- as.character(dat$miRNA) %>% as.numeric()
dat

dbClearResult(m_res)
dbClearResult(mi_res)

p <- ggscatterstats(
    data = dat,
    x = miRNA,
    y = mRNA,
    xlab = "x",
    ylab = "y",
    messages = FALSE)





