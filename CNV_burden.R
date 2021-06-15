library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)


# copy number segments 
dat <- fread("rawdata/pancan/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.xena.gz")
anno_dat <- fread("rawdata/pancan/pancan_Curated clinical data") %>% 
  dplyr::filter(as.numeric(substr(sample,14,15)) < 10)   # remove no-tumor sample
#----------------------------------------CNA score-----------------------------------------

## 根据文献，矫正和计算 burden DOI: https://doi.org/10.7554/eLife.50267

pr_dat <- readRDS("result/purity.rds") %>% 
  tibble::rownames_to_column(var = "sampleID") %>% 
  dplyr::select(sampleID,TumorPurity)

dat_c <- dplyr::inner_join(dat,pr_dat,by = "sampleID")

r.lim <- 0.4 # min purity accepted
loss.lim <- log2((1/2)*r.lim) # min loss value accepted

dat_c$seg.mean <- map2(dat_c$value,dat_c$TumorPurity,function(s,p){
  if (p < r.lim) {
    p <- r.lim
  }
  inside.log <- (2^s + (p - 1)) / p
  if (inside.log <= 2^(loss.lim)) {inside.log <- 2^(loss.lim)}
  log2(inside.log)
}) %>% unlist()


######
l3 <- read.table("rawdata/cytobands_level3_pq.csv", header = TRUE, sep = "\t")
l4 <- read.table("rawdata/cytobands_level4_chrom.csv", header = TRUE, sep = "\t")

l3$length <- l3$end - l3$start
l4$length <- l4$end - l4$start

filt <- dat_c %>% 
  rename("loc.end" = "end","loc.start" = "start") %>% 
  as.data.frame()

filt$chr <- gsub("chr", "", as.character(filt$chr))

filt$chr[which(filt$chr == "X")] <- 23
filt$chr[which(filt$chr == "Y")] <- 24
filt$chr <- as.numeric(filt$chr)
filt$length <- filt$loc.end - filt$loc.start
filt$BAF <- 0.5

filt$score <- NA
filt$classified <- NA
filt$type <- NA
filt$intensity <- NA
filt$weight <- NA
filt$comments <- NA
chrom.percent <- 0.9
arm.percent <- 0.5

min.length <- 100000
max.dist.segm <- 1000000

low.cutoff.up <- 0.2
medium.cutoff.up <- round(log2(3 / 2), 2)
high.cutoff.up <- round(log2(4 / 2), 2)

low.cutoff.dw <- -0.2
medium.cutoff.dw <- round(log2(1 / 2), 2)
high.cutoff.dw <- round(log2(0.6 / 2), 2)

chrom.percent <- 0.9
arm.percent <- 0.5

focal.percent.low <- 0.05
focal.percent.medium <- 0.15
focal.percent.high <- 0.3

min.baf <- 0.2

for (i in 1:nrow(filt)) {
  
  #### remove "chr"
  chr <- gsub("chr", "", as.character(filt[i, "chr"]))
  
  
  ### check if it is a chrom-level SCNA
  if (filt[i, "length"] > chrom.percent * l4[which(l4$chr == chr), "length"]) {
    w <- 0
    
    ## classify gains
    if (filt[i, "seg.mean"] >= low.cutoff.up & filt[i, "seg.mean"] < medium.cutoff.up) w <- 1
    if (filt[i, "seg.mean"] >= medium.cutoff.up & filt[i, "seg.mean"] < high.cutoff.up) w <- 2
    if (filt[i, "seg.mean"] >= high.cutoff.up) w <- 3
    
    ## classify losses
    if (filt[i, "seg.mean"] <= (low.cutoff.dw) & filt[i, "seg.mean"] > (medium.cutoff.dw)) w <- 1
    if (filt[i, "seg.mean"] <= (medium.cutoff.dw) & filt[i, "seg.mean"] > (high.cutoff.dw)) w <- 2
    if (filt[i, "seg.mean"] <= (high.cutoff.dw)) w <- 3
    
    if (!is.na(filt[i, "BAF"])) {
      ## classify upd
      if (filt[i, "seg.mean"] >= (low.cutoff.dw) & filt[i, "seg.mean"] <= (low.cutoff.up) & (filt[i, "BAF"] > (0.5 + min.baf) | filt[i, "BAF"] < (0.5 - min.baf))) {
        w <- 2
        filt$type[i] <- "CN-LOH"
      }
    }
    
    
    
    if (w != 0) {
      filt$classified[i] <- "chromosomal"
      if (is.na(filt$type[i])) {
        if (filt[i, "seg.mean"] > 0) {
          filt$type[i] <- "Gain"
        } else {
          filt$type[i] <- "Loss"
        }
      }
      filt$score[i] <- w
      filt$intensity[i] <- w
    }
    
    # skip to new cnv
    next
  }
  
  
  ### check if it is an arm-level SCNA
  arms <- l3[which(l3$chr == chr), ]
  centromer <- arms[1, "end"]
  l.p <- (min(c(centromer, filt[i, "loc.end"])) - filt[i, "loc.start"]) / arms[1, "length"]
  l.q <- (filt[i, "loc.end"] - (max(c(centromer, filt[i, "loc.start"])))) / arms[2, "length"]
  
  if (l.p > arm.percent & l.q > arm.percent) {
    print("WARNING: both arms significant!")
    filt$comments[i] <- "Both arms significant!" # but only one is counted
  }
  
  if (l.p > arm.percent) {
    w <- 0
    
    ## classify gains
    if (filt[i, "seg.mean"] >= low.cutoff.up & filt[i, "seg.mean"] < medium.cutoff.up) w <- 1
    if (filt[i, "seg.mean"] >= medium.cutoff.up & filt[i, "seg.mean"] < high.cutoff.up) w <- 2
    if (filt[i, "seg.mean"] >= high.cutoff.up) w <- 3
    
    ## classify losses
    if (filt[i, "seg.mean"] <= (low.cutoff.dw) & filt[i, "seg.mean"] > (medium.cutoff.dw)) w <- 1
    if (filt[i, "seg.mean"] <= (medium.cutoff.dw) & filt[i, "seg.mean"] > (high.cutoff.dw)) w <- 2
    if (filt[i, "seg.mean"] <= (high.cutoff.dw)) w <- 3
    
    if (!is.na(filt[i, "BAF"])) {
      ## classify upd
      if (filt[i, "seg.mean"] >= (low.cutoff.dw) & filt[i, "seg.mean"] <= (low.cutoff.up) & (filt[i, "BAF"] > (0.5 + min.baf) | filt[i, "BAF"] < (0.5 - min.baf))) {
        w <- 2
        filt$type[i] <- "CN-LOH"
      }
    }
    
    if (w != 0) {
      filt$classified[i] <- "arm"
      if (is.na(filt$type[i])) {
        if (filt[i, "seg.mean"] > 0) {
          filt$type[i] <- "Gain"
        } else {
          filt$type[i] <- "Loss"
        }
      }
      filt$score[i] <- w
      filt$intensity[i] <- w
    }
    
    # skip to new CNV
    next
  }
  
  if (l.q > arm.percent) {
    w <- 0
    
    ## classify gains
    if (filt[i, "seg.mean"] >= low.cutoff.up & filt[i, "seg.mean"] < medium.cutoff.up) w <- 1
    if (filt[i, "seg.mean"] >= medium.cutoff.up & filt[i, "seg.mean"] < high.cutoff.up) w <- 2
    if (filt[i, "seg.mean"] >= high.cutoff.up) w <- 3
    
    ## classify losses
    if (filt[i, "seg.mean"] <= (low.cutoff.dw) & filt[i, "seg.mean"] > (medium.cutoff.dw)) w <- 1
    if (filt[i, "seg.mean"] <= (medium.cutoff.dw) & filt[i, "seg.mean"] > (high.cutoff.dw)) w <- 2
    if (filt[i, "seg.mean"] <= (high.cutoff.dw)) w <- 3
    
    ## classify upd
    if (!is.na(filt[i, "BAF"])) {
      if (filt[i, "seg.mean"] >= (low.cutoff.dw) & filt[i, "seg.mean"] <= (low.cutoff.up) & (filt[i, "BAF"] > (0.5 + min.baf) | filt[i, "BAF"] < (0.5 - min.baf))) {
        w <- 2
        filt$type[i] <- "CN-LOH"
      }
    }
    
    if (w != 0) {
      filt$classified[i] <- "arm"
      if (is.na(filt$type[i])) {
        if (filt[i, "seg.mean"] > 0) {
          filt$type[i] <- "Gain"
        } else {
          filt$type[i] <- "Loss"
        }
      }
      filt$score[i] <- w
      filt$intensity[i] <- w
    }
    # skip to new CNV
    next
  }
  
  
  ### check if it is a focal-level SCNA
  if (filt[i, "seg.mean"] > low.cutoff.up | filt[i, "seg.mean"] < low.cutoff.dw) {
    ww <- 0
    
    
    ### identify percentages of arms
    arms <- l3[which(l3$chr == chr), ]
    centromer <- arms[1, "end"]
    l.p <- (min(c(centromer, filt[i, "loc.end"])) - filt[i, "loc.start"]) / arms[1, "length"]
    l.q <- (filt[i, "loc.end"] - (max(c(centromer, filt[i, "loc.start"])))) / arms[2, "length"]
    l.tot <- max(l.p, l.q)
    
    if (l.p > 0 & l.q > 0) {
      print("WARNING: focal SCNA including centromer!")
      filt$comments[i] <- "Focal SCNA including centromer" # only the highest weight is counted
    }
    
    
    ## weight of focal SCNA according to its length
    if (l.tot <= focal.percent.low) ww <- 1
    if (l.tot <= focal.percent.medium & l.tot > focal.percent.low) ww <- 2
    if (l.tot <= focal.percent.high & l.tot > focal.percent.medium) ww <- 3
    if (l.tot > focal.percent.high) ww <- 4
    
    w <- 0
    
    ## classify gains
    if (filt[i, "seg.mean"] >= low.cutoff.up & filt[i, "seg.mean"] < medium.cutoff.up) w <- 1
    if (filt[i, "seg.mean"] >= medium.cutoff.up & filt[i, "seg.mean"] < high.cutoff.up) w <- 2
    if (filt[i, "seg.mean"] >= high.cutoff.up) w <- 3
    
    ## classify losses
    if (filt[i, "seg.mean"] <= (low.cutoff.dw) & filt[i, "seg.mean"] > (medium.cutoff.dw)) w <- 1
    if (filt[i, "seg.mean"] <= (medium.cutoff.dw) & filt[i, "seg.mean"] > (high.cutoff.dw)) w <- 2
    if (filt[i, "seg.mean"] <= (high.cutoff.dw)) w <- 3
    
    
    if (w != 0) {
      filt$classified[i] <- "focal"
      if (filt[i, "seg.mean"] > 0) {
        filt$type[i] <- "Gain"
      } else {
        filt$type[i] <- "Loss"
      }
      filt$score[i] <- w * ww
      filt$intensity[i] <- w
      filt$weight[i] <- ww
    }
    
    
    # skip to new cnv
    next
  }
  
  
  
  ## classify CN-LOH
  if (!is.na(filt[i, "BAF"])) {
    
    
    ### identify percentages of arms
    arms <- l3[which(l3$chr == chr), ]
    centromer <- arms[1, "end"]
    l.p <- (min(c(centromer, filt[i, "loc.end"])) - filt[i, "loc.start"]) / arms[1, "length"]
    l.q <- (filt[i, "loc.end"] - (max(c(centromer, filt[i, "loc.start"])))) / arms[2, "length"]
    l.tot <- max(l.p, l.q)
    
    if (l.p > 0 & l.q > 0) {
      print("WARNING: focal SCNA including centromer!")
      filt$comments[i] <- "Focal SCNA including centromer" # only the highest weight is counted
    }
    
    
    ww <- 0
    
    ## weight of focal SCNA according to its length
    if (l.tot <= focal.percent.low) ww <- 1
    if (l.tot <= focal.percent.medium & l.tot > focal.percent.low) ww <- 2
    if (l.tot <= focal.percent.high & l.tot > focal.percent.medium) ww <- 3
    if (l.tot > focal.percent.high) ww <- 4
    
    w <- 0
    
    if (filt[i, "seg.mean"] >= (low.cutoff.dw) & filt[i, "seg.mean"] <= (low.cutoff.up) & (filt[i, "BAF"] > (0.5 + min.baf) | filt[i, "BAF"] < (0.5 - min.baf))) {
      w <- 2
    }
    
    
    if (w != 0) {
      filt$classified[i] <- "focal"
      filt$type[i] <- "CN-LOH"
      filt$score[i] <- w * ww
      filt$intensity[i] <- w
      filt$weight[i] <- ww
    }
    
    # skip to new cnv
    next
  }
}


score_res <- filter(filt,!is.na(score)) %>% 
  group_by(sampleID,classified) %>% 
  summarise(score  = sum(score)) %>% 
  tidyr::pivot_wider(names_from = classified, 
                     values_from = `score`,
                     values_fill = 0) %>% 
  mutate(BCS = arm + chromosomal, FCS = focal)

normBCS <- with(score_res, (BCS - mean(BCS))/sd(BCS))
normFCS <- with(score_res, (FCS - mean(FCS))/sd(FCS))

score_res$GCS <- normBCS + normFCS
burden_res_anno <- dplyr::inner_join(score_res,anno_dat,by=c("sampleID" = "sample"))


saveRDS(burden_res_anno[,1:11],"result/diff/cnv_burden.rds")

ggplot(burden_res_anno,aes(x = `cancer type abbreviation` , y = BCS)) +  # 验证成功
  geom_boxplot()+
  geom_signif()


test_dat <- burden_res_anno[burden_res_anno$`cancer type abbreviation` == "LUAD",]

wilcox.test(FCS ~ gender , data=test_dat)$p.value


#----------------------------------Region profile -------------------------------------------------

hg19_1m <- read.table("rawdata/CNApp/autosomes_hg19_by_1Mb.txt")


mean_region_fun <- function(chr,reg_start,reg_end,segdat){
  
  dat <- segdat[segdat$chr == chr,]
  per_v <- c()
  reg.mean_v <- c()
  
  case1 <- which( # ..*|*|..
    dat$start < reg_start & (dat$end > reg_start & dat$end <= reg_end)
  )

  if(length(case1)>=1){
    for(i in case1){
      per <- (dat[i,"end"] - reg_start)/(reg_end - reg_start)
      reg.mean <- dat[i,"seg.mean"]
      per_v <- c(per_v,per)
      reg.mean_v <- c(reg.mean_v , reg.mean)
    }
  }

  case2 <- which( # ..|**|..
    dat$start >= reg_start & dat$end <= reg_end
  )

  if(length(case2)>=1){
    for(i in case2){
      per <- (dat[i,"end"] - dat[i,"start"])/(reg_end - reg_start)
      reg.mean <- dat[i,"seg.mean"]
      per_v <- c(per_v,per)
      reg.mean_v <- c(reg.mean_v , reg.mean)
    }
  }
  
  case3 <- which( # ..|*|*..
    (dat$start >= reg_start & dat$start < reg_end) & dat$end > reg_end
  )

  if(length(case3)>=1){
    for(i in case3){
      per <- (reg_end - dat[i,"start"])/(reg_end - reg_start)
      reg.mean <- dat[i,"seg.mean"]
      per_v <- c(per_v,per)
      reg.mean_v <- c(reg.mean_v , reg.mean)
    }
  }

  case4 <- which( # ..*||*..
    dat$start < reg_start & dat$end > reg_end
  )

  if(length(case4)>=1){
    for(i in case4){
      per <- 1
      reg.mean <- dat[i,"seg.mean"]
      per_v <- c(per_v,per)
      reg.mean_v <- c(reg.mean_v , reg.mean)
    }
  }

  #if(sum(per_v)>1){print("erro")}
  
  if(length(per_v)>=1){
    mean_region <-  sum(reg.mean_v * per_v)/sum(per_v)
  }else{
    mean_region <- 0 
  }
  
  mean_region
}

##
library(parallel)

seg_file <- dat %>% rename("seg.mean" = "value") %>% 
  dplyr::mutate(chr = paste0("chr",chr)) %>% 
  dplyr::filter(chr %in% paste0("chr",1:22)) %>% 
  as.data.frame()

cl <- makeCluster(15, type="FORK")

seg_mean_res <- parLapply(cl,unique(seg_file$sampleID), function(x){
  print(x)
  segdat <- seg_file[seg_file$sampleID == x,]
  seg_mean <- pmap(list(chr = hg19_1m[,1],reg_start = hg19_1m[,2],reg_end = hg19_1m[,3]),
       mean_region_fun,segdat = segdat) %>% unlist()
  seg_mean
}) 

stopCluster(cl)

smr <- do.call(rbind,seg_mean_res) %>% t() %>% as.data.frame()
colnames(smr) <- unique(seg_file$sampleID)
smr_anno <- anno_dat[anno_dat$sample %in% colnames(smr) ,] %>% 
  dplyr::filter(gender %in% c("FEMALE","MALE")) %>% 
  dplyr::select(sample,`cancer type abbreviation`,gender)
smr <- smr[,match(smr_anno$sample,colnames(smr))]

sl <- list(seg_mean_res = smr,smr_anno = smr_anno,chr_anno = hg19_1m[,1])
saveRDS(object = sl,file = "result/diff/cnv_region_profile.rds")




  