library(IOBR)
IOBR_signature <- list("TME-associated" = names(signature_tme),
                       "Tumor-metabolism" = names(signature_metabolism), 
                       "Tumor-intrinsic"= names(signature_tumor)
                         )
IOBR_signature_gene <- signature_collection
IOBR_signature_citation <- signature_collection_citation


save(IOBR_signature,IOBR_signature_gene,IOBR_signature_citation,sig_group,
     file = "process/result/imm/IOBR_signature.RData")
