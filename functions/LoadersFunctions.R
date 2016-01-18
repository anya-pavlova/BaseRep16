########################################################
## Loading functions                                  ##
########################################################


#---LoadAlphaDiv: loading alpha diversity---
LoadAlphaDiv <- function(inpdir)
{
  inpdir
  AlphaDivTbl <- (read.table (inpdir,row.names=1,  header=T,sep="\t",stringsAsFactors=F ))
  AlphaDivTbl <- t(AlphaDivTbl)
  AlphaDivTbl <- as.data.frame(AlphaDivTbl[-c(1, 2),])
  rownames(AlphaDivTbl)<-gsub("X","", rownames(AlphaDivTbl))
  colnames(AlphaDivTbl) <- c("AlphaDiversity")
  return(AlphaDivTbl)
}
#---end LoadAlphaDiv---

#---Loading case and control---
Load <- function(FamInpCase, GenInpCase, SpeInpCase, OtuInpCase, FamInpCtrl, GenInpCtrl, SpeInpCtrl, OtuInpCtrl)
  #Load <- function(FamInpCase, FamInpCtrl, GenInpCase, GenInpCtrl, SpeInpCase, SpeInpCtrl, OtuInpCase, OtuInpCtrl)
{
  FamilyCase <- read_qiime_sum_feats (FamInpCase)  
  GenusCase <- read_qiime_sum_feats (GenInpCase) 
  SpeciesCase <- read_qiime_sum_feats (SpeInpCase)   
  OtuCase <- read_qiime_otu_table_no_tax (OtuInpCase) 
  
  FamilyCtrl <- read_qiime_sum_feats (FamInpCtrl)  
  GenusCtrl <- read_qiime_sum_feats (GenInpCtrl) 
  SpeciesCtrl <- read_qiime_sum_feats (SpeInpCtrl)   
  OtuCtrl <- read_qiime_otu_table_no_tax (OtuInpCtrl) 
  
  # percent OTU
  OtupCase <- 100 * OtuCase / rowSums(OtuCase)
  OtupCtrl <- 100 * OtuCtrl / rowSums(OtuCtrl)
  
  list(FamilyCase, GenusCase, SpeciesCase, OtuCase, OtupCase ,FamilyCtrl, GenusCtrl, SpeciesCtrl, OtuCtrl, OtupCtrl)
}

