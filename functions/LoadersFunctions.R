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
  FamilyCase <- read_qiime_sum_feats (FamCaseOtuTbl)  
  GenusCase <- read_qiime_sum_feats (GenCaseOtuTbl) 
  SpeciesCase <- read_qiime_sum_feats (SpeCaseOtuTbl)   
  OtuCase <- read_qiime_otu_table_no_tax (OtuCaseTbl) 
  
  FamilyCtrl <- read_qiime_sum_feats (FamCtrlOtuTbl)  
  GenusCtrl <- read_qiime_sum_feats (GenCtrlOtuTbl) 
  SpeciesCtrl <- read_qiime_sum_feats (SpeCtrlOtuTbl)   
  OtuCtrl <- read_qiime_otu_table_no_tax (OtuCtrlTbl) 
  
  # percent OTU
  OtupCase <- 100 * OtuCase / rowSums(OtuCase)
  OtupCtrl <- 100 * OtuCtrl / rowSums(OtuCtrl)
  
  list(FamilyCase, GenusCase, SpeciesCase, OtuCase, OtupCase ,FamilyCtrl, GenusCtrl, SpeciesCtrl, OtuCtrl, OtupCtrl)
 }








