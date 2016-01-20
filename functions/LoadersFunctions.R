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
  return(data.matrix(AlphaDivTbl))
}
#---end LoadAlphaDiv---

PathWay <- ReadIni("/home/anna/metagenome/HSN/Patway/Pathway.ini") 
PathWayCase <- PathWay$Case
PathWayCtrl <- PathWay$Control
MetaTable <- read.csv(PathWayCase$MetaCaseCsv)

#---Loading case and control (family, genus, species,otu, meta data, alpha diversity)---
Load <- function (PatWay.ini)
{
  
  Family <- UniteMatrices(read_qiime_sum_feats (PathWay$Case$FamCaseOtuTbl), read_qiime_sum_feats (PathWay$Ctrl$FamCtrlOtuTbl))
  Genus <- UniteMatrices(read_qiime_sum_feats (PathWay$Case$GenCaseOtuTbl), read_qiime_sum_feats (PathWay$Ctrl$GenCtrlOtuTbl))
  Species <- UniteMatrices(read_qiime_sum_feats (PathWay$Case$SpeCaseOtuTbl), read_qiime_sum_feats (PathWay$Ctrl$SpeCtrlOtuTbl))
  Otu <- UniteMatrices(read_qiime_otu_table_no_tax (PathWay$Case$OtuCaseTbl), read_qiime_otu_table_no_tax (PathWay$Ctrl$OtuCtrlTbl))
  AlphaDiv <- UniteMatrices(LoadAlphaDiv(PathWay$Case$AlphaDivCase), LoadAlphaDiv(PathWay$Ctrl$AlphaDivCtrl))
  
  # percent OTU
  Otup <- 100 * Otu / rowSums(Otu)
  
  ######### insert load alpha diversity as vector in list
  # load meta data
  MetaTable <- rbind(read.csv(PathWay$Case$MetaCaseCsv), read.csv(PathWay$Ctrl$MetaCtrlCsv))
  
  list(Family=Family, Genus=Genus, Species=Species, Otu=Otu, Otup=Otup, Meta=MetaTable, AlphaDiv=AlphaDiv)
}


#--- loading case and control ---
# передавать вектор с путями в правильном (заданном) порядке в функцию, чтобы не ошибаться в порядке перечисления
TotalTable <- Load (PathWay)
#--- end loading case and control ---



###########################################################
### Loading case and control, alpha ad beta diversity  ###
###########################################################





##### check rowSums=100%
#rowSums(TotalTable$Family)

##### check rownames matching for all tables and metadata
#setdiff(rownames(FamilyCtrl), MetaCtrl$samples_name)
#setdiff(rownames(GenusCtrl), MetaCtrl$samples_name)
#setdiff(rownames(SpeciesCtrl), MetaCtrl$samples_name)

#rownames(TotalTable$Family)



