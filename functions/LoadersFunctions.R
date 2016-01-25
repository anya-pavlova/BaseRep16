########################################################
## Loading functions                                  ##
########################################################

#--- Load pathway ---
PathWay <- ReadIni("Patway/Pathway.ini") 
#--- end loading pathway ---
family <- as.data.frame(PathWay$Case$FamCaseOtuTbl)

#---Loading case and control (family, genus, species,otu, meta data, alpha diversity)---
Load <- function(PathWay)
{
  
  Family <- UniteMatrices(read_qiime_sum_feats(PathWay$Case$FamCaseOtuTbl), read_qiime_sum_feats (PathWay$Ctrl$FamCtrlOtuTbl))
  Genus <- UniteMatrices(read_qiime_sum_feats (PathWay$Case$GenCaseOtuTbl), read_qiime_sum_feats (PathWay$Ctrl$GenCtrlOtuTbl))
  Species <- UniteMatrices(read_qiime_sum_feats (PathWay$Case$SpeCaseOtuTbl), read_qiime_sum_feats (PathWay$Ctrl$SpeCtrlOtuTbl))
  Otu <- UniteMatrices(read_qiime_otu_table_no_tax (PathWay$Case$OtuCaseTbl), read_qiime_otu_table_no_tax (PathWay$Ctrl$OtuCtrlTbl))
  AlphaDiv <- UniteMatrices(LoadAlphaDiv(PathWay$Case$AlphaDivCase), LoadAlphaDiv(PathWay$Ctrl$AlphaDivCtrl))
  
  # percent OTU
  Otup <- 100 * Otu / rowSums(Otu)
  
  ######### insert load alpha diversity as vector in list
  # load meta data
  MetaTable <- rbind(read.csv(PathWay$Case$MetaCaseCsv, stringsAsFactors=F), read.csv(PathWay$Ctrl$MetaCtrlCsv, stringsAsFactors=F))

  #TODO(Anna) create checking ...
  
  list(Family=Family, Genus=Genus, Species=Species, Otu=Otu, Otup=Otup, Meta=MetaTable, AlphaDiv=AlphaDiv)
}
#--- end loading case and control (family, genus, species,otu, meta data, alpha diversity)---






##### check rowSums=100%
#rowSums(TotalTable$Family)

##### check rownames matching for all tables and metadata
#setdiff(rownames(FamilyCtrl), MetaCtrl$samples_name)
#setdiff(rownames(GenusCtrl), MetaCtrl$samples_name)
#setdiff(rownames(SpeciesCtrl), MetaCtrl$samples_name)

#rownames(TotalTable$Family)



