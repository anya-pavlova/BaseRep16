########################################################
## Loading functions                                  ##
########################################################
library(futile.logger)

# Loading case and control
Load <- function(pathway)
{
  flog.info("start load table")
  family <- UniteMatrices(read_qiime_sum_feats(pathway$Case$FamCaseOtuTbl), read_qiime_sum_feats (pathway$Ctrl$FamCtrlOtuTbl))
  genus <- UniteMatrices(read_qiime_sum_feats (pathway$Case$GenCaseOtuTbl), read_qiime_sum_feats (pathway$Ctrl$GenCtrlOtuTbl))
  species <- UniteMatrices(read_qiime_sum_feats (pathway$Case$SpeCaseOtuTbl), read_qiime_sum_feats (pathway$Ctrl$SpeCtrlOtuTbl))
  otu <- UniteMatrices(read_qiime_otu_table_no_tax (pathway$Case$OtuCaseTbl), read_qiime_otu_table_no_tax (pathway$Ctrl$OtuCtrlTbl))
  alphaDiv <- UniteMatrices(LoadAlphaDiv(pathway$Case$AlphaDivCase), LoadAlphaDiv(pathway$Ctrl$AlphaDivCtrl))
  
  # percent OTU
  otup <- 100 * otu / rowSums(otu)
  
  ######### insert load alpha diversity as vector in list
  # load meta data
  metaTable <- rbind(read.csv(pathway$Case$MetaCaseCsv, stringsAsFactors=F), read.csv(pathway$Ctrl$MetaCtrlCsv, stringsAsFactors=F))

  flog.info("finished load table")
  #TODO(Anna) create checking ...
  list(family=family, genus=genus, species=species, otu=otu, otup=otup, meta=metaTable, alphaDiv=alphaDiv)
}
#--- end loading case and control ---






##### check rowSums=100%
#rowSums(TotalTable$family)

##### check rownames matching for all tables and metadata
#setdiff(rownames(FamilyCtrl), MetaCtrl$samples_name)
#setdiff(rownames(GenusCtrl), MetaCtrl$samples_name)
#setdiff(rownames(SpeciesCtrl), MetaCtrl$samples_name)

#rownames(TotalTable$family)



