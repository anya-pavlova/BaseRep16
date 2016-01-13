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
  colnames(AlphaDivTbl) <- c("Alpha Diversity")
  return(AlphaDivTbl)
}
#---end LoadAlphaDiv---