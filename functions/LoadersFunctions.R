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