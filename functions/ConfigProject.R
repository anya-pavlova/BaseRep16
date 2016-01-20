#####################################################
## configuration data                             ###
#####################################################

#--- general graphs directory ---

GraphsDir <- "/home/anna/metagenome/HSN/Graphs"

#--- case ---
DataCase <- "/home/anna/metagenome/HSN/data/qiime/Case"

# output case directory 
OutdirCase <- "/home/anna/metagenome/HSN/out/Case"

# case fam, g, sp, otu files
FamCaseFile <- "otu_table_L5.txt"
GenCaseFile <- "otu_table_L6.txt"
SpeCaseFile <- "otu_table_L7.txt"
OtuCaseFile <- "otu_table_MOD.txt"

# case fam, g, sp, otu tablesn
FamCaseOtuTbl <- (paste(DataCase, FamCaseFile, sep = "/")) 
GenCaseOtuTbl <- (paste(DataCase, GenCaseFile, sep = "/")) 
SpeCaseOtuTbl <- (paste(DataCase, SpeCaseFile, sep = "/")) 
OtuCaseTbl <- (paste(DataCase, OtuCaseFile, sep = "/")) 

# case alpha diversity
AlphaDivCase <- "data/qiime/Case/alpha_collated/chao1.txt"

#case meta data
MetaCaseCsv <- "/home/anna/metagenome/HSN/MetaDataCase.csv"

#---control---
DataControl <- "/home/anna/metagenome/HSN/data/qiime/Control"

#output control directory 

OutdirControl <- "/home/anna/metagenome/HSN/out/Control"

# control fam, g, sp, otu files
FamCtrlFile <- "otu_table_L5.txt"
GenCtrlFile <- "otu_table_L6.txt"
SpeCtrlFile <- "otu_table_L7.txt"
OtuCtrlFile <- "otu_table_MOD.txt"

# control fam, g, sp, otu tables
FamCtrlOtuTbl <- (paste(DataControl, FamCtrlFile, sep = "/")) 
GenCtrlOtuTbl <- (paste(DataControl, GenCtrlFile, sep = "/")) 
SpeCtrlOtuTbl <- (paste(DataControl, SpeCtrlFile, sep = "/")) 
OtuCtrlTbl <- (paste(DataControl, OtuCtrlFile, sep = "/")) 

# control alpha diversity
AlphaDivCtrl <- "data/qiime/Case/alpha_collated/chao1.txt"

#control meta data
MetaCtrlCsv <- "/home/anna/metagenome/HSN/MetaDataCtrl.csv"

#---sequencing statistics---
SeqStatTblFile <- "/home/anna/metagenome/HSN/StatTable.txt"

ReadIni <- function(IniFilename) 
{ 
  connection <- file(IniFilename) 
  Lines  <- readLines(connection) 
  close(connection) 
  
  Lines <- chartr("[]", "==", Lines)  # change section headers 
  
  connection <- textConnection(Lines) 
  d <- read.table(connection, as.is = TRUE, sep = "=", fill = TRUE) 
  close(connection) 
  
  L <- d$V1 == ""                    # location of section breaks 
  d <- subset(transform(d, V3 = V2[which(L)[cumsum(L)]])[1:3], 
              V1 != "") 
  
  ToParse  <- paste("INI.list$", d$V3, "$",  d$V1, " <- '", 
                    d$V2, "'", sep="") 
  
  INI.list <- list() 
  eval(parse(text=ToParse)) 
  
  return(INI.list) 
} 


PathWay <- ReadIni("/home/anna/metagenome/HSN/Patway/Pathway.ini") 

Family <- read_qiime_sum_feats (PathWay$Case$FamCaseOtuTbl)
