
cur <- "~/metagenome/HSN"
#cur <- ""
setwd(cur)


#---Congig---
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

#---control---
DataControl <- "/home/anna/metagenome/HSN/data/qiime/Control"

#output control directory 

OutdirControl <- "/home/anna/metagenome/HSN/out"

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


#--- Qiime functions ---
read_qiime_sum_feats <- function(filename)
{
  f <- t(read.table(filename, header=T, row.names=1, sep="\t"))
  rownames(f) <- gsub("X", "", rownames(f))
  f <- data.matrix(f)
  f <- 100 * f / rowSums(f)
  f
}

read_qiime_otu_table_no_tax <- function(filename)
{   
  f <- read.table(filename, header=T, row.names=1, sep="\t")
  # drop taxonomy
  f <- f[, -ncol(f)]
  f <- t(f)
  rownames(f) <- gsub("X", "", rownames(f))
  f
}

read_qiime_single_alpha_rar <- function(filename)
{   
  f <- read.table(filename, header=T, row.names=1, sep="\t", colClasses = "character")
  f <- f[,-c(1,2)]
  colnames(f) <- gsub("X", "", colnames(f))
  f <- t(f)
  colnames(f) <- "
  chao1"
  f <- f[which(f[,"chao1"] != "n/a"),,drop=F]
  f <- data.frame(f, stringsAsFactors = F)
  f[,"chao1"] <- as.numeric(f[,"chao1"])
  f
}


LoadCase <- function()
  #Load <- function(FamInpCase, FamInpCtrl, GenInpCase, GenInpCtrl, SpeInpCase, SpeInpCtrl, OtuInpCase, OtuInpCtrl)
{
  FamilyCase <- read_qiime_sum_feats ("/home/anna/metagenome/HSN/data/qiime/Case//otu_table_L5.txt")  
  GenusCase <- read_qiime_sum_feats ("/home/anna/metagenome/HSN/data/qiime/Case//otu_table_L6.txt") 
  SpeciesCase <- read_qiime_sum_feats ("/home/anna/metagenome/HSN/data/qiime/Case//otu_table_L7.txt")   
  OtuCase <- read_qiime_otu_table_no_tax ("/home/anna/metagenome/HSN/data/qiime/Case//otu_table_MOD.txt") 
  # percent OTU
  OtupCase <- 100 * OtuCase / rowSums(OtuCase)
  
  list(FamilyCase, GenusCase, SpeciesCase, OtuCase, OtupCase)
}

LoadCtrl <- function()
{
  FamilyCtrl <- read_qiime_sum_feats ("/home/anna/metagenome/HSN/data/qiime/Control/otu_table_L5.txt")  
  GenusCtrl <- read_qiime_sum_feats ("/home/anna/metagenome/HSN/data/qiime/Control/otu_table_L6.txt") 
  SpeciesCtrl <- read_qiime_sum_feats ("/home/anna/metagenome/HSN/data/qiime/Control/otu_table_L7.txt")   
  OtuCtrl <- read_qiime_otu_table_no_tax ("/home/anna/metagenome/HSN/data/qiime/Control/otu_table_MOD.txt") 
   
  OtupCtrl <- 100 * OtuCtrl / rowSums(OtuCtrl)
  list(FamilyCtrl, GenusCtrl, SpeciesCtrl, OtuCtrl, OtupCtrl)
}




