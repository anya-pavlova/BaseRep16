#####################################################
## configuration data                             ###
#####################################################

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
alpha_div_case <- "data/qiime/Case/alpha_collated/chao1.txt"

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

# control alpha diversity
alpha_div_ctrl <- "data/qiime/Control/alpha_collated/chao1_MOD.txt"

#---sequencing statistics---
SeqStatTblFile <- "/home/anna/metagenome/HSN/StatTable.txt"



