#####################################################
## configuration data                             ###
#####################################################

#--- case ---
DataCase <- "/home/anna/metagenome/ParkinsonAndSklerosis/data/qiime/Case"

# output case directory 
OutdirCase <- "/home/anna/metagenome/ParkinsonAndSklerosis/out/Case"

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

#fam_case_otu_table <- (paste(data_case, fam_case_file, sep = "/")) 
#g_case_otu_table <- (paste(data_case, g_case_file, sep = "/")) 
#sp_case_otu_table <- (paste(data_case, sp_case_file, sep = "/")) 
#otu_case_otu_table <- (paste(data_case, otu_case_file, sep = "/")) 

# case alpha diversity
alpha_div_case <- "data/qiime/Case/alpha_collated/chao1.txt"


#---control---
DataControl <- "/home/anna/metagenome/ParkinsonAndSklerosis/data/qiime/Control"

#output control directory 

OutdirControl <- "/home/anna/metagenome/ParkinsonAndSklerosis/out"

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

#fam_ctrl_otu_table <- (paste(data_case, fam_ctrl_file, sep = "/")) 
#g_ctrl_otu_table <- (paste(data_case, g_ctrl_file, sep = "/")) 
#sp_ctrl_otu_table <- (paste(data_case, sp_ctrl_file, sep = "/")) 
#otu_ctrl_otu_table <- (paste(data_case, otu_ctrl_file, sep = "/")) 

# control alpha diversity
alpha_div_ctrl <- "data/qiime/Control/alpha_collated/chao1_MOD.txt"

#---sequencing statistics---
SeqStatTblFile <- "/home/anna/metagenome/ParkinsonAndSklerosis/StatTable.txt"



