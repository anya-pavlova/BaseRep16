#####################################################
## configuration data                             ###
#####################################################

#---case---
data_case <- "/home/anna/metagenome/BaseRep16/data/qiime/Case"

#output case directory 
OutdirCase <- "/home/anna/metagenome/BaseRep16/out/Case"

# case fam, g, sp, otu files
fam_case_file <- "otu_table_L5.txt"
g_case_file <- "otu_table_L6.txt"
sp_case_file <- "otu_table_L7.txt"
otu_case_file <- "otu_table_MOD.txt"

# case fam, g, sp, otu tablesn
fam_case_otu_table <- (paste(data_case, fam_case_file, sep = "/")) 
g_case_otu_table <- (paste(data_case, g_case_file, sep = "/")) 
sp_case_otu_table <- (paste(data_case, sp_case_file, sep = "/")) 
otu_case_otu_table <- (paste(data_case, otu_case_file, sep = "/")) 

# case alpha diversity
alpha_div_case <- "data/qiime/Case/alpha_collated/chao1.txt"


#---control---
data_control <- "/home/anna/metagenome/BaseRep16/data/qiime/Control"

#output control directory 

OutdirControl <- "/home/anna/metagenome/BaseRep16/out"

# control fam, g, sp, otu files
fam_ctrl_file <- "otu_table_L5.txt"
g_ctrl_file <- "otu_table_L6.txt"
sp_ctrl_file <- "otu_table_L7.txt"
otu_ctrl_file <- "otu_table_MOD.txt"

# control fam, g, sp, otu tables
fam_ctrl_otu_table <- (paste(data_case, fam_ctrl_file, sep = "/")) 
g_ctrl_otu_table <- (paste(data_case, g_ctrl_file, sep = "/")) 
sp_ctrl_otu_table <- (paste(data_case, sp_ctrl_file, sep = "/")) 
otu_ctrl_otu_table <- (paste(data_case, otu_ctrl_file, sep = "/")) 

# control alpha diversity
alpha_div_ctrl <- "data/qiime/Control/alpha_collated/chao1_MOD.txt"


