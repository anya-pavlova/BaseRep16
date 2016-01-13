###########################################################
### Loading case and control, alpha ad beta diversity  ###
###########################################################

#otu - operational taxonomic unit 

#---loading case and control alpha diversity---      
AlphaDivCase <- LoadAlphaDiv(alpha_div_case)
#---end loading alpha diversity

#---loading case---                       
load_case <- function()
{  
  # NON-RARIFIED
  family_case <- read_qiime_sum_feats (fam_case_otu_table)  
  genus_case <- read_qiime_sum_feats (g_case_otu_table) 
  species_case <- read_qiime_sum_feats (sp_case_otu_table)   
  otu_case <- read_qiime_otu_table_no_tax (otu_case_otu_table) 
    
  # percent OTU
  otup_case <- 100 * otu_case / rowSums(otu_case)
  
  list(family_case, genus_case, species_case,otu_case, otup_case)
}


LoadOtuTable <- function(FamilyInpdir, GenusInpdir, SpeciesInpdir, OtuIntdir)
{
  family <- read_qiime_sum_feats (FamilyInpdir)  
  genus <- read_qiime_sum_feats (GenusInpdir) 
  species <- read_qiime_sum_feats (SpeciesInpdir)   
  otu <- read_qiime_otu_table_no_tax (OtuIntdir) 
  otup <- 100 * otu / rowSums(otu)
  l <- list(family, genus, species,otu, otup)  
  family_case <- l[[1]]; genus_case <- l[[2]]; species_case <- l[[3]]; otu_case <- l[[4]]; otup_case <- l[[5]]; 
}

LoadOtuTable (fam_case_otu_table, g_case_otu_table, sp_case_otu_table, otu_case_otu_table)

l <- load_case(); family_case <- l[[1]]; genus_case <- l[[2]]; species_case <- l[[3]]; otu_case <- l[[4]]; otup_case <- l[[5]]; rm(l)
#---end loading case---


#---loading control--- 
load_control <- function()
{  
  # NON-RARIFIED
  family_ctrl <- read_qiime_sum_feats (fam_ctrl_otu_table)
  genus_ctrl <- read_qiime_sum_feats (g_ctrl_otu_table) 
  species_ctrl <- read_qiime_sum_feats (sp_ctrl_otu_table) 
  otu_ctrl <- read_qiime_otu_table_no_tax (otu_ctrl_otu_table)
  
  # percent OTU
  otup_ctrl <- 100 * otu_ctrl / rowSums(otu_ctrl)
  
  list(species_ctrl, genus_ctrl, family_ctrl, otu_ctrl, otup_ctrl)
}

l1 <- load_control(); family_ctrl <- l1[[1]]; genus_ctrl <- l1[[2]]; species_ctrl <- l1[[3]]; 
      otu_ctrl <- l1[[4]]; otup_ctrl <- l1[[5]]; rm(l1)
#---end loading case and control---

#---preparation metadate just case
meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par))))
nrow(meta)

#---preparation metadate with control
meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par)), rep("control",length(tags_ctrl))))


tags_par <-rownames(sp)[-grep("os$|^sp|^park|^15.1|^15.2|^15.3|^19.1|^19.2", rownames(sp))]
#tags_ctrl <- sort(rownames(sp_ctrl))
#sp<-sp[tags_par,]
length(tags_par)
length(rownames(sp))
length(grep("os$|^sp|^15.1|^15.2|^15.3|^19.1|^19.2", rownames(sp)))


#sp_all<-unite_feat_matrices(sp,sp_ctrl)
#ff <- sp_all[c(tags_par, tags_ctrl),] #sp_all?
#ff <- sp[, which(colSums(sp) > 0)]

#rownames(ff) <- paste(meta_ourn[rownames(ff), "id_timed"],  meta_ourn[rownames(ff), "FIO"], sep=" |") 


load_alpha_rar_case <- function()
{
  alpha_case <- read_qiime_single_alpha_rar(alpha_div_case) 
  alpha_rar_case <- rbind(alpha_case)
  alpha_rar_case
}

f <- read.table(alpha_div_case, skip=2, header=T, row.names=1,sep="\t", colClasses = "character")

f <- read.table(alpha_div_case, header=T, row.names=1,sep="\t", colClasses = "character",stringsAsFactors=F)
f <- t(f)
f.df <- as.data.frame(f)
f.df<-as.data.frame(f.df[-c(1:2),])

rownames(f.df)<-gsub("X","", rownames(f.df))
colnames(f.df)<-"alpha_diversity"
f.df$alpha_diversity<-as.numeric(as.character(f.df$alpha_diversity))

write.table(f_d, '/home/anna/metagenome/BaseRep16/out/alpha_div.txt',quote=F,sep='\t')

mean_alpha <- mean(f.df$alpha_diversity)
sd_alpha <- sd(f.df$alpha_diversity)


cairo_pdf('/home/anna/metagenome/parkinson_diseas/BR_park/graphs/alpha_boxplot.jpg', width = 10,  height = 10)
#pdf('/home/anna/metagenome/parkinson_diseas/BR_park/graphs', pagecentre = FALSE)
boxplot(f.df$alpha_diversity)
dev.off()
#load_alpha_rar_control <- function()
#{
#  alpha_control <- read_qiime_single_alpha_rar(".....................") 
#  alpha_rar_control <- rbind(alpha_control)
#  x <- alpha_rar_control
#}

#alpha_control_tbl <- read.table ("/home/anna/metagenome/parkinson_diseas/BR_park/data/qiime/MIPT_Tomsk_2015_07_30/alpha_collated/chao1_MOD.txt", 
#                                 comment.char = "", quote ="", as.is = T, skip = 2,header = F)


#alpha_diversity_table <- read.table("/home/anna/metagenome/Viliuisk/data/qiime/MIPT_2015_10_08/alpha_collated//t2.txt",
#                                    comment.char = "", quote ="", as.is = T, skip = 2,
#                                    header = F)

#load_alpha_rar_control()
#boxplot(alpha_control)



# read data from QIIME






###################################################### 
## beta diversity                                   ##
######################################################

# названия образцов интересующих нас для исследования и названия образцов контроля
#tags_ctrl <- sort(rownames(sp_ctrl))
#tags_vil <- sort(rownames(sp)[grep("^V", rownames(sp))])
#tags_Y89 <- sort(rownames(sp)[grep("^Y8|^Y9", rownames(sp))]) #for Y8 and Y9, yana and Borya 

#tags_par <-rownames(sp)[-grep("^xx.|W$", rownames(sp))] #нужны все образы кроме: xx-*, *-W

tags_par <-rownames(sp)[-grep("os$|^sp|^park|^15.1|^15.2|^15.3|^19.1|^19.2", rownames(sp))]
#tags_ctrl <- sort(rownames(sp_ctrl))
#sp<-sp[tags_par,]
length(tags_par)
length(rownames(sp))

length(grep("os$|^sp|^15.1|^15.2|^15.3|^19.1|^19.2", rownames(sp)))


############## CHOOSE TAX LEVEL##############################
#sp, g, fam

ff <- fam[tags_par,]
ff_par <- fam[tags_par, ]
dd_par <- bcdist(ff_par)

ff_par <- ff_par[, which(colMaxs(ff_par) > 0)]
#ff <- sp[, which(colSums(sp) > 0)]

#rownames(ff) <- paste(meta_ourn[rownames(ff), "id_timed"],  meta_ourn[rownames(ff), "FIO"], sep=" |") 

# BC dist

dd <- bcdist(ff)
dd
#dd.df <- as.data.frame(dd)


# UniFrac
#w_unifrac.m <- read.table("data/weighted_unifrac_merged_otu_tables.txt", header=T, row.names=1, sep="\t", colClasses = "character", as.is = T)
#colnames(w_unifrac.m) <- str_replace(colnames(w_unifrac.m), pattern = "^X", replacement = "")
#w_unifrac.m <- w_unifrac.m[c(tags_ea1, tags_ea2, tags_ab1, tags_ab2), c(tags_ea1, tags_ea2, tags_ab1, tags_ab2)]
#dd <- as.dist(w_unifrac.m)

mdd <- data.matrix(dd)
hr <- hclust(dd, method="ward.D") #The "ward" method has been renamed to "ward.D"; note new "ward.D2"
groups <- cutree(hr, k=2)


########################################################################################
# TEST CODE 
########################################################################################


makeMDS<-function(dis, outdir,meta, colFact,plot_name,m_width = 5, m_height=5)
{   
    myMDS <- isoMDS(dis, k=2)
    outdirMDS<-paste(outdir,'/MDS_',plot_name,'.pdf',sep='')
    #outdirMDS<-paste(outdir,'/MDS_',plot_name,'.jpg',sep='')
    outdirMDS<-file.path(outdirMDS)
    dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]),cols = meta[,colFact], lab = rownames(myMDS$points))
    plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, color=factor(cols), label=lab)) + geom_point(size=0.3) + geom_text(size=2) + theme_bw()		
    ggsave(outdirMDS, plotMDS, width = m_width, height=m_height)  
}

#meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par)), rep("control",length(tags_ctrl))))
meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par))))
nrow(meta)


makeMDS(dd, "/home/anna/metagenome/parkinson_diseas/BR_park/graphs", meta, 2, "case")
meta[,2]
nrow(meta)
dis<-dd
outdir="/home/anna/metagenome/parkinson_diseas/BR_park/graphs"
#length(as.vector(myMDS$points[,1]))

#nrow(sp_ctrl)
#nrow(sp_all)

## just case



meta <- cbind (1:nrow(ff_par), c(rep("case",length(tags_par))))
makeMDS(dd_par,"/home/anna/metagenome/parkinson_diseas/BR_park/graphs", meta, 2, "case")

#############make heatmap
#makeHeat<-function(featVect, meta, colRow,outdir, graph_name)
makeHeat<-function(featVect, meta, outdir, graph_name)
{
  #cols<-rainbow (length(levels(as.factor(meta[,colRow]))))
  #cols<-cbind(levels(as.factor(meta[,colRow])), cols)
  #rownames(cols)<-cols[,1]
  #colVec<-c()
  #for (i in rownames(featVect))
  #{
  #  colVec<-append(colVec, cols[meta[i,colRow],2])
  #}
  outdirHeat<-paste(outdir,'/Heatmap_',colRow,'_',graph_name,'_','.pdf',sep='')
  outdirHeat<-file.path(outdirHeat)
  pdf(outdirHeat, pagecentre = FALSE)
  heatmap.2(as.matrix(featVect), margin=c(20,12), distfun=distSpear, cexRow=.1,col=cols.gentleman(500), cexCol=.6, 
            hclustfun=hclustW, trace='none') #,RowSideColors=colVec)
  dev.off()
}
#########end

rownames(meta)<-rownames(ff)
colSums(ff)
###choose top features with total % of abundance across all samples
####input: feature vectors (table), percantage
####output: table with chosen features

chooseTOPfeature<-function(dat,perc)
{
  summof<-sum(colSums(dat))
  sorted<-sort(colSums(dat), decreasing = TRUE)
  count <-0
  num<-0
  name_list<-c()  
  for (i in sorted) #in sorted by total abundance features, take those making %
  {
    count<-count+((i/summof)*100)
    num <-num+1
    if (count>perc) 
    {
      break
    }
    else 
    {
      name_list<-append(name_list, names(sorted[num]))
    }  
  }	
  y<-dat[,name_list]
}
####end

meta <- cbind (1:nrow(ff_par), c(rep("case",length(tags_par))))
rownames(meta)<-rownames(ff_par)



ff_top <- chooseTOPfeature(ff,85)
ff_par_top <- chooseTOPfeature(ff_par,85)
library(stringr)
#colnames(ff_top)<-str_extract(colnames(ff_top), "(?<=o__).+$")

#colnames(ff_par_top)<-str_extract(colnames(ff_par_top), "(?<=o__).+$")

colnames(ff_top) <- paste("o_", unlist(data.frame(strsplit(colnames(ff_top), "o_"))[2,]), sep="")

ff_par_top<-ff_par_top[,order(colSums(ff_par_top), decreasing = T)]

## HeatMap (family)

#colnames(ff_par_top) <- paste("f_", unlist(data.frame(strsplit(colnames(ff_par_top), "f_"))[2,]), sep="")

drop_unecessary_tax_names <- function(data, level="genera") {
  if (level == "genera") {
    pattern <- "g__.+"
  } else if (level == "phylum") {
    pattern <- "p__.+"
  }
  new.names <- str_extract(colnames(data), pattern)
  good.ids <- which(!is.na(new.names))
  data <- data[,good.ids]
  colnames(data) <- new.names[good.ids]
  return(data)
}

makeHeat(ff_par_top, meta, outdir, "g")

# ff_par.sub <- drop_unecessary_tax_names(ff_par, level = "genera")
# ff_par_top.sub <- chooseTOPfeature(ff_par.sub, perc = 85)
# 
# 
# new.tax.names <- str_replace(colnames(ff_par), pattern = "\\;\\s\\w__\\;.+", replacement = "")
# new.tax.names <- str_replace(new.tax.names, pattern = "\\;\\s\\w__$", replacement = "")
# colnames(ff_par)<-new.tax.names
# 
# tail(colnames(ff_par))
# tail(new.tax.names)
# 
# makeHeat(ff_par_top.sub, meta, 2, outdir, "genus")

## HeatMap (genus)

outdir="/home/anna/metagenome/parkinson_diseas/BR_park/graphs2"

meta <- cbind (1:nrow(ff_par), c(rep("case",length(tags_par))))
rownames(meta)<-rownames(ff_par)
ff_top <- chooseTOPfeature(ff,85)
ff_par_top <- chooseTOPfeature(ff_par,85)
library(stringr)
#colnames(ff_top)<-str_extract(colnames(ff_top), "(?<=o__).+$")
colnames(ff_par_top)<-str_extract(colnames(ff_par_top), "(?<=o__).+$")
ff_par_top<-ff_par_top[,order(colSums(ff_par_top), decreasing = T)]


#ff_par_top.df <- as.data.frame(ff_par_top)
#ff_par_top_sort_df_R <- sort(ff_par_top.df$f__Ruminococcaceae; g__; s__)

outdir
colnames(ff_par_top) <- paste("g_", unlist(data.frame(strsplit(colnames(ff_par_top), "g_"))[2,]), sep="")
makeHeat(ff_par_top, meta, 2, outdir)


## HeatMap (species)

outdir="/home/anna/metagenome/parkinson_diseas/BR_park/graphs3"

meta <- cbind (1:nrow(ff_par), c(rep("case",length(tags_par))))
rownames(meta)<-rownames(ff_par)
ff_top <- chooseTOPfeature(ff,85)
ff_par_top <- chooseTOPfeature(ff_par,85)

#colnames(ff_top)<-str_extract(colnames(ff_top), "(?<=o__).+$")
colnames(ff_par_top)<-str_extract(colnames(ff_par_top), "(?<=o__).+$")
ff_par_top<-ff_par_top[,order(colSums(ff_par_top), decreasing = T)]

outdir
colnames(ff_par_top) <- paste("s_", unlist(data.frame(strsplit(colnames(ff_par_top), "s_"))[2,]), sep="")
makeHeat(ff_par_top, meta, 2, outdir)

############ Статистика сквенирования 
stat_table <- read.table("/home/anna/metagenome/parkinson_diseas/BR_park/stat_table.txt",
                                    comment.char = "", quote ="", as.is = T, header = T)
#среднее значение + отклонение (initial, after_filt, after_rare, mapped)

#initial
mean_init <- mean(stat_table$initial)
sd_init <- sd(stat_table$initial)

#after_filt
mean_af <- mean(stat_table$after_filt)
sd_af <- sd(stat_table$after_filt)

#after_rar
mean_ar <- mean(stat_table$after_rare)
sd_ar <- sd(stat_table$after_rare)

#mapped
mean_map <- mean(stat_table$mapped)
sd_map <- sd(stat_table$mapped)

stat_df = data.frame (mean_init, sd_init, mean_ar, sd_ar, mean_map, sd_map)
write.table(stat_df, '/home/anna/metagenome/parkinson_diseas/BR_park/out/stat_df.txt',quote=F,sep='\t')


boxplot(stat_table$initial)

###### mean abundances
mean.abund.df<-as.data.frame(colMeans(ff_par))
colnames(mean.abund.df)<-"mean.percentage.abundances"
mean.abund.df$standard.deviation<-colSds(ff_par)
mean.abund.df<-mean.abund.df[order(mean.abund.df$mean.percentage.abundances, decreasing=T),]
write.table(mean.abund.df, '/home/anna/metagenome/parkinson_diseas/BR_park/out/mean_abundances.txt',quote=F,sep='\t')

##########################################################################################
# end of testing
##########################################################################################

### making normal txt-files
#alpha_diversity_table <- read.table("/home/anna/metagenome/Viliuisk/data/qiime/MIPT_2015_10_08/alpha_collated//t2.txt",
#                                    comment.char = "", quote ="", as.is = T, skip = 2,
#                                    header = F)
#head(alpha_diversity_table)
#
#con <- file("/home/anna/metagenome/Viliuisk/data/qiime/MIPT_2015_10_08/alpha_collated//t2.txt", "r")
#x <- sapply(readLines(con), strsplit, split = "\t")
#close(con)
#
#w.d <- data.frame(t(as.numeric(x[[2]][-1])))
#colnames(w.d) <- x[[1]][-(1)]
#head(w.d)



