###########################################################
### Loading case and control, alpha ad beta diversity  ###
###########################################################

#otu - operational taxonomic unit 

#--- loading case and control alpha diversity ---      
AlphaDivCase <- LoadAlphaDiv(alpha_div_case)
#--- end loading alpha diversity ---

#--- loading case ---                       
LoadCase <- function()
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

l <- LoadCase(); family_case <- l[[1]]; genus_case <- l[[2]]; species_case <- l[[3]]; otu_case <- l[[4]]; otup_case <- l[[5]]; rm(l)
#--- end loading case ---


#--- loading control --- 
LoadControl <- function()
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

l1 <- LoadControl(); family_ctrl <- l1[[1]]; genus_ctrl <- l1[[2]]; species_ctrl <- l1[[3]]; 
      otu_ctrl <- l1[[4]]; otup_ctrl <- l1[[5]]; rm(l1)
#--- end loading control ---


#--- preparation metadate ---
#case
#meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par))))
#nrow(meta)
#case+control
#meta <- cbind (1:nrow(ff), c(rep("case",length(tags_par)), rep("control",length(tags_ctrl))))
#--- end of preparation metadate ---

###################################################### 
## beta diversity                                   ##
######################################################

############## CHOOSE TAX LEVEL##############################
#sp, g, fam

ff <- family_case[TagsCaseFam] 

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


###### mean abundances
mean.abund.df<-as.data.frame(colMeans(ff_par))
colnames(mean.abund.df)<-"mean.percentage.abundances"
mean.abund.df$standard.deviation<-colSds(ff_par)
mean.abund.df<-mean.abund.df[order(mean.abund.df$mean.percentage.abundances, decreasing=T),]
write.table(mean.abund.df, '/home/anna/metagenome/parkinson_diseas/BR_park/out/mean_abundances.txt',quote=F,sep='\t')

##########################################################################################
# end of testing
##########################################################################################




