rm(list=ls(all=TRUE)); gc()
################################################
# 16S rRNA metagenomes                         #
################################################

current_directory <- "/home/anna/metagenome/BaseRep16"
#current_directory <- ""
setwd(current_directory)

source("lib/functions.R")

#group.one, g.one  (case)
#group.two, g.tow (control)

message("loading pathway")
path.to.config <- "src/pathway.ini"
pathway <- ReadIni(path.to.config) 

path.to.constants <- "src/constants.ini"
constants <- ReadIni(path.to.constants)
as.numeric(constants$constants$perc)
#constants <- lapply(constants, as.numeric)
flog.info("end of loading pathway")

message("working with data")
total.table <- Load(pathway)

table.top <- list(group.one = TopFeatures(total.table, type = "Type.1", 
                                          subset = "case", 
                                          perc = as.numeric(constants$constants$perc)),
                  group.two = TopFeatures(total.table, type = "Type.1", 
                                     subset = "control", perc = as.numeric(constants$constants$perc)))

table.top.total <- list(family.top = UniteMatrices(table.top$group.one$family.top, 
                                                   table.top$group.two$family.top), 
                        genus.top = UniteMatrices(table.top$group.one$family.top, 
                                                  table.top$group.two$family.top), 
                        species.top = UniteMatrices(table.top$group.one$family.top, 
                                                    table.top$group.two$family.top), 
                        otu.top = UniteMatrices(table.top$group.one$family.top, 
                                                table.top$group.two$family.top))
#isn't samples name (?)
WriteTable(table.top.total$family.top, pathway$group.one$outdir,"family.table.top") 
WriteTable(Genusgroup.one, Outdirgroup.one, "Genusgroup.one") 
WriteTable(Speciesgroup.one, Outdirgroup.one, "Speciesgroup.one") 
WriteTable(Otugroup.one, Outdirgroup.one, "Otugroup.one")

message("start working with alpha diversity")
alphaDivgroup.one <- (total.table$AlphaDiv[
  which(rownames(total.table$AlphaDiv) 
        %in% total.table$Meta[which(total.table$Meta[,"Type.1"] 
                                    %in% "group.one"),"samples_name"]),])
]

WriteTable(alphaDivgroup.one, pathway$group.one$Outdirgroup.one, "AlphaDivgroup.oneTbl"))
AlphaDivMean <- lapply(AlphaDivgroup.oneL[match('AlphaDiversity', names(AlphaDivgroup.oneL))], mean)
AlphaDivSd <- lapply(AlphaDivgroup.oneL[match('AlphaDiversity', names(AlphaDivgroup.oneL))], sd)
AlphaDivMeanSd <- c(AlphaDivMean, AlphaDivSd)
names(AlphaDivMeanSd) <- c("mean", "sd")

WriteTable (AlphaDivMeanSd, Outdirgroup.one, "AlphaDivMeanAndSd") 

cairo_pdf((paste(Outdirgroup.one, '/Graphs/alpha_boxplot.pdf', sep = "/")), width = 10,  height = 10)
boxplot(alphaDivgroup.one$AlphaDiversity)
dev.off()

#group.two
alphaDivgroup.one <- (total.table$AlphaDiv[
  which(rownames(total.table$AlphaDiv) 
        %in% total.table$Meta[which(total.table$meta[,"Type.1"] 
                                    %in% "group.two"),"samples_name"]),])
]

message("end start working with alpha diversity")

message("start sequencing statistics")
seqStatTable <- read.table(pathway$statistics$SeqStatTblFile, 
                           comment.char = "", quote ="", as.is = T, header = T)

seqStatList <- as.list(seqStatTable)
seqStatMean <- lapply(seqStatList[-match('samp_name', names(seqStatList))], mean)
seqStatSd <- lapply(seqStatList[-match('samp_name', names(seqStatList))], sd)
names(seqStatMean) <- c("mean_init", "mean_AF", "mean_AR", "mean_Mapp")
names(seqStatSd) <- c("sd_init", "sd_AF", "sd_AR", "sd_Mapp")
seqStatMeanSd <- c(seqStatMean, seqStatSd)
seqStatMeanSd.df <- t(as.data.frame(seqStatMeanSd))

WriteTable(seqStatMeanSd.df, pathway$group.one$Outdirgroup.one, "SeqStatMeanAndSd")
message("end sequencing statistics")

#--- make MDS ---
#будем ли отрисовывать по топовым представленнос
#Family
MDS<-function(dis, outdir, meta, colFact, plot_name,m_width = 5, m_height=5)
{   
  myMDS <- isoMDS(dis, k=2)
  outdirMDS<-paste(outdir,'/MDS_',plot_name,'.pdf',sep='')
  outdirMDS<-file.path(outdirMDS)
  dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]),cols = meta[,colFact], lab = rownames(myMDS$points))
  plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, color=factor(cols), label=lab)) + geom_point(size=0.3) + geom_text(size=2, x=X + 20, y=Y) + theme_bw()    
  ggsave(outdirMDS, plotMDS, width = m_width, height=m_height)  
}

ddHSN<-bcdist(total.table$family)
MDS(ddHSN, pathway$group.one$GraphsDir, 
    total.table$meta[which(rownames(total.table$meta)%in%total.table$meta[whic%in% "case"))] )
MDS(ddHSN, pathway$group.one$GraphsDir, 
    total.table$meta[which(rownames(total.table$meta)%in%total.table$meta[which]%in% "case"),"samples_name"]),], "Type.1", "Family")




ddHSN<-bcdist(total.table$family)
MDS(ddHSN, pathway$group.one$GraphsDir, 
    total.table$meta[which(rownames(total.table$meta) 
                         %in% total.table$meta[which]
                         %in% "group.one"),"samples_name"]),], "Type.1", "Family")
MDS(ddHSN, )

total.table$meta

MDS(ddHSN, pathway$group.one$GraphsDir, 
    total.table$meta[which(rownames(total.table$meta) 
                           %in% total.table$meta[which(total.table$meta[,"Type.1"] 
                                                       %in% "case"),"samples_name"]),],"Type.1", "Family")

total.table$meta[which(total.table$meta[,"Type.1"] == "case"),]$samples_name

message("make MDS")

# make MDS
# input: distance matrix, output directory, meta matrix,  
# output: MDS plot



#Genus

ddHSN<-bcdist(total.table$genus)
MDS(ddHSN, pathway$group.one$GraphsDir, pathway$group.one$meta.group.one.csv, "Type.1", "Genus")

#Species
ddHSN<-bcdist(UniteMatrices(Speciesgroup.one, SpeciesCtrl))
MDS(ddHSN, pathway$group.one$GraphsDir , total.table$MetaTable, "Type.1", "Species")
#--- end ---

#---stastistics---
###make  distance in Spearman correlation  - for internal use
####input: feature vectors (table)
####output: class dist object
distSpear<-function (x) 
{
  result <- 1-cor(t(x), method='spearman', use='pair')
  as.dist(result)
}
###make JS distance - for internal use
####input: feature vectors (table)
####output: class dist object
distJS<-function(x)
{
  result<-c()
  input<-c()
  covs<-x
  covs<-covs+0.000001 # add pseudocount
  covs_norm<-covs/rowSums(covs) #normalize input by 1
  input = t(covs_norm)
  heads<-list(colnames(input))
  result<-matrix(0,nrow = ncol(input), ncol = ncol(input),dimnames = c(heads, heads)) #make result distance table
  kld1<-0
  kld2<-0
  for (i in 1:ncol(input))
  {
    for (p in 1:ncol(input)) 
    {
      xj<-input[,i]
      yj<-input[,p]
      mj<-(xj+yj)/2
      kld1<-kld1 + xj*(log(xj/mj)) #compute Kullback-Leibler coefficients
      kld2<-kld2 + yj*(log(yj/mj))
      result[i,p]<-sqrt(sum((kld1/2) ,(kld2/2))) #compute Jenson-Shannon distance
      kld1<-0
      kld2<-0
    }  
  }
  js_d<-as.dist(result)
  return (js_d)
}
#######end


###make all distances
####input: feature vectors (table), vector with names of metrics ('JS','Spear','Eu','Man','Can','BC') chosen
####output: a list with class dist objects, names of the objects - like metrics chosen
allDist<-function(x, metric=c('JS','Spear','Eu','Man','Can','BC'))
{
  dists<-list()
  num<-0
  for (i in (metric)) 
  {
{
  num<-num+1
  if (i == 'JS')
  {
    dis<-distJS(x)
  }
  if (i == 'Spear')
  {
    dis<-distSpear(x)
  }
  if (i == 'Eu')
  {
    dis<-dist(x)
  }
  if (i == 'Man')
  {
    dis<-dist(x, method='manhattan')
  }
  if (i == 'Can')
  {
    dis<-dist(x, method='canberra')
  }
  if (i == 'BC')
  {
    dis<-bcdist(x)
  }
  dists[[num]]<-dis
}	
  }	
names(dists)<-metric
return (dists)
}