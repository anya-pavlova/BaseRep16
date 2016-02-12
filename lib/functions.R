################################################
# functions                                    #
################################################

#loading libraries
library(futile.logger)
library(lint)
library(stringr)
library(ecodist)
library(data.table)
library(ape)
library(ggplot2)
library(gplots)
library(scales)
library(MASS)
library(matrixStats)
library(gridExtra)
library(grid)
library(phyloseq) #source("https://bioconductor.org/biocLite.R"), biocLite("phyloseq")
library(reshape)

message("loading UniteMatrices, ReadIni, WriteTable functions")
##################
# Functions: 
# UniteMatrix, ReadIni, WriteTable, QIIME functions,     

# UniteMatrices 
# input: two matrix  
# output: united matrix
UniteMatrices <- function(t1, t2)
{
  uc <- sort(union(colnames(t1), colnames(t2)))
  ur <- c(rownames(t1), rownames(t2))
  t <- matrix(0, nrow = nrow(t1) + nrow(t2), ncol = length(uc))
  colnames(t) <- uc
  rownames(t) <- ur
  t[rownames(t1), colnames(t1)] <- t1
  t[rownames(t2), colnames(t2)] <- t2
  identical(t1, t[rownames(t1), colnames(t1)])
  identical(t2, t[rownames(t2), colnames(t2)])
  t
}

# ReadIni
# input: .ini file 
# output: list
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

# WriteTable
# input: some structure
# output: txt file with data frome input file
WriteTable <- function (someStructure, outdir, type)
{
  filename<-paste(outdir, "/", type, '.txt', sep="")
  write.table(someStructure, filename, quote=T, sep='\t')
}

message("loading QIIME functions")
# QIIME functions 
# input: tables frome qiime
# output: matrix with data from input file
ReadQiimeSumFeats <- function(filename)
{
  f <- t(read.table(filename, header=T, row.names=1, sep="\t"))
  rownames(f) <- gsub("X", "", rownames(f))
  f <- data.matrix(f)
  f <- 100 * f / rowSums(f)
  f
}

ReadQiimeOtuTableNoTax <- function(filename)
{   
  f <- read.table(filename, header=T, row.names=1, sep="\t")
  # drop taxonomy
  f <- f[, -ncol(f)]
  f <- t(f)
  rownames(f) <- gsub("X", "", rownames(f))
  f
}

ReadQiimeSingleAlphaRar <- function(filename)
{   
  f <- read.table(filename, header=T, row.names=1, sep="\t", 
                  colClasses = "character")
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

message("loading alpha diversity, Load, chooseTOPfeature, TopFeatures functions")
# loading alpha diversity
LoadAlphaDiv <- function(inpdir)
{
  inpdir
  AlphaDivTbl <- (read.table (inpdir, row.names=1,  header=T,sep="\t",
                              stringsAsFactors=F ))
  AlphaDivTbl <- t(AlphaDivTbl)
  AlphaDivTbl <- as.data.frame(AlphaDivTbl[-c(1, 2),])
  rownames(AlphaDivTbl)<-gsub("X","", rownames(AlphaDivTbl))
  colnames(AlphaDivTbl) <- c("AlphaDiversity")
  return(data.matrix(AlphaDivTbl))
}

# Loading group.one and control
# input: pathway with relative representation matrix  
# output: list with matrix

Load <- function(pathway)
{
  family <- UniteMatrices(ReadQiimeSumFeats(pathway$group.one$fam.group.one.otu.tbl), 
                          ReadQiimeSumFeats (pathway$group.two$fam.group.tow.otu.tbl))
  genus <- UniteMatrices(ReadQiimeSumFeats (pathway$group.one$gen.group.one.otu.tbl), 
                         ReadQiimeSumFeats (pathway$group.two$gen.group.tow.otu.tbl))
  species <- UniteMatrices(ReadQiimeSumFeats (pathway$group.one$spe.group.one.otu.tbl), 
                           ReadQiimeSumFeats (pathway$group.two$spe.group.tow.otu.tbl))
  otu <- UniteMatrices(ReadQiimeOtuTableNoTax (pathway$group.one$otu.group.one.tbl), 
                       ReadQiimeOtuTableNoTax (pathway$group.two$otu.group.tow.tbl))
  alpha.div <- UniteMatrices(LoadAlphaDiv(pathway$group.one$alpha.div.group.one), 
                             LoadAlphaDiv(pathway$group.two$alpha.div.group.tow))
  meta.table <- rbind(read.csv(pathway$group.one$meta.group.one.csv, stringsAsFactors=F),
                      read.csv(pathway$group.two$meta.group.tow.csv, stringsAsFactors=F))
  # percent OTU
  otup <- 100 * otu / rowSums(otu)
  #ToDo(Anna) create checking
  list(family=family, genus=genus, species=species, otu=otu, otup=otup, 
       meta=meta.table, alpha.div=alpha.div)
}

##### check rowSums=100%
#rowSums(total.table$family)

##### check rownames matching for all tables and metadata
#setdiff(rownames(Familygroup.two), Metagroup.two$samples_name)
#setdiff(rownames(Genusgroup.two), Metagroup.two$samples_name)
#setdiff(rownames(Speciesgroup.two), Metagroup.two$samples_name)

#rownames(total.table$family)


# chooseTOPfeature 
# choose top features with total % of abundance across all samples
# input: feature vectors (table), percantage
# output: table with chosen features
chooseTOPfeature<-function(data,perc)
{
  summof<-sum(colSums(data))
  sorted<-sort(colSums(data), decreasing = TRUE)
  count <-0
  num<-0
  name_list<-c()  
  for (i in sorted) #in sorted by total abundance features, take those making %
  {
    count <- count+((i/summof)*100)
    num <- num+1
    if (count > perc) 
    {
      break
    }
    else 
    {
      name_list <- append(name_list, names(sorted[num]))
    }    
  }  
  y <- data[,name_list]
  return(y)
}


# TopFeatures 
# choose top features for group.one and control separately
# input: feature vectors (table), percantage
# output: table with chosen features
TopFeatures <- function(total.table, type, subset, perc)
{
  
  family.top <- chooseTOPfeature(total.table$family[which(rownames(total.table$family)
                                                        %in% total.table$meta[which(total.table$meta[,type] 
                                                                                   %in% subset),"samples_name"]),], perc=perc)
  genus.top <- chooseTOPfeature(total.table$genus[which(rownames(total.table$genus)
                                                      %in% total.table$meta[which(total.table$meta[,type]
                                                                                 %in% subset),"samples_name"]),], perc=perc)
  species.top <- chooseTOPfeature(total.table$species[which(rownames(total.table$species) 
                                                          %in% total.table$meta[which(total.table$meta[,type] 
                                                                                     %in% subset),"samples_name"]),], perc=perc)
  otu.top <- chooseTOPfeature(total.table$otu[which(rownames(total.table$otu)
                                                  %in% total.table$meta[which(total.table$meta[,type] 
                                                                             %in% subset),"samples_name"]),], perc=perc)
    
  list(family.top = family.top, genus.top = genus.top, species.top = species.top, otu.top = otu.top)

 }

#---stastistics ----------------------------------------------------------------
message("loading statistics functions")

# make  distance in Spearman correlation  - for internal use
# input: feature vectors (table)
# output: class dist object
distSpear<-function (x) 
{
  result <- 1-cor(t(x), method='spearman', use='pair')
  as.dist(result)
}

# good colours for heatplot
cols.gentleman <- function(ncol=500) {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, 'RdBu'))(ncol)
  return(rev(hmcol))
}

# make  distance in Spearman correlation  - for internal use
# input: feature vectors (table)
# output: class dist object
distSpear<-function (x) 
{
  result <- 1-cor(t(x), method='spearman', use='pair')
  as.dist(result)
}

# make JS distance - for internal use
# input: feature vectors (table)
# output: class dist object
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


# make all distances
# input: feature vectors (table), vector with names of metrics ('JS','Spear','Eu','Man','Can','BC') chosen
# output: a list with class dist objects, names of the objects - like metrics chosen
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


