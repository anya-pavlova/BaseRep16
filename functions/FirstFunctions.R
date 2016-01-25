################################################################
## UniteMatrices, ReadIni, LoadAlphaDiv, read_qiime functions ##
################################################################

#--- UniteMatrices ---
#-input: two matrix  
#-output: united matrix
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
#--- end UniteMatrices ---

#--- ReadIni ---
#-input: .ini file 
#-output: list
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
#--- end ReadIni ---

#-----WriteTable: write TOP features in .txt file-----
#TRA - table of the relative abundance
WriteTable <- function (TRA, outdir, type)
{
  filename<-paste(outdir,"/",type, '.txt', sep="")
  write.table(TRA, filename, quote=F, sep='\t')
}
#--- end WriteTable ---

#---LoadAlphaDiv: loading alpha diversity---
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
#---end LoadAlphaDiv---

#--- QIIME functions ---

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

#--- chooseTOPfeature ---

#-choose top features with total % of abundance across all samples
#-input: feature vectors (table), percantage
#-output: table with chosen features
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
  y <- dat[,name_list]
  return(y)
}
#--- end chooseTOPfeature ---