##################################
# functions for NeededSampleSize #
##################################

library(HMP)
library(futile.logger)

#adding features from another cohort and alphabetic ordering of features
merge.features <- function (group1,group2,eps){
  dataCase <- group1
  dataCntrl <- group2
  
  featCase <- colnames(dataCase)
  featCntrl <- colnames(dataCntrl)
  newfeatCntrl <- setdiff(featCase,featCntrl)
  newfeatCase <- setdiff(featCntrl,featCase)
  if(is.null(newfeatCase) && is.null(newfeatCntrl)) {
    dataCase <- dataCase[,order(colnames(dataCase))]
    dataCntrl <- dataCntrl[,order(colnames(DataControl))]
  }
  else{
    n<-nrow(dataCase)
    for (i in 1:length(newfeatCase)) {
      dataCase <- cbind(dataCase,rep(eps,n))
    }
    colnames(dataCase) <- c(featCase, newfeatCase)
    dataCase <- dataCase[,order(colnames(dataCase))]
    
    n <- nrow(dataCntrl)
    for (i in 1:length(newfeatCntrl)) {
      dataCntrl <- cbind(dataCntrl,rep(eps,n))
      i<-i+1
    } 
    colnames(dataCntrl) <- c(featCntrl, newfeatCntrl)
    dataCntrl <- dataCntrl[,order(colnames(dataCntrl))]
  }  
  list(dataCase,dataCntrl)
}

# Calculates the statistics, p-value and power of HMP comparisson of 2 groups,
# input: two matrices of counts with equal number of reads per sample, 
#        MC is the number of  Monte-Carlo experiments for HMP. In practice this should be at least 1,000.
# output: list containing staticstics value, p-value and power of MP comparisson of 2 groups
HMPcompare2gr <- function (gr1, gr2, MC = 1000){
  dataCase <- gr1
  dataCntrl <- gr2
  nrsCase <- rowSums(gr1)
  nrsCntrl <- rowSums(gr2)
  
  #comparing 2 groups with HMP
  mygroup <- list(dataCase, dataCntrl)
  compareCaseCntrl2 <- Xmcupo.sevsample(mygroup, ncol(gr1))
  print(compareCaseCntrl2)
  
  #calculating power
  fit.Case <- DM.MoM(dataCase)
  fit.Cntrl <- DM.MoM(dataCntrl)
  
  Nrs <- list(nrsCase, nrsCntrl)
  pi_2gr <- rbind(fit.Case$pi, fit.Cntrl$pi)
  group.theta <- rbind(fit.Case$theta,fit.Cntrl$theta)
  
  power <- MC.Xmcupo.statistics(Nrs, MC, fit.Cntrl$pi, pi_2gr, group.theta, "ha", 0.05)
  
  #summering results
  compareCaseCntrl2[3]  <- power
  names(compareCaseCntrl2)[3] <- "power"
  compareCaseCntrl2
}

#calculates the power of two groups comparison with HMP-package for different number for samples in on of them
#input: matrices of counts with equal number of reads per sample, step for increasing size of Case
#output: matrix with the first column corresponding
power.curve <- function (Case, Cntrl, step = 1, MC = 1500){
  # data <- merge.features(Case,Cntrl, 1e-10)
  group1 <- Case
  group2 <- Cntrl
  n1 <- nrow(group1)
  n2 <- nrow(group2)
  nc <- ncol(group1)
  s <- 3
  power <- c()
  nsamples <- c()
  i <- 1
  #takes n samples form each group and calculates power for them
  while (s <= n1)
  {
    usedrows1 <- c(1:s)
    usedrows2 <- c(1:n2)
    gr1 <- group1[usedrows1,]
    gr2 <- group2[usedrows2,]
    
    #checks wether there are unpresented in subgroups features
    for (l in 1:nc) if (sum(gr1[,l]) == 0) gr1[1,l] <- 1e-15
    for (l in 1:nc) if (sum(gr2[,l]) == 0) gr2[1,l] <- 1e-15
    
    #gr2[,which(colSums(gr2) == 0)] <- 1e-15
    
    
    power[i] <- HMPcompare2gr(gr1,gr2, MC)[[3]]
    nsamples[i] <- s 
    print (c(nsamples[i], power[i]))
    i <- i+1
    s <- s + step
  }
  if(s != n1 + step)
  {
    power[i] <- HMPcompare2gr(group1,group2, MC)[[3]]
    nsamples[i] <- n1
  }
  
  res <- matrix(c(nsamples, power), ncol = 2)
  colnames(res) <- c("number of samples", "power")
  plot(res[,1], res[,2])
  res
}

flog.info("supplementary functions")
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

flog.info("loading QIIME functions")
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

# Loading case and control
# input: pathway with relative representation matrix  
# output: list with matrix
Load <- function(pathway)
{
  flog.info("start load table")
  family <- UniteMatrices(ReadQiimeSumFeats(pathway$Case$FamCaseOtuTbl), 
                          ReadQiimeSumFeats (pathway$Ctrl$FamCtrlOtuTbl))
  genus <- UniteMatrices(ReadQiimeSumFeats (pathway$Case$GenCaseOtuTbl), 
                         ReadQiimeSumFeats (pathway$Ctrl$GenCtrlOtuTbl))
  species <- UniteMatrices(ReadQiimeSumFeats (pathway$Case$SpeCaseOtuTbl), 
                           ReadQiimeSumFeats (pathway$Ctrl$SpeCtrlOtuTbl))
  otu <- UniteMatrices(ReadQiimeOtuTableNoTax (pathway$Case$OtuCaseTbl), 
                       ReadQiimeOtuTableNoTax (pathway$Ctrl$OtuCtrlTbl))
  alphaDiv <- UniteMatrices(LoadAlphaDiv(pathway$Case$AlphaDivCase), 
                            LoadAlphaDiv(pathway$Ctrl$AlphaDivCtrl))
  
  # percent OTU
  otup <- 100 * otu / rowSums(otu)
  # load meta data
  metaTable <- rbind(read.csv(pathway$Case$MetaCaseCsv, stringsAsFactors=F), 
                     read.csv(pathway$Ctrl$MetaCtrlCsv, stringsAsFactors=F))
  
  flog.info("finished load table")
  #TODO(Anna) create checking ...
  list(family=family, genus=genus, species=species, otu=otu, otup=otup, 
       meta=metaTable, alphaDiv=alphaDiv)
}

### Test functions

#calculates needed number of samples in group1 to compare it with group2 with some needed power, when all data is known.
samp_size <- function(group1, group2, needed_power = 0.8, MC = 1500){
  samp_pow <- power.curve(group1, group2, 1, MC)
  n_samp1 <- nrow(samp_pow)
  if (samp_pow[n_samp1,2] < needed_power) stop ("not enough samples to make a decision: power does not exeeds needed value")
  
  #flaggs the beginning of power curve's reliable part 
  flag <- n_samp1
  i <- n_samp1-1
  while (samp_pow[i,2] <= samp_pow[i+1,2] && i > 0){
    flag <- i+1
    i <- i-1
  }
  if (flag >= n_samp1-1) stop("not enough samples to make a desicion: unreliable power curve")
  
  #find the minimal sufficient number of samples
  i <- n_samp1
  while (samp_pow[i]>= needed_power && i>=flag) i <- i-1
  res <- samp_pow[i,1]
}
# test data for adding families 
#dataCase <- cbind(rep(1,5),rep(2,5),rep(3,5))
#dataCntrl <- cbind (rep(4,5),rep(5,5),rep(6,5),rep(7,5))
#colnames(dataCase) <- c("asdtsdg","bafgadfh","csdhs")
#colnames(dataCntrl) <- c("baadg","d","asdtsdg","e")

