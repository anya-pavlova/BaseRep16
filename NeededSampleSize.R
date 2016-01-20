library(HMP)

#loading data
dataCase <-FamilyCase#/100
dataCase <- dataCase[order(as.numeric(rownames(dataCase))),]
dataCntrl <- FamilyCntrl#/100

#

#adding families from another cohort
fam_case <- colnames(FamilyCase)
fam_cntrl <- colnames(FamilyCntrl)
newfam_cntrl <- setdiff(fam_case,fam_cntrl)
newfam_case <- setdiff(fam_cntrl,fam_case)

n<-nrow(dataCase)
for (i in 1:length(newfam_case)) {
  dataCase <- cbind(dataCase,rep(1e-10,n))
} 
colnames(dataCase) <- c(fam_case, newfam_case)
alph_order_case <- order(colnames(dataCase))
dataCase <- dataCase[,alph_order_case]

n <- nrow(dataCntrl)
for (i in 1:length(newfam_cntrl)) {
  dataCntrl <- cbind(dataCntrl,rep(1e-10,n))
  i<-i+1
} 
colnames(dataCntrl) <- c(fam_cntrl, newfam_cntrl)
alph_order_cntrl <- order(colnames(dataCntrl))
dataCntrl <- dataCntrl[,alph_order_cntrl]

################################################################################
## Test data for adding families ###############################################
################################################################################
#FamilyCase <- cbind(rep(1,5),rep(2,5),rep(3,5))
#FamilyCntrl <- cbind (rep(4,5),rep(5,5),rep(6,5),rep(7,5))
#colnames(FamilyCase) <- c("asdtsdg","bafgadfh","csdhs")
#colnames(FamilyCntrl) <- c("baadg","d","asdtsdg","e")
################################################################################


###Compare Case and Controll  with HMP package####
mygroup <- list(dataCase, dataCntrl)

readsCntrl <- read_qiime_otu_table_no_tax("/home/vera/Documents/metagenom/SampleSize/data/otu_table_control.txt")
statCntrl <- c()
for(i in 1:nrow(readsCntrl)) statCntrl[i] <- sum(readsCntrl[i,])
for (i in 1:ncol(statCase)) colnames(statCase)[i] <- as.character(statCase[1,i])
names(statCntrl) <- rownames(readsCntrl)
Nrs.Cntrl <- statCntrl
names(Nrs.Cntrl) <- c()

statCase <- read.table("/home/vera/Documents/metagenom/SampleSize/data/StatTable_case.txt", header = TRUE)
statCase <- statCase[,-1]
statCase <- statCase[order(as.numeric(statCase[,1])),]
Nrs.Case <- as.numeric(statCase[,4])#здесь нет проверки, что образцы в нужном порядке, но на тестовых данных это так

fit.Case <- DM.MoM(dataCase)
fit.Cntrl <- DM.MoM(dataCntrl)

dataAll <- rbind(dataCase,dataCntrl)
fit.All <- DM.MoM(dataAll)

Nrs <- list(Nrs.Case, Nrs.Cntrl)
group.alphap <- rbind(fit.Case$gamma, fit.Cntrl$gamma)
group.theta <- rbind(fit.Case$theta,fit.Cntrl$theta)
pi_2gr <- rbind(fit.Case$pi, fit.Cntrl$pi)

#сравнение DM-моделей для двух групп
#compareCaseCntrl <- Xdc.sevsample(mygroup, ncol(dataCase))
#compareCaseCntrl
#MC.Xdc.statistics(Nrs, 10000, group.alphap, 2, type = "ha", siglev = 0.05, est = "mom")

#обобщённый тест Вальда для сравнения групповых векторов представленности
compareCaseCntrl2 <- Xmcupo.sevsample(mygroup, ncol(dataCase))
compareCaseCntrl2
MC.Xmcupo.statistics(Nrs, 10000, fit.All$pi, pi_2gr, group.theta, "ha", 0.05)

