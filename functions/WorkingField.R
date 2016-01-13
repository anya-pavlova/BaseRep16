########################################################
## working with alpha, beta -diversity                ##
########################################################

AlphaDivCase <- LoadAlphaDiv(alpha_div_case)
write.table(AlphaDivCase, '/home/anna/metagenome/BaseRep16/out/Case/AlphaDiveTbl.txt',quote=F,sep='\t')

MeanAlphaDiv <- mean(AlphaDivCase$AlphaDiversity)
write.table(AlphaDivCase, '/home/anna/metagenome/BaseRep16/out/Case/AlphaDivMean.txt',quote=F,sep='\t')

SdAlphaDiv <- sd(AlphaDivCase$AlphaDiversity)
write.table(AlphaDivCase, '/home/anna/metagenome/BaseRep16/out/Case/AlphaDivSd.txt',quote=F,sep='\t')

cairo_pdf('/home/anna/metagenome/BaseRep16/out/Case/Graphs/alpha_boxplot.pdf', width = 10,  height = 10)
boxplot(AlphaDivCase$AlphaDiversity)
dev.off()

#-----end alpha diversity-----