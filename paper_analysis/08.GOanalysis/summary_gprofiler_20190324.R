setwd("E:\\GoogleDrive\\Research\\2020\\Avian_vocal_learning_20201020\\Final_Draft\\GitHub\\VocalLearningBirds_2020\\Vocallearningbirds_2020\\8.GOanalysis\\")

library(ggplot2)
data = read.table("summary_gprofiler_20190323_gene_and_gobp_over1.txt",header=T,sep='\t')
attach(data)
summary(data)


pdf("Fig6_20190323.pdf",width=5,height=5)


cor_value_cntGene_cntGO <- cor(CNT_gene, CNT_gobp)
lm_cntGene_cntGO <- lm(CNT_gene~CNT_gobp)
lm.sum = summary(lm_cntGene_cntGO)
tmp_r2 = lm.sum$adj.r.squared
tmp_p = lm.sum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                   list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == Tmp_Pvalue),
                   list(Tmp_Pvalue = format(tmp_p, digits = 2)))[2]
ggplot(data = data, aes(x = CNT_gene, y = CNT_gobp))+#, size = Average_p)) + 
  geom_point(alpha = 0.2, col = "red") + geom_smooth(method="lm") + annotate(x=125, y=40, label=rp[1], geom="text", size=5) + annotate(x=125, y=37, label=rp[2], geom="text", size=5) 
  #scale_size_discrete(range = c(1,10))



cor_value_cntGene_meanP <- cor(CNT_gene, Average_p)
lm_cntGene_p <- lm(CNT_gene~Average_p)
lm.sum = summary(lm_cntGene_p)
tmp_r2 = lm.sum$adj.r.squared
tmp_p = lm.sum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                   list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == Tmp_Pvalue),
                   list(Tmp_Pvalue = format(tmp_p, digits = 2)))[2]

ggplot(data = data, aes(x = CNT_gene, y = Average_p)) +#, size = CNT_gobp)) + 
  geom_point(alpha = 0.2, col = "red") + geom_smooth(method="lm") + annotate(x=750, y = 0.04, label=rp[1], geom="text", size=5)+ annotate(x=750, y = 0.035, label=rp[2], geom="text", size=5)
#scale_size_discrete(range = c(1,10))


dev.off()
