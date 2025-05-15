setwd("E:\\GoogleDrive\\Research\\2020\\Avian_vocal_learning_20201020\\Final_Draft\\GitHub\\VocalLearningBirds_2020\\Vocallearningbirds_2020\\12.HumanSAV\\")
library(ggplot2)


data = read.table("summary_humanTAAS_withNonhuman_Type_random_core.txt",header=T,sep='\t')
attach(data)
summary(data)
summary(data$HumanTAAS.HumanOrtho[data$Types=="Control"])

pdf("Fig8_20201107.pdf",width=7,height=2.5)
p10 <- ggplot(data,aes(x = Types, y = HumanTAAS.HumanOrtho)) +
  geom_boxplot(colour = "blue", size=0.5)
p10+ geom_jitter(colour = "black", size=0.5)
summary(p10)
dev.off()
detach(data)

