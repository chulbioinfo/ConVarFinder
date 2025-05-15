# Please set workding direcotry into "6.CorrelationPlots" in zipfile
setwd("./")


# Library
#install.packages("rlang")
#install.packages("car")
library(car)


# Fig4ab
# input data
data <- read.table("1000ctrl/0.summary_Phy_CSAV_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_phylo_info = cbind.data.frame(ProOriBL,ProTerBL,DistanceTN,DistanceTB)
colnames(tmp_phylo_info)<-c("POB","PTB","DTN","DTB")
tmp_data_TAAS = cbind.data.frame(TV_C, TV_D, TV_A)
colnames(tmp_data_TAAS)<-c("iCSAV","dCSAV","CSAV")
detach(data)

data <- read.table("1000ctrl/0.summary_Phy_CSCV_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_data_TCC= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_TCC)<-c("iCSCV","dCSCV","CSCV")
detach(data)

data <- read.table("1000ctrl/0.summary_Phy_CSNV_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_data_TSNV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_TSNV)<-c("iCSNV","dCSNV","CSNV")
detach(data)


tmp_data <- cbind.data.frame(tmp_data_TAAS$iCSAV,
                             tmp_data_TAAS$dCSAV,
                             tmp_data_TAAS$CSAV,
                             tmp_data_TCC$iCSCV,
                             tmp_data_TCC$dCSCV,
                             tmp_data_TCC$CSCV,
                             tmp_data_TSNV$iCSNV,
                             tmp_data_TSNV$dCSNV,
                             tmp_data_TSNV$CSNV,
                             tmp_phylo_info$POB,
                             tmp_phylo_info$PTB,
                             tmp_phylo_info$DTN,
                             tmp_phylo_info$DTB
)
colnames(tmp_data)<-c("iCSAV","dCSAV","CSAV","iCSCV","dCSCV","CSCV","iCSNV","dCSNV","CSNV","POB","PTB","DTN","DTB")
rownames(tmp_data)<-data$nTargets


R_plot <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var, tmp_data){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n", xaxt='n', yaxt='n',ann=FALSE)
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = ""
    if(tmp_p < 0.05){
      tmp_tmp_p = "*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "***"
        }
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend("left",legend=rp, bty='n', text.col = c("black","red"),cex=1.5)
}

correlation_plot <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var, tmp_data){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n", yaxt='n', xaxt='n')
  #grid(lwd = 0.3)
  points(tmp_ind_var,
         tmp_dep_var,
         col="grey",pch=20)

  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line

  points(tmp_ind_var[rownames(tmp_data)=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[rownames(tmp_data)=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)

  outlierTest(dep_ind_var.lm,n.max=3)
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm,n.max=3)$rstudent))))
  points(tmp_ind_var[outlier_list], tmp_dep_var[outlier_list], col="black",pch=4)
  #print(cbind(as.vector(data$nTargets[outlier_list]),as.vector(tmp_ind_var[outlier_list]),as.vector(tmp_dep_var[outlier_list])))
}

correlation_plot_with_XYax <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var, tmp_data){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  #grid(lwd = 0.3)
  points(tmp_ind_var,
         tmp_dep_var,
         col="grey",pch=20)

  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line

  points(tmp_ind_var[rownames(tmp_data)=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[rownames(tmp_data)=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)

  outlierTest(dep_ind_var.lm,n.max=3)
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm,n.max=3)$rstudent))))
  points(tmp_ind_var[outlier_list], tmp_dep_var[outlier_list], col="black",pch=4)
  #print(cbind(as.vector(data$nTargets[outlier_list]),as.vector(tmp_ind_var[outlier_list]),as.vector(tmp_dep_var[outlier_list])))
}

correlation_plot_with_Xax <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var, tmp_data){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n", yaxt='n')
  #grid(lwd = 0.3)
  points(tmp_ind_var,
         tmp_dep_var,
         col="grey",pch=20)

  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line

  points(tmp_ind_var[rownames(tmp_data)=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[rownames(tmp_data)=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)

  outlierTest(dep_ind_var.lm,n.max=3)
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm,n.max=3)$rstudent))))
  points(tmp_ind_var[outlier_list], tmp_dep_var[outlier_list], col="black",pch=4)
  #print(cbind(as.vector(data$nTargets[outlier_list]),as.vector(tmp_ind_var[outlier_list]),as.vector(tmp_dep_var[outlier_list])))
}

correlation_plot_with_Yax <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var, tmp_data){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n",xaxt='n')
  #grid(lwd = 0.3)
  points(tmp_ind_var,
         tmp_dep_var,
         col="grey",pch=20)

  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line

  points(tmp_ind_var[rownames(tmp_data)=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[rownames(tmp_data)=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)

  outlierTest(dep_ind_var.lm,n.max=3)
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm,n.max=3)$rstudent))))
  points(tmp_ind_var[outlier_list], tmp_dep_var[outlier_list], col="black",pch=4)
  #print(cbind(as.vector(data$nTargets[outlier_list]),as.vector(tmp_ind_var[outlier_list]),as.vector(tmp_dep_var[outlier_list])))
}


# correlation test matrix
data_matrix = matrix(c(1:(13*13)),nrow=13,ncol=13)
for(i in 1:length(colnames(tmp_data))){
  data_list = list()
  for(j in 1:length(colnames(tmp_data))){
    if(i==j){
      data_matrix[i,j] = NA
    }else{
      dep_ind_var.cor = cor.test(tmp_data[,i],tmp_data[,j],method="spearman")
      if(i<j){
        tmp_p = dep_ind_var.cor$p.value
        data_matrix[i,j] = tmp_p
      }else{
        tmp_rho = dep_ind_var.cor$estimate
        data_matrix[i,j] = tmp_rho
      }
    }
  }
}
colnames(data_matrix)<-colnames(tmp_data)
rownames(data_matrix)<-colnames(tmp_data)

data_matrix
write.csv(data_matrix, file = "FigS6a.csv")



# regression test matrix
data_matrix = matrix(c(1:(13*13)),nrow=13,ncol=13)
for(i in 1:length(colnames(tmp_data))){
  data_list = list()
  for(j in 1:length(colnames(tmp_data))){
    if(i==j){
      data_matrix[i,j] = NA
    }else{
      dep_ind_var.lm = lm(tmp_data[,j] ~ tmp_data[,i])
      lm.sum = summary(dep_ind_var.lm)
      if(i<j){
        tmp_p = lm.sum$coefficients[2,4]
        data_matrix[i,j] = tmp_p
      }else{
        tmp_r2 = lm.sum$adj.r.squared
        data_matrix[i,j] = tmp_r2
      }
    }
  }
}
colnames(data_matrix)<-colnames(tmp_data)
rownames(data_matrix)<-colnames(tmp_data)

data_matrix
write.csv(data_matrix, file = "FigS6c.csv")





pdf("./Fig4_1000randomctrl_v20201022.pdf",width=16,height=16)
par(mfrow=c(length(colnames(tmp_data)),length(colnames(tmp_data))), mar = c(2.2, 2.2, 0, 0))
attach(tmp_data)
for(i in 1:length(colnames(tmp_data))){
  for(j in 1:length(colnames(tmp_data))){
    if(i==j){
      #par(oma=c(0,0,0,0))
      hist(tmp_data[,i],breaks=24,col="darkgrey", border="darkgrey",main="",xlab = "", ylab = "", axes = TRUE, plot = TRUE)
      legend("topleft",legend = colnames(tmp_data)[i],bty='n',text.col ="blue",cex=1.5)
    }else{
      if(i<j){
        #par(oma=c(0,0,0,0))
        R_plot(tmp_data[,j],tmp_data[,i],"","","",tmp_data)
      }else{
        if(j==1){
          if(i==length(colnames(tmp_data))){
            #par(oma=c(2,2,0,0))
            correlation_plot_with_XYax(tmp_data[,i],tmp_data[,j],"","","",tmp_data)
          }else{
            #par(oma=c(0,2,0,0))
            correlation_plot_with_XYax(tmp_data[,i],tmp_data[,j],"","","",tmp_data)
          }
        }else{
          if(i==length(colnames(tmp_data))){
            #par(oma=c(2,0,0,0))
            correlation_plot_with_XYax(tmp_data[,i],tmp_data[,j],"","","",tmp_data)
          }else{
            #par(oma=c(0,0,0,0))
            correlation_plot_with_XYax(tmp_data[,i],tmp_data[,j],"","","",tmp_data)
          }
        }
      }
    }
  }
}
dev.off()






# Fig4cd
# input data

data <- read.table("61ctrl/0.summary_Phy_CSAV_20201011.txt",sep='\t',header = T)
attach(data)
tmp_phylo_info = cbind.data.frame(ProOriBL,ProTerBL,DistanceTN,DistanceTB)
colnames(tmp_phylo_info)<-c("POB","PTB","DTN","DTB")
tmp_data_TAAS = cbind.data.frame(TV_C, TV_D, TV_A)
colnames(tmp_data_TAAS)<-c("iCSAV","dCSAV","CSAV")
detach(data)

data <- read.table("61ctrl/0.summary_Phy_CSCV_20201011.txt",sep='\t',header = T)
attach(data)
tmp_data_TCC= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_TCC)<-c("iCSCV","dCSCV","CSCV")
detach(data)

data <- read.table("61ctrl/0.summary_Phy_CSNV_20201011.txt",sep='\t',header = T)
attach(data)
tmp_data_TSNV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_TSNV)<-c("iCSNV","dCSNV","CSNV")
detach(data)


tmp_data <- cbind.data.frame(tmp_data_TAAS$iCSAV,
                             tmp_data_TAAS$dCSAV,
                             tmp_data_TAAS$CSAV,
                             tmp_data_TCC$iCSCV,
                             tmp_data_TCC$dCSCV,
                             tmp_data_TCC$CSCV,
                             tmp_data_TSNV$iCSNV,
                             tmp_data_TSNV$dCSNV,
                             tmp_data_TSNV$CSNV,
                             tmp_phylo_info$POB,
                             tmp_phylo_info$PTB,
                             tmp_phylo_info$DTN,
                             tmp_phylo_info$DTB
)
colnames(tmp_data)<-c("iCSAV","dCSAV","CSAV","iCSCV","dCSCV","CSCV","iCSNV","dCSNV","CSNV","POB","PTB","DTN","DTB")
rownames(tmp_data)<-data$nTargets


# correlation test matrix
data_matrix = matrix(c(1:(13*13)),nrow=13,ncol=13)
for(i in 1:length(colnames(tmp_data))){
  data_list = list()
  for(j in 1:length(colnames(tmp_data))){
    if(i==j){
      data_matrix[i,j] = NA
    }else{
      dep_ind_var.cor = cor.test(tmp_data[,i],tmp_data[,j],method="spearman")
      if(i<j){
        tmp_p = dep_ind_var.cor$p.value
        data_matrix[i,j] = tmp_p
      }else{
        tmp_rho = dep_ind_var.cor$estimate
        data_matrix[i,j] = tmp_rho
      }
    }
  }
}
colnames(data_matrix)<-colnames(tmp_data)
rownames(data_matrix)<-colnames(tmp_data)

data_matrix
write.csv(data_matrix, file = "FigS6b.csv")



# regression test matrix
data_matrix = matrix(c(1:(13*13)),nrow=13,ncol=13)
for(i in 1:length(colnames(tmp_data))){
  data_list = list()
  for(j in 1:length(colnames(tmp_data))){
    if(i==j){
      data_matrix[i,j] = NA
    }else{
      dep_ind_var.lm = lm(tmp_data[,j] ~ tmp_data[,i])
      lm.sum = summary(dep_ind_var.lm)
      if(i<j){
        tmp_p = lm.sum$coefficients[2,4]
        data_matrix[i,j] = tmp_p
      }else{
        tmp_r2 = lm.sum$adj.r.squared
        data_matrix[i,j] = tmp_r2
      }
    }
  }
}
colnames(data_matrix)<-colnames(tmp_data)
rownames(data_matrix)<-colnames(tmp_data)

data_matrix
write.csv(data_matrix, file = "FigS6d.csv")


pdf("./FigS5_61corectrl_v20201022.pdf",width=16,height=16)
par(mfrow=c(length(colnames(tmp_data)),length(colnames(tmp_data))), mar = c(2.2, 2.2, 0, 0))
attach(tmp_data)
for(i in 1:length(colnames(tmp_data))){
  for(j in 1:length(colnames(tmp_data))){
    if(i==j){
      #par(oma=c(0,0,0,0))
      hist(tmp_data[,i],breaks=24,col="darkgrey", border="darkgrey", main="",xlab = "", ylab = "", axes = TRUE, plot = TRUE)
      legend("topleft",legend = colnames(tmp_data)[i],bty='n',text.col ="blue",cex=1.5)
    }else{
      if(i<j){
        #par(oma=c(0,0,0,0))
        R_plot(tmp_data[,j],tmp_data[,i],"","","",tmp_data)
      }else{
        if(j==1){
          if(i==length(colnames(tmp_data))){
            #par(oma=c(2,2,0,0))
            correlation_plot_with_XYax(tmp_data[,i],tmp_data[,j],"","","",tmp_data)
          }else{
            #par(oma=c(0,2,0,0))
            correlation_plot_with_XYax(tmp_data[,i],tmp_data[,j],"","","",tmp_data)
          }
        }else{
          if(i==length(colnames(tmp_data))){
            #par(oma=c(2,0,0,0))
            correlation_plot_with_XYax(tmp_data[,i],tmp_data[,j],"","","",tmp_data)
          }else{
            #par(oma=c(0,0,0,0))
            correlation_plot_with_XYax(tmp_data[,i],tmp_data[,j],"","","",tmp_data)
          }
        }
      }
    }
  }
}
dev.off()
