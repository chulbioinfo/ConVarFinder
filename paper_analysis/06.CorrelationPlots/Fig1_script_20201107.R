# Please set workding direcotry into "6.CorrelationPlots" in zipfile
setwd("./")

#install.packages("rlang")
#install.packages("car")
library("car")


# Random control sets #########################################################

# Loading data
data <- read.table("1000ctrl/0.summary_Phy_CSAV_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_phylo_info = cbind.data.frame(ProTerBL,DistanceTN,DistanceTB,ProOriBL)
colnames(tmp_phylo_info)<-c("PTB","DTN","DTB","POB")
tmp_data_CSAV = cbind.data.frame(nTargets,ctrl,TV_A, TV_C, TV_D)
colnames(tmp_data_CSAV)<-c("nTargets","ctrl","CSAV","iCSAV","dCSAV")
detach(data)

data <- read.table("1000ctrl/0.summary_Phy_CSCV_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_data_CSCV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_CSCV)<-c("iCSCV","dCSCV","CSCV")
detach(data)

data <- read.table("1000ctrl/0.summary_Phy_CSNV_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_data_CSNV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_CSNV)<-c("iCSNV","dCSNV","CSNV")
detach(data)

data <- cbind.data.frame(tmp_data_CSAV$nTargets,
                                       tmp_data_CSAV$ctrl,
                                       tmp_phylo_info,
                                       tmp_data_CSAV$iCSAV,
                                       tmp_data_CSCV$iCSCV,
                                       tmp_data_CSNV$iCSNV,
                                       tmp_data_CSAV$dCSAV,
                                       tmp_data_CSCV$dCSCV,
                                       tmp_data_CSNV$dCSNV,
                                       tmp_data_CSAV$CSAV,
                                       tmp_data_CSCV$CSCV,
                                       tmp_data_CSNV$CSNV
)
colnames(data)<-c("nTargets","ctrl","PTB","DTN","DTB","POB","iCSAV","iCSCV","iCSNV","dCSAV","dCSCV","dCSNV","CSAV","CSCV","CSNV")
rownames(data)<-data$nTargets
#####################################################################



# Figure 1d
Fig_1d <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  grid(lwd = 0.3)
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)

  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line

  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)

  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm, cutoff = Inf, n.max=Inf)$rstudent))))
  tmp_matrix = cbind(as.vector(outlier_list),as.vector(data$nTargets[outlier_list]),as.vector(outlierTest(dep_ind_var.lm, cutoff = Inf, n.max=Inf)$bonf.p))
  print(tmp_matrix[tmp_matrix[,2]=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"])
  #print(tmp_matrix)

  outlierTest(dep_ind_var.lm,n.max=3)
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm,n.max=3)$rstudent))))
  points(tmp_ind_var[outlier_list], tmp_dep_var[outlier_list], col="black",pch=4)
  print(cbind(as.vector(data$nTargets[outlier_list]),as.vector(tmp_ind_var[outlier_list]),as.vector(tmp_dep_var[outlier_list])))

  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = "1"
    if(tmp_p < 0.05){
      tmp_tmp_p = "0.05*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "0.01**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "0.001***"
        }
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) < Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    "topleft",
    legend=rp, bty='n', text.col = "black")
  legend("bottomright",
         legend = c("6 vocal learners","Random control sets", "Outliers"),
         col=c("red","grey","black"),
         pch=c(15,20,4),
         text.col = c("red","grey","black")
  )
}
#par(mfrow=c(1,1),mar = c(2, 2, 2, 2))

pdf("./Fig1d_iCSAV_dCSAV_randomctrl_v20201011.pdf",width=6,height=6)
layout(mat = matrix(c(1),1,1,byrow=T))
Fig_1d(data$dCSAV,data$iCSAV,"iCSAV/dCSAV of random controls (n=1,001)", "Different amino acid convergences (dCSAV)", "Identical amino acid convergences (iCSAV)")
dev.off()
#


# Figure 1e
Fig_1e <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  grid(lwd = 0.3)
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)

  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line

  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)

  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm, cutoff = Inf, n.max=Inf)$rstudent))))
  tmp_matrix = cbind(as.vector(outlier_list),as.vector(data$nTargets[outlier_list]),as.vector(outlierTest(dep_ind_var.lm, cutoff = Inf, n.max=Inf)$bonf.p))
  print(tmp_matrix[tmp_matrix[,2]=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"])
  #print(tmp_matrix)

  outlierTest(dep_ind_var.lm,n.max=3)
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm,n.max=3)$rstudent))))
  points(tmp_ind_var[outlier_list], tmp_dep_var[outlier_list], col="black",pch=4)
  print(cbind(as.vector(data$nTargets[outlier_list]),as.vector(tmp_ind_var[outlier_list]),as.vector(tmp_dep_var[outlier_list])))

  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = "1"
    if(tmp_p < 0.05){
      tmp_tmp_p = "0.05*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "0.01**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "0.001***"
        }
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) < Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    "topleft",
    legend=rp, bty='n', text.col = "black")
  legend("bottomright",
         legend = c("6 vocal learners","Random control sets", "Outliers"),
         col=c("red","grey","black"),
         pch=c(15,20,4),
         text.col = c("red","grey","black")
  )
}
#par(mfrow=c(1,1),mar = c(2, 2, 2, 2))

pdf("./Fig1e_iCSAV_POB_randomctrl_v20201011.pdf",width=6,height=6)
layout(mat = matrix(c(1),1,1,byrow=T))
Fig_1e(data$POB,data$iCSAV, "iCSAV/POB of random controls (n=1,001)", "Product of origin branch lengths (POB)", "Identical amino acid convergences (iCSAV)")
dev.off()

# Core control sets
############ Data loading #######################################

data <- read.table("61ctrl/0.summary_Phy_CSAV_20201011.txt",sep='\t',header = T)
attach(data)
tmp_phylo_info = cbind.data.frame(ProTerBL,DistanceTN,DistanceTB,ProOriBL)
colnames(tmp_phylo_info)<-c("PTB","DTN","DTB","POB")
tmp_data_CSAV = cbind.data.frame(nTargets,ctrl,TV_A, TV_C, TV_D)
colnames(tmp_data_CSAV)<-c("nTargets","ctrl","CSAV","iCSAV","dCSAV")
detach(data)

data <- read.table("61ctrl/0.summary_Phy_CSCV_20201011.txt",sep='\t',header = T)
attach(data)
tmp_data_CSCV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_CSCV)<-c("iCSCV","dCSCV","CSCV")
detach(data)

data <- read.table("61ctrl/0.summary_Phy_CSNV_20201011.txt",sep='\t',header = T)
attach(data)
tmp_data_CSNV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_CSNV)<-c("iCSNV","dCSNV","CSNV")
detach(data)

data <- cbind.data.frame(tmp_data_CSAV$nTargets,
                         tmp_data_CSAV$ctrl,
                         tmp_phylo_info,
                         tmp_data_CSAV$iCSAV,
                         tmp_data_CSCV$iCSCV,
                         tmp_data_CSNV$iCSNV,
                         tmp_data_CSAV$dCSAV,
                         tmp_data_CSCV$dCSCV,
                         tmp_data_CSNV$dCSNV,
                         tmp_data_CSAV$CSAV,
                         tmp_data_CSCV$CSCV,
                         tmp_data_CSNV$CSNV
)
colnames(data)<-c("nTargets","ctrl","PTB","DTN","DTB","POB","iCSAV","iCSCV","iCSNV","dCSAV","dCSCV","dCSNV","CSAV","CSCV","CSNV")
rownames(data)<-data$nTargets
#################################################################

# Fig.1f
Fig_1f <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  grid(lwd = 0.3)
  # control sets
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)

  # regression analysis
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line

  # set of birds of prey
  points(tmp_ind_var[data$nTargets=="FALPE,CARCR,TYTAL,HALLE,HALAL,CATAU"],
         tmp_dep_var[data$nTargets=="FALPE,CARCR,TYTAL,HALLE,HALAL,CATAU"],
         col="purple",pch=17)

  # set of waterbirds
  points(tmp_ind_var[data$nTargets=="PELCR,EGRGA,NIPNI,PHACA,FULGL,PYGAD,APTFO,GAVST,PHALE,EURHE,CHAVO,BALRE,PHORU,PODCR,ANAPL"],
         tmp_dep_var[data$nTargets=="PELCR,EGRGA,NIPNI,PHACA,FULGL,PYGAD,APTFO,GAVST,PHALE,EURHE,CHAVO,BALRE,PHORU,PODCR,ANAPL"],
         col="blue",pch=16)

  # set of closest control set
  points(tmp_ind_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="orange",pch=18)

  # set of vocal learners
  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)

  # outliers
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm, cutoff = Inf, n.max=Inf)$rstudent))))
  tmp_matrix = cbind(as.vector(outlier_list),as.vector(data$nTargets[outlier_list]),as.vector(outlierTest(dep_ind_var.lm, cutoff = Inf, n.max=Inf)$bonf.p))
  print(tmp_matrix[tmp_matrix[,2]=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"])
  #print(tmp_matrix)

  outlierTest(dep_ind_var.lm,n.max=3)
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm,n.max=3)$rstudent))))
  points(tmp_ind_var[outlier_list], tmp_dep_var[outlier_list], col="black",pch=4)
  print(cbind(as.vector(data$nTargets[outlier_list]),as.vector(tmp_ind_var[outlier_list]),as.vector(tmp_dep_var[outlier_list])))

  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = "1"
    if(tmp_p < 0.05){
      tmp_tmp_p = "0.05*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "0.01**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "0.001***"
        }
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) < Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    "bottomright",
    legend=rp, bty='n', text.col = "black")
  legend("topleft",
         legend = c("6 vocal learners","5 vocal learners+Swift","6 birds of prey","15 waterbirds","Core control sets", "Outliers"),
         col=c("red","orange","purple","blue","grey","black"),
         pch=c(15,18,17,16,20,4),
         text.col = c("red","orange","purple","blue","grey","black")
  )
}
#par(mfrow=c(1,1),mar = c(2, 2, 2, 2))

pdf("./Fig1f_iCSAV_dCSAV_corectrl_v20201011.pdf",width=6,height=6)
layout(mat = matrix(c(1),1,1,byrow=T))
Fig_1f(data$dCSAV,data$iCSAV,"iCSAV/dCSAV of core controls (n=62)", "Different amino acid convergences (dCSAV)", "Identical amino acid convergences (iCSAV)")
dev.off()

# Fig.1g
Fig_1g <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  grid(lwd = 0.3)
  # control sets
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)

  # regression analysis
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line

  # set of birds of prey
  points(tmp_ind_var[data$nTargets=="FALPE,CARCR,TYTAL,HALLE,HALAL,CATAU"],
         tmp_dep_var[data$nTargets=="FALPE,CARCR,TYTAL,HALLE,HALAL,CATAU"],
         col="purple",pch=17)

  # set of waterbirds
  points(tmp_ind_var[data$nTargets=="PELCR,EGRGA,NIPNI,PHACA,FULGL,PYGAD,APTFO,GAVST,PHALE,EURHE,CHAVO,BALRE,PHORU,PODCR,ANAPL"],
         tmp_dep_var[data$nTargets=="PELCR,EGRGA,NIPNI,PHACA,FULGL,PYGAD,APTFO,GAVST,PHALE,EURHE,CHAVO,BALRE,PHORU,PODCR,ANAPL"],
         col="blue",pch=16)

  # set of closest control set
  points(tmp_ind_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="orange",pch=18)

  # set of vocal learners
  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)

  # outliers
  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm, cutoff = Inf, n.max=Inf)$rstudent))))
  tmp_matrix = cbind(as.vector(outlier_list),as.vector(data$nTargets[outlier_list]),as.vector(outlierTest(dep_ind_var.lm, cutoff = Inf, n.max=Inf)$bonf.p))
  print(tmp_matrix[tmp_matrix[,2]=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"])
  #print(tmp_matrix)


  outlierTest(dep_ind_var.lm,n.max=3)

  outlier_list <- unlist(as.numeric(names(unlist(outlierTest(dep_ind_var.lm,n.max=3)$rstudent))))
  points(tmp_ind_var[outlier_list], tmp_dep_var[outlier_list], col="black",pch=4)
  print(cbind(as.vector(data$nTargets[outlier_list]),as.vector(tmp_ind_var[outlier_list]),as.vector(tmp_dep_var[outlier_list])))
  #print(cbind(as.vector(data$nTargets),as.vector(tmp_ind_var),as.vector(tmp_dep_var)))

  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = "1"
    if(tmp_p < 0.05){
      tmp_tmp_p = "0.05*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "0.01**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "0.001***"
        }
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) < Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    "bottomright",
    legend=rp, bty='n', text.col = "black")
  legend("topleft",
         legend = c("6 vocal learners","5 vocal learners+Swift","6 birds of prey","15 waterbirds","Core control sets", "Outliers"),
         col=c("red","orange","purple","blue","grey","black"),
         pch=c(15,18,17,16,20,4),
         text.col = c("red","orange","purple","blue","grey","black")
  )
}

pdf("./Fig1g_iCSAV_POB_corectrl_v20201011.pdf",width=6,height=6)
layout(mat = matrix(c(1),1,1,byrow=T))
Fig_1g(data$POB,data$iCSAV,"iCSAV/POB of core controls (n=62)", "Product of origin branch lengths (POB)", "Identical amino acid convergences (iCSAV)")
dev.off()
