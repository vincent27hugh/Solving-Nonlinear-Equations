# main_Jun01T1_Plot.py

# Python PLOT program for Task 1 of Jun 01,2017
# June 1,2017
# May 22-24, 2017
# 20170411 PM
# Mar18,2017
# Task #20170221
# Related to May16/Oct17,2016
# edited in Feb21,2017
rm(list=ls())
setwd("/Users/huwei/Dropbox/temp/20170628 R Lang")
# path0<-getwd()
# if (!is.null(path0)) setwd(path0)
# use data frame to save and call parameters value
value = c(1.0,0.2,0.6,0.06,0.72,0.1,0.36,0.012,0.06,0.53,1.355,18.0,25.0,0.08)
myPara = data.frame(value)
row.names(myPara)=c("pstar","b","phi","sigma","beta","lambda","c_f","r","c_p","delta","A","B1","B2","mu")
#         value
# pstar   1.000
# b       0.200
# phi     0.600
# sigma   0.060
# beta    0.720
# lambda  0.100
# c_f     0.360
# r       0.012
# c_p     0.060
# delta   0.530
# A       1.355
# B1     18.000
# B2     25.000
# mu      0.080
# call the value:
# myPara["pstar","value"]

# variable vector for different cases
lvt <- list(pstar = seq(0.18,2.5,by=0.01),
            b = seq(0.01,1.0,by=0.01),
            phi = seq(0.3,0.95,by=0.01),
            sigma = seq(0.01,1.2,by=0.01),
            beta = seq(0.4,0.99,by=0.01),
            lambda = seq(0.02,0.63,by=0.01),
            c_f = seq(0.12,1.3,by=0.01),
            r = seq(0.003,0.15,by=0.001),
            c_p = seq(0.01,0.2,by=0.001),
            delta = seq(0.21,0.7,by=0.01))
# extract part of object
# example 
# test = "delta"
# vt = lvt[[test]]
lstr_var <- list(pstar = "p^*",
            b = "b",
            phi = "$\\phi$",
            sigma = "$\\sigma$",
            beta = "$\\beta$",
            lambda = "$\\lambda$",
            c_f = "$c_f$",
            r = "r",
            c_p = "$c_p$",
            delta = "$\\delta$")
# extract part of object
# example 
# test = "delta"
# vt = lvt[[test]]

#######################
taskcode<-"Jun01T1"
num_para<-"V27"
str_para0<-paste("paraJun01",num_para,sep = "_")

cell_typen<- c("I","II","III","IV","O")
# 10 cases
cell_case<-c("pstar","b","phi","sigma","beta","lambda","c_f","r","c_p","delta")

v_tt <- 3 # 1:length(cell_typen)
v_cc <- (10)

tt = v_tt
typen=cell_typen[tt]

for (cc in v_cc) {
  casecode=cell_case[cc]
  
  if (typen == "I") {
    epsilon_u <- sqrt(3)
  } else if (typen == "II") {
    epsilon_u <- sqrt(2)/(sqrt(pi-2))
  } else if (typen == "III") {
    epsilon_u <- 1
  } else if (typen == "IV") {
    epsilon_u <- (3-sqrt(3))/2
  } else if (typen == "O") {
    epsilon_u <- (pi/4)/(sin(pi/4))
  }
  
  ### variable vector
  vt = lvt[[casecode]]
  
  str_fig = paste(taskcode,typen,casecode,sep="_")
  
  # Data Outpu
  temp<-read.table(file=paste(str_fig,str_para0,"epsilon_d.txt",sep="_"),sep="\n")
  epsilon_d = temp[[1]]
  
  temp<-read.table(file=paste(str_fig,str_para0,"epsilon_c.txt",sep="_"),sep="\n")
  epsilon_c = temp[[1]]
  
  temp<-read.table(file=paste(str_fig,str_para0,"theta.txt",sep="_"),sep="\n")
  theta = temp[[1]]
  
  temp<-read.table(file=paste(str_fig,str_para0,"alpha.txt",sep="_"),sep="\n")
  alpha = temp[[1]]
  
  temp<-read.table(file=paste(str_fig,str_para0,"error.txt",sep="_"),sep="\n")
  error = temp[[1]]
  
  temp<-read.table(file=paste(str_fig,str_para0,"flag.txt",sep="_"),sep="\n")
  flag = temp[[1]]
  
  # you need to install package:
  # install.packages("latex2exp")
  library("latex2exp")
  #
  flag1<-flag==1
  flag2<-alpha>0&alpha<1
  flag3<-theta>0
  pflag<-flag1&flag2&flag3
  #
  png(filename=paste(str_fig,str_para0,"sol_flag.png",sep="_"),width=800,height=700)
  par(mfrow=c(2,2))
  plot(vt[pflag],epsilon_d[pflag],type="l",lwd=3,
       xlab=list(TeX(lstr_var[casecode]),cex=1.5),
       ylab=list(TeX("$\\epsilon_d$"),cex=1.5),
       main=list(paste(taskcode,typen,sep=","),cex=1.75),
       xlim=c(min(vt[pflag]),max(vt[pflag])))
  plot(vt[pflag],epsilon_c[pflag],type="l",lwd=3,
       xlab=list(TeX(lstr_var[casecode]),cex=1.5),
       ylab=list(TeX("$\\epsilon_c$"),cex=1.5),
       main=list(paste(taskcode,typen,sep=","),cex=1.75),
       xlim=c(min(vt[pflag]),max(vt[pflag])))
  plot(vt[pflag],theta[pflag],type="l",lwd=3,
       xlab=list(TeX(lstr_var[casecode]),cex=1.5),
       ylab=list(TeX("$\\theta$"),cex=1.5),
       main=list(paste(taskcode,typen,sep=","),cex=1.75),
       xlim=c(min(vt[pflag]),max(vt[pflag])))
  plot(vt[pflag],alpha[pflag],type="l",lwd=3,
       xlab=list(TeX(lstr_var[casecode]),cex=1.5),
       ylab=list(TeX("$\\alpha$"),cex=1.5),
       main=list(paste(taskcode,typen,sep=","),cex=1.75),
       xlim=c(min(vt[pflag]),max(vt[pflag])))
  dev.off()
}