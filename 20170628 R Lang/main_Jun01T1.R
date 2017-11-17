# R Program for Task 1 of Jun 01, 2017
# edit in Jun, 28, 2017
# June 1,2017
# May 22-24, 2017
# 20170411 PM
# Mar18,2017
# Task #20170221
# Related to May16/Oct17,2016
# edited in Feb21,2017
rm(list=ls())
setwd("/Users/huwei/Dropbox/temp/20170628 R Lang")
#path0<-getwd()
#if (!is.null(path0)) setwd(path0)
# install.packages("nleqslv")
library("nleqslv")
# use data frame to save and call parameters value
value = c(1.0,0.2,0.6,0.06,0.72,0.1,0.36,0.012,0.06,0.53,1.355,18.0,25.0,0.08)
orgPara = data.frame(value)
row.names(orgPara)=c("pstar","b","phi","sigma","beta","lambda","c_f","r","c_p","delta","A","B1","B2","mu")
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
########
fun_q_theta <- function(myPara, x) {
   myPara["A","value"]* x **(-myPara["B1","value"]/myPara["B2","value"])
}
## if you want the so-called 'error function'
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
## (see Abramowitz and Stegun 29.2.29)
## and the so-called 'complementary error function'

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
###############################
for (cc in v_cc){
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
  
  epsilon_d <-vector()
  epsilon_c <-vector()
  theta <-vector()
  alpha <-vector()
  error<-vector()
  flag<-vector()
  
  for (i in 1:length(vt)) {
    myPara<-orgPara
    myPara[casecode,"value"]<-vt[i]
    # define some functions
    # define F(x)
    fun_F_x <- function(x) {
      if (typen == "O") {
        temp1 = (pi/4.0)/sin(pi/4.0) - x
        temp2 = (pi/2.0)/sin(pi/2.0) - ((pi/4.0)/sin(pi/4.0))**2.0
        temp3 = 1.0 + (temp1/sqrt(temp2))**(-4.0)
        return(1.0 - (temp3)**(-1.0))
      } else if (typen == "I"){
        if(x < -1.0*epsilon_u){
          return(0)
        } else if (x > epsilon_u){
          return(1.0)
        } else {
          return((x + epsilon_u)/(2.0*epsilon_u))}
      } else if (typen == "II") {
        temp=sqrt(1.0/pi)-x * sqrt(.5-1.0/pi)
        return(1.0-erf(temp))
      } else if (typen == 'III') {
        return(.5-.5 * erf((log(-x+1.0)+.5*log(2.0))/sqrt(2.0*log(2.0))))
      } else if (typen == "IV") {
        return (((3.0-2.0*x)/np.sqrt(3.0))**(-3.0))
      }
    }
    # integrate of F(x)
    fun_int_F<-function(a,b) {
      integrand <- function(x) {1.0-fun_F_x(x)}
      temp = integrate(integrand,lower = a, upper = b)
      return(temp$value)
    }
    
    # Nonlinear Equations to be solved
    Equation<-function(x){
      y<-numeric(length(x))
      
      epsilon_d <- x[1]
      epsilon_c <- x[2]
      theta <- x[3]
      alpha <- x[4] 
      
      q_theta = fun_q_theta(myPara=myPara,x=theta)
      
      int_Fdu = fun_int_F(epsilon_d,epsilon_u)
      int_Fcu = fun_int_F(epsilon_c,epsilon_u)
      
      part11 = epsilon_d+myPara["lambda","value"]*int_Fdu/(myPara["r","value"]+myPara["lambda","value"]+theta*q_theta)
      part12 = (myPara["b","value"]-myPara["pstar","value"])/myPara["sigma","value"]
      
      temp211 = myPara['delta',"value"]*(myPara['r',"value"]+myPara['lambda',"value"]-myPara['phi',"value"]*(1-fun_F_x(epsilon_c)))*(epsilon_c-epsilon_d)
      temp212 = myPara['r',"value"]+myPara['lambda',"value"]+theta*q_theta
      temp213 = (myPara['lambda',"value"]/(myPara['r',"value"]+myPara['lambda',"value"])-myPara['lambda',"value"]*myPara['delta',"value"]/temp212)*int_Fcu
      part21 = epsilon_c-temp211/((1-myPara['phi',"value"])*temp212)+temp213
      
      temp2211 = ((myPara['beta',"value"]+myPara['phi',"value"]*(1.0-myPara['beta',"value"]))*myPara['c_f',"value"]*(1.0-alpha))/((1-myPara['beta',"value"])*(1-myPara['phi',"value"]))
      temp2212 = (myPara['beta',"value"]*myPara['c_p',"value"]*alpha)/(1-myPara['beta',"value"])
      temp221 = theta*(temp2211 + temp2212)
      part22 = myPara['delta',"value"]*epsilon_d+((1-myPara['delta',"value"])*(myPara['b',"value"]-myPara['pstar',"value"])+temp221)/myPara['delta',"value"]
      
      part31 = 1/q_theta
      
      temp321 = (1-myPara['beta',"value"])*(1-myPara['phi',"value"])/myPara['c_f',"value"]
      temp322 = myPara['sigma',"value"]*(epsilon_u-epsilon_c)/(myPara['r',"value"]+myPara['lambda',"value"])+myPara['delta',"value"]*myPara['sigma',"value"]*(epsilon_c-epsilon_d)/((1-myPara['phi',"value"])*(myPara['r',"value"]+myPara['lambda',"value"]+theta*q_theta))
      part32 = temp321*temp322
      
      part41 = 1/q_theta
      
      temp421 = (1-myPara['beta',"value"])*myPara['delta',"value"]*myPara['sigma',"value"]*(epsilon_u-epsilon_d)
      temp422 = myPara['c_p',"value"]*(myPara['r',"value"]+myPara['lambda',"value"]+theta*q_theta)
      part42 = temp421/temp422
      
      y[1]<-part11-part12
      y[2]<-part21-part22
      y[3]<-part31-part32
      y[4]<-part41-part42
      
      y
    }
    # operator to solve the nonlinear equations
    
    # A numeric vector with an initial guess of the root of the function.
    xstart<-c(-5.0,-1.0,1.0,0.9)
    # trace: trace :Non-negative integer. A value of 1 will give a detailed report of the progress of the iteration.
    
    sol<-nleqslv(xstart,Equation,method="Broyden",control = list(trace=1))
    
    epsilon_d[i]<-sol$x[1]
    epsilon_c[i]<-sol$x[2]
    theta[i]<-sol$x[3]
    alpha[i]<-sol$x[4]
    
    # Euclidean norm of function values
    error[i]<-norm(sol$fvec,type="2")
    # termination code as integer. The values returned are:
    # 1
    # Function criterion is near zero. Convergence of function values has been achieved.
    # 2
    # x-values within tolerance. This means that the relative distance between two consecutive x-values is smaller than xtol but that the function value criterion is still larger than ftol. Function values may not be near zero; therefore the user must check if function values are acceptably small.
    # 3
    # No better point found. This means that the algorithm has stalled and cannot find an acceptable new point. This may or may not indicate acceptably small function values.
    # 4
    # Iteration limit maxit exceeded.
    # 5
    # Jacobian is too ill-conditioned.
    # 6
    # Jacobian is singular.
    # 7
    # Jacobian is unusable.
    # -10
    # User supplied Jacobian is most likely incorrect.
    flag[i]<- sol$termcd
  }
  # Data Outpu
  write.table(epsilon_d,file=paste(str_fig,str_para0,"epsilon_d.txt",sep="_"),sep="\n",row.name=FALSE,col.names = FALSE)
  # temp<-read.table(file=paste(str_fig,str_para0,"epsilon_d.txt",sep="_"),sep="\n")
  # epsilon_d = temp[[1]]
  write.table(epsilon_c,file=paste(str_fig,str_para0,"epsilon_c.txt",sep="_"),sep="\n",row.name=FALSE,col.names = FALSE)
  write.table(theta,file=paste(str_fig,str_para0,"theta.txt",sep="_"),sep="\n",row.name=FALSE,col.names = FALSE)
  write.table(alpha,file=paste(str_fig,str_para0,"alpha.txt",sep="_"),sep="\n",row.name=FALSE,col.names = FALSE)
  write.table(error,file=paste(str_fig,str_para0,"error.txt",sep="_"),sep="\n",row.name=FALSE,col.names = FALSE)
  write.table(flag,file=paste(str_fig,str_para0,"flag.txt",sep="_"),sep="\n",row.name=FALSE,col.names = FALSE)
}