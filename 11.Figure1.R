######################################################################################
# Code to generate Figure 1
# Create Return Level Plots for the 4 GEV models
#
######################################################################################
### Setup 
rm(list=ls())
# Load Grinsted Storm Surge Data
load("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/Grinstedmaxblockdata.RData") # Load Grinsted Storm Surge Index
# Load MCMC Output for each positive Control sample
setwd("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/MCMC_Output")
load("Grinsted_Stat_Output.RData");run2_stat<-run2;rm(run2) # Data from Stationary GEV
load("Grinsted_FullNS_Output.RData");run2_ns<-run2;rm(run2)#   Non-Stationary
load("Grinsted_Mu_Output.RData");run2_mu<-run2;rm(run2) # Mu-only Non-stationary
load("Grinsted_MuSigma_Output.RData");run2_musigma<-run2;rm(run2) # Mu and Sigma Non-Stationary
dataset=newblockmax$surge # Set Grinsted SSI 
tempset=newblockmax$relabase # Set temperature
burnin=6000 # Set Burnin for the MCMC Output

## Determine 90% Credible Interval for Each parameter
testdat<-list() # Initial List
#Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceStat.R")
dat<-run2_stat[[1]][-(1:burnin),] # Load Data
testdat[[1]]<-t(apply(dat,2,quantile,probs=c(0.05,0.16,0.5,0.84,0.95))) # Get 90% Credible Interval
testdat[[1]][,3]<-apply(dat,2,mean) # Replace Median with Mean
#Mu-only Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceMu.R")
dat<-run2_mu[[1]][-(1:burnin),]
testdat[[2]]<-t(apply(dat,2,quantile,probs=c(0.05,0.16,0.5,0.84,0.95)))
testdat[[2]][,3]<-apply(dat,2,mean)# Replace Median with Mean
#Mu Sigma-only Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceMuSigma.R")
dat<-run2_musigma[[1]][-(1:burnin),]
testdat[[3]]<-t(apply(dat,2,quantile,probs=c(0.05,0.16,0.5,0.84,0.95)))
testdat[[3]][,3]<-apply(dat,2,mean)# Replace Median with Mean
#Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSource.R")
dat<-run2_ns[[1]][-(1:burnin),]
testdat[[4]]<-t(apply(dat,2,quantile,probs=c(0.05,0.16,0.5,0.84,0.95)))
testdat[[4]][,3]<-apply(dat,2,mean) # Replace Median with Mean

########################################################################################
# Update Parmeters with respect to temperature. We will calculate the new GEV distribution
# parameters after we factor in temperature data. 

par_stat<-list()# Initial List

###Stationary Model 
t0=0 # Baseline Temperature: Updated Parameters will be at baseline temp. 
mu0<-testdat[[1]][1,] # Mu 
sigma0<-exp(testdat[[1]][2,]) # SIgma 
xi0<-testdat[[1]][3,] # Xi
par_stat[[1]]<-cbind(mu0,sigma0,xi0) # Combine
t1=1 # One degree warmer
mu1<-testdat[[1]][1,]
sigma1<-exp(testdat[[1]][2,])
xi1<-testdat[[1]][3,]
par_stat[[1]]<-cbind(par_stat[[1]],cbind(mu1,sigma1,xi1)) # combine

###Mu-Only Non-Stationary Model 
t0=0 # Baseline Temperature: Updated Parameters will be at baseline temp. 
mu0<-testdat[[2]][1,]*(1+t0*testdat[[2]][4,] ) # Mu 
sigma0<-exp(testdat[[2]][2,]) # Sigma
xi0<-testdat[[2]][3,] # Xi
par_stat[[2]]<-cbind(mu0,sigma0,xi0)  # Combine
t1=1 # Warmer Temperature 
mu1<-testdat[[2]][1,]*(1+t1*testdat[[2]][4,])# Mu 
sigma1<-exp(testdat[[2]][2,])# Sigma
xi1<-testdat[[2]][3,]# Xi
par_stat[[2]]<-cbind(par_stat[[2]],cbind(mu1,sigma1,xi1)) # Combine


###Mu and Sigma-Only Non-Stationary Model 
t0=0# Baseline Temperature: Updated Parameters will be at baseline temp. 
mu0<-testdat[[3]][1,]*(1+t0*testdat[[3]][4,])
sigma0<-exp(testdat[[3]][2,]*(1+t0*testdat[[3]][5,]))
xi0<-testdat[[3]][3,]
par_stat[[3]]<-cbind(mu0,sigma0,xi0)
t1=1# Warmer Temperature 
mu1<-testdat[[3]][1,]*(1+t1*testdat[[3]][4,])
sigma1<-exp(testdat[[3]][2,]*(1+t1*testdat[[3]][5,]))
xi1<-testdat[[3]][3,]
par_stat[[3]]<-cbind(par_stat[[3]],cbind(mu1,sigma1,xi1))

###Non-Stationary
t0=0# Baseline Temperature: Updated Parameters will be at baseline temp. 
mu0<-testdat[[4]][1,]*(1+t0*testdat[[4]][4,])
sigma0<-exp(testdat[[4]][2,]*(1+t0*testdat[[4]][5,]))
xi0<-testdat[[4]][3,]*(1+t0*testdat[[4]][6,])
par_stat[[4]]<-cbind(mu0,sigma0,xi0)
t1=1# Warmer Temperature 
mu1<-testdat[[4]][1,]*(1+t1*testdat[[4]][4,])
sigma1<-exp(testdat[[4]][2,]*(1+t1*testdat[[4]][5,]))
xi1<-testdat[[4]][3,]*(1+t1*testdat[[4]][6,])
par_stat[[4]]<-cbind(par_stat[[4]],cbind(mu1,sigma1,xi1))

############################################################################################
#
# Calculate the Return Periods for a grid of SSI Levels
#

level<-c(1:5000) # Calculate Return Periods for SSI levels 1-400. 
year100<-100*52 # 100 Year Storm. Multiply by weeks
retnew0<-list();retnew1<-list() # Create a list for the baseline and warmer temperatures
#Stationary Model 
retnew0[[1]]<-apply(par_stat[[1]],1,function(x) benreturnyear(mu=x[1],sigma=x[2],xi=x[3],x=level))/52
retnew1[[1]]<-apply(par_stat[[1]],1,function(x) benreturnyear(mu=x[4],sigma=x[5],xi=x[6],x=level))/52
#Mu-Only Non-Stationary Model 
retnew0[[2]]<-apply(par_stat[[2]],1,function(x) benreturnyear(mu=x[1],sigma=x[2],xi=x[3],x=level))/52
retnew1[[2]]<-apply(par_stat[[2]],1,function(x) benreturnyear(mu=x[4],sigma=x[5],xi=x[6],x=level))/52
#Mu and Sigma-Only Non-Stationary Model 
retnew0[[3]]<-apply(par_stat[[3]],1,function(x) benreturnyear(mu=x[1],sigma=x[2],xi=x[3],x=level))/52
retnew1[[3]]<-apply(par_stat[[3]],1,function(x) benreturnyear(mu=x[4],sigma=x[5],xi=x[6],x=level))/52
#Non-Stationary
retnew0[[4]]<-apply(par_stat[[4]],1,function(x) benreturnyear(mu=x[1],sigma=x[2],xi=x[3],x=level))/52
retnew1[[4]]<-apply(par_stat[[4]],1,function(x) benreturnyear(mu=x[4],sigma=x[5],xi=x[6],x=level))/52

# Calculate 100year storm levels
vect100yearstorm<-vector("numeric")
for(i in 1:4){
  vect100yearstorm[i]<-apply(par_stat[[i]],1,function(x) benreturn(mu=x[1],sigma=x[2],xi=x[3],t=year100))[3]
}


############################################################################################
#
# Plot the Return Periods/Levels for each model
#
dev.off()
fontsize<-2.8
par(mfrow=c(2,2),mar=c(5, 6, 2, 1) + 0.1,mgp=c(4, 1.5, 0))# Plotting Parameters

labelvect<-c("(a)","(b)","(c)","(d)")
for(i in 1:4){
  plot(y=level,x=retnew0[[i]][,3],typ="l",lwd=3,col="blue",log=c("xy"),xaxt="n",yaxt="n",
       xlim=c(0.02,500),ylim=c(2,4000),
       ylab="Storm Surge Index",
       xlab="Return period (years)",#main=labelvect[i], 
       font.lab=1 ,cex.lab=fontsize,cex.main=fontsize,cex.axis=fontsize)
  
  axis(2, at=c(2, 5, 10, 25,60,150,500,2000), 
       labels=TRUE,cex.axis=fontsize)
  labels1<-c(0.02, 0.1, 0.3, 1,3,10,30,50,100,250,500)
  axis(1, at=c(0.02, 0.1, 0.3, 1,3,10,30,50,100,250,500),
       labels=labels1,cex.axis=fontsize)

  
  lines(y=level,x=retnew1[[i]][,3],col="red",lwd=3)

  
  polydatx1<-c(retnew0[[i]][,1],rev(retnew0[[i]][,5]))
  polydaty1<-c((level),rev(level))
  polygon(y=polydaty1, x=polydatx1,
          col=rgb(0, 0, 1,0.25), border = NA)
  
  polydatx2<-c(retnew1[[i]][,1],rev(retnew1[[i]][,5]))
  polydaty2<-c((level),rev(level))
  polygon(y=polydaty2, x=polydatx2,
          col=rgb(1, 0, 0,0.25), border = NA)
  abline(h=113.396956,lwd=2,lty=3)
  
  points(x=((length(dataset)+1)/(1:length(dataset)))/52, 
         y=sort(dataset,decreasing=TRUE),pch=4,
         lwd=2,cex=2)
  abline(v=100,lty=1,lwd=2,col="gray")
  abline(h=vect100yearstorm[i],lty=1,lwd=2,col="gray")
  
  mtext(labelvect[i], 3, adj=-0.1,las=1, font=2,cex=1.75)


  if(i==1){
  legend(x=0.007,y=7000, 
         legend = c(
                    'Base Temperature', '90% Cred. Interval (Base)', 
                    '1-Degree Warmer', '90% Cred. Interval (Warmer)',
                    'Observations','Hurricane Katrina','100-Year Event'), 
         bty = "n",
         seg.len=rep(0.6,7),
         lwd = c(3,NA,3,NA,NA,3,3),
         lty = c(1,NA,1,NA,NA,3,1),
         pch=c(NA,NA,NA,NA,4,NA,NA),
         col = c("blue", "blue", "red","red",1,1,"gray"),
         density=c(0,NA,0,NA,0,0,0),
         fill = c("blue", rgb(0, 0, 1,0.25), "red",rgb(1, 0, 0,0.25),1,1,1),
         border = c(NA,NA,NA,NA,NA,NA,NA),
         cex = 1.9, x.intersp = 0.6, y.intersp = 0.75,ncol=2,text.font=1,text.width=rep(1.5,7)
  )
  }
}

# save as 1900 by 1200