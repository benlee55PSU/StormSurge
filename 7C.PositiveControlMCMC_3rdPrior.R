#######################################################################################
# Bayesian Inference via MCMC
# Positive Control Sample
# Fit a GEV distribution on data generated from a fully non-stationary GEV distribution
# Using Narrow Normal Prior Distributions on the Parameters
# (1) Stationary
# (2) Mu Only Non-stationary
# (3) Mu Sigma Non-Stationary
# (4) Fully Non-Stationary
###########################################################################################
#Initialize Data
rm(list=ls())
setwd("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/")
load("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/FinalPositiveControlSample45.Rdata")
maindat<-finaltempdat45[[14]]
blockmaxdat<-cbind(maindat$year ,maindat$tempset ,maindat$datset1 )
blockmaxdat<-blockmaxdat[blockmaxdat[,1]%in%1923:2010,]
datset=blockmaxdat[,3] # Set Dataset
tempset=blockmaxdat[,2] # Set temperature Set

###########################################################################################
# (1) Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/Prior2SourceStat.R")
# Initial Conditions
start<-c(rep(1,3)) # Start Value
errvect<-c(0.0625, 0.03125, 0.06255)
uvect<-rep(0,3) # Does not matter
sdvect<-rep(1,3)  # Does not matter

#Testing for resuming MCMC
statrun1<-MCMC_breaks_resume(10000,2500,res=FALSE,resdat=NA,datset,tempset,
                             start,errvect,uvect,sdvect,filetitle="Test1") # First Run
statrun2<-MCMC_breaks_resume(10000,2500,res=TRUE,resdat=statrun1,datset,tempset,
                             start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Run Diagnostics Again
BM_MCMC(statrun2,int=500,burn=0,bline=6000)
Mean_MCMC(statrun2,int=500,burn=0,bline=6000)

burnin=1000 # Burnin
MULTnsgevplots(statrun2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(statrun2,burn=burnin) #RejectionRate
CredIntervalsGEV(statrun2,burn=burnin) # Credible Intervals
save(statrun2,datset,tempset,file="Prior3_PositiveControl_Stat_Output.RData")


###########################################################################################
# (2) Mu Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/Prior2SourceMu.R")
# Initial Conditions
start<-c(rep(1,4)) # Start Value
errvect<-c(0.0625, 0.03125 ,0.0625, 0.0625)
uvect<-c(rep(0,4)) # Does not matter
sdvect<-c(rep(1,4))  # Does not matter
#Testing for resuming MCMC
murun1<-MCMC_breaks_resume(10000,2500,res=FALSE,resdat=NA,datset,tempset,
                           start,errvect,uvect,sdvect,filetitle="Test1") # First Run
murun2<-MCMC_breaks_resume(10000,2500,res=TRUE,resdat=murun1,datset,tempset,
                           start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Run Diagnostics Again
BM_MCMC(murun2,int=500,burn=0,bline=6000)
Mean_MCMC(murun2,int=500,burn=0,bline=6000)
burnin=1000 # Burnin
MULTnsgevplots(murun2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(murun2,burn=burnin) #RejectionRate
CredIntervalsGEV(murun2,burn=burnin) # Credible Intervals
save(murun2,datset,tempset,file="Prior3_PositiveControl_Mu_Output.RData")


###########################################################################################
# (3) Mu Sigma Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/Prior2SourceMuSigma.R")
# Initial Conditions
start<-c(rep(1,5)) # Start Value
errvect<-c(0.0625, 0.03125, 0.0625,0.0625, 0.25)
uvect<-c(rep(0,5)) # Does not matter
sdvect<-c(rep(1,5))  # Does not matter
#Testing for resuming MCMC
run1<-MCMC_breaks_resume(10000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(10000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Run Diagnostics Again
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=1000 # Burnin
MULTnsgevplots(run2,rm.burn=FALSE,burn=burnin) #Trace Plots Without Burn-In
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,datset,tempset,file="Prior3_PositiveControl_MuSigma_Output.RData")
###########################################################################################
# (4) Fully Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/Prior2Source.R")
# Initial Conditions
start<-c(rep(1,6)) # Start Value
errvect<-c(0.0625, 0.03125, 0.0625,0.0625, 0.25 ,0.25)
uvect<-c(rep(0,6)) # Does not matter
sdvect<-c(rep(1,6))  # Does not matter
#Testing for resuming MCMC
run1<-MCMC_breaks_resume(10000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(10000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Run Diagnostics Again
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=1000 # Burnin
MULTnsgevplots(run2,rm.burn=FALSE,burn=burnin) #Trace Plots Without Burn-In
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,datset,tempset,file="Prior3_PositiveControl_FullNS_Output.RData")
