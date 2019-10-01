#######################################################################################
# Bayesian Inference via MCMC
# Observed Grinsted Data
# Fit a GEV distribution on data generated from a fully non-stationary GEV distribution
# Using Uniform Prior Distributions on the Parameters
# (1) Stationary
# (2) Mu Only Non-stationary
# (3) Mu Sigma Non-Stationary
# (4) Fully Non-Stationary
###########################################################################################
#Initialize Data
rm(list=ls())
setwd("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/")
load("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/Grinstedmaxblockdata.RData")
datset=newblockmax$surge # Set Dataset
tempset=newblockmax$relabase # Set temperature Set
###########################################################################################
# (1) Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceStat.R")
# Initial Conditions
start<-c(rep(1,3)) # Start Value
errvect<-c(0.0625, 0.03125, 0.06255) # Proposal Distribution SD
uvect<-c(rep(NA,3)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(NA,3))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors
#MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2")
# Diagnostics
# Run Diagnostics Again
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=10000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,file="Grinsted_Stat_Output.RData") # Save Data
load("Grinsted_Stat_Output.RData")
cred.table
###########################################################################################
# (2) Mu Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceMu.R")
# Initial Conditions
start<-c(rep(1,4)) # Start Value
errvect<-c(0.0625 ,0.03125, 0.0625 ,0.0625) # SD for Proposal 
uvect<-c(rep(NA,4)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(NA,4))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors
#MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Diagnostics
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=1000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals

save(run2,file="Grinsted_Mu_Output.RData") # Save Data
load("Grinsted_Mu_Output.RData")
cred.table
###########################################################################################
# (3) Mu Sigma Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceMuSigma.R")
# Initial Conditions
start<-c(rep(1,5)) # Start Value
errvect<-c(0.0625, 0.03125, 0.0625,0.0625, 0.25) # SD for Proposal
uvect<-c(rep(NA,5)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(NA,5))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors

##############################################################################################
#Testing for resuming MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Run Diagnostics Again
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=6000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,file="Grinsted_MuSigma_Output.RData")
load("Grinsted_MuSigma_Output.RData")
cred.table
###########################################################################################
# (4) Fully Non-Stationary
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSource.R") # Load Source File
# Initial Conditions
start<-c(rep(1,6)) # Start Value
errvect<-c(0.0625, 0.03125, 0.0625,0.0625, 0.25 ,0.25) # Proposal SD
uvect<-c(rep(NA,6)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(NA,6))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors

#Testing for resuming MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Run Diagnostics 
burnin=6000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals

#save(run2,file="Grinsted_FullNS_Output.RData") # Save Data
load("Grinsted_FullNS_Output.RData")
cred.table
###########################################################################################

