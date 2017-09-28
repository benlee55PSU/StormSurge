###################################################################################
#
# This file generates samples from a non-stationary GEV distribution under 4 RCP scenarios. 
# WE use the parameter estimates from Grinsted et al. (2013)
# THe fExtremes package is used to generate the GEV data. 
######################################################################################
# Step 1: Load Temperature anomaly data for past and future 
load("Box Sync/FinalPaper/DataRepository/GEVPostiveControlHistoricalRCP.RData")
## Sample Generation for Positive Control
# Initialize
library(fExtremes)
#Generate Positive Control Sample
set.seed(11) # Set Seed
temps_26<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,2],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # Dataframe for sample

temps_45<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,3],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # Dataframe for sample

temps_60<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,4],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # Dataframe for sample

temps_85<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,5],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # Dataframe for sample


mu0<-2.41 ;sigm0<-0.48 ;xi0<-0.54 ;amu<-0.13;asigma<-0.49; axi<-0.22 # Initialize Grinsted Parameters
trueparmat<-matrix(c(mu0,sigm0,xi0,amu,asigma,axi),nrow=1,ncol=6) # True Parameters
colnames(trueparmat)<-c("mu_0_true","sigma_0_true","xi_0_true","amu_true","asigma_true","axi_true")

# Create Containers for Sample Data
samplemat_26<-matrix(NA,nrow=nrow(temps_26),30)# Matrix to hold sample data
samplemat_45<-matrix(NA,nrow=nrow(temps_45),30)# Matrix to hold sample data
samplemat_60<-matrix(NA,nrow=nrow(temps_60),30)# Matrix to hold sample data
samplemat_85<-matrix(NA,nrow=nrow(temps_85),30)# Matrix to hold sample data

# Specify GEV Parameters under the different RCP scenarios
for(i in 1:nrow(temps_26)){
  set.seed(i)# Set Seed
  
  # RCP2.6
  temps_26[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_26[i,2])#mu
  temps_26[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_26[i,2]))#sigma
  temps_26[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_26[i,2])#xi
  samplemat_26[i,]<-rgev(30, xi=temps_26[i,5], mu=temps_26[i,3], beta=temps_26[i,4]) # Sample Generation

  # RCP4.5
  temps_45[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_45[i,2])#mu
  temps_45[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_45[i,2]))#sigma
  temps_45[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_45[i,2])#xi
  samplemat_45[i,]<-rgev(30, xi=temps_45[i,5], mu=temps_45[i,3], beta=temps_45[i,4]) # Sample Generation
  
  # RCP6.0
  temps_60[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_60[i,2])#mu
  temps_60[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_60[i,2]))#sigma
  temps_60[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_60[i,2])#xi
  samplemat_60[i,]<-rgev(30, xi=temps_60[i,5], mu=temps_60[i,3], beta=temps_60[i,4]) # Sample Generation
  
  # RCP8.5
  temps_85[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_85[i,2])#mu
  temps_85[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_85[i,2]))#sigma
  temps_85[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_85[i,2])#xi
  samplemat_85[i,]<-rgev(30, xi=temps_85[i,5], mu=temps_85[i,3], beta=temps_85[i,4]) # Sample Generation
}

temps_26<-cbind(temps_26,samplemat_26) # Create master data frame
temps_45<-cbind(temps_45,samplemat_45) # Create master data frame
temps_60<-cbind(temps_60,samplemat_60) # Create master data frame
temps_85<-cbind(temps_85,samplemat_85) # Create master data frame


for (i in 1:30){
  colnames(temps_26)[i+5]<-paste("S",i,sep="") # Name each column (sample)
  colnames(temps_45)[i+5]<-paste("S",i,sep="") # Name each column (sample)
  colnames(temps_60)[i+5]<-paste("S",i,sep="") # Name each column (sample)
  colnames(temps_85)[i+5]<-paste("S",i,sep="") # Name each column (sample)
}

#Plot Time series of extreme events under each RCP Scenario
plot.ts(temps_26[,6],xlim=c(0,10000),ylim=c(1,1000)) # Time series of events under RCP 2.6
plot.ts(temps_45[,8],xlim=c(0,10000),ylim=c(1,1000)) # Time series of events under RCP 4.5
plot.ts(temps_60[,10],xlim=c(0,10000),ylim=c(1,1000)) # Time series of events under RCP 6.0
plot.ts(temps_85[,12],xlim=c(0,10000),ylim=c(1,1000)) # Time series of events under RCP 8.5

# Save resulting file
save(temps_26,temps_45,temps_60,temps_85,trueparmat,
     file="Box Sync/FinalPaper/DataRepository/PositiveControlSample_Learning.Rdata") #Save
rm(list=ls())


# Save in list form for RCP4.5
load("Box Sync/FinalPaper/DataRepository/PositiveControlSample_Learning.Rdata")
tempdat<-list(tempset=temps_45$tempset1, #Global Average Temperature
              year=temps_45$year, # YEar
              datset1=temps_45$S1, # Non-stationary Sample 1
              datset2=temps_45$S2, # Non-stationary Sample 2
              datset3=temps_45$S3, # Non-stationary Sample 3
              datset4=temps_45$S4, # Non-stationary Sample 4
              datset5=temps_45$S5, # Non-stationary Sample 5
              datset6=temps_45$S6, # Non-stationary Sample 6
              datset7=temps_45$S7, # Non-stationary Sample 7
              datset8=temps_45$S8, # Non-stationary Sample 8
              datset9=temps_45$S9, # Non-stationary Sample 9
              datset10=temps_45$S10) # Non-stationary Sample 10

# Save Final File
save(tempdat,file="Box Sync/FinalPaper/DataRepository/PositiveControlSample45.Rdata")
# Load in Final File
load("Box Sync/FinalPaper/DataRepository/PositiveControlSample45.Rdata")

# Specify stop points for each of the learning sample
end<-c(seq(10,210,by=10)*52,11544) # 10 year intervals
# Create the final Learning Sample under the RCP4.5 Scenario
tempdatlist<-list()
for(i in 1:length(end)){
  tempdatlist[[i]]<-temps_45[1:end[i],]
}
lapply(tempdatlist,dim) # Check Dimensions of each sample

finaltempdat45<-list()# FInal Learning Sample
for(i in 1:length(tempdatlist)){
  finaltempdat45[[i]]<-list(tempset=tempdatlist[[i]]$tempset1,
                            year=tempdatlist[[i]]$year,
                            datset1=tempdatlist[[i]]$S1)
}
# Save the final Learning Sample
save(finaltempdat45,file="Box Sync/FinalPaper/DataRepository/FinalPositiveControlSample45.Rdata")
