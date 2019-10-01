############################################################################################
#Generate Positive Control Sample for Learning
#Non-Stationary GEV Data
# Generated Using GLobal Average Temperature
##############################################################################################
# Load Data and restrict to years 1923-to-2100
load("Box Sync/FinalPaper/DataRepository/GEVPostiveControlHistoricalRCP.RData")
finalglobalrcp<-finalglobalrcp[which(finalglobalrcp[,1]%in%1923:2100),]

## Sample Generation for Positive Control
# Initialize
library(fExtremes)

#Generate Positive Control Sample
set.seed(11) # Set Seed
temps_26<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,2],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # RCP2.6

temps_45<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,3],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # RCP4.5

temps_60<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,4],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # RCP6.0

temps_85<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,5],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # RCP8.5

# Parameters from Grinsted et al. (2013)
mu0<-2.41 ;sigm0<-0.48 ;xi0<-0.54 ;amu<-0.13;asigma<-0.49; axi<-0.22 # Initialize Grinsted Parameters
trueparmat<-matrix(c(mu0,sigm0,xi0,amu,asigma,axi),nrow=1,ncol=6) # True Parameters
colnames(trueparmat)<-c("mu_0_true","sigma_0_true","xi_0_true","amu_true","asigma_true","axi_true")

# Create Matrices to hold sample data
samplemat_26<-matrix(NA,nrow=nrow(temps_26),30)# RCP2.6
samplemat_45<-matrix(NA,nrow=nrow(temps_45),30)# RCP4.5
samplemat_60<-matrix(NA,nrow=nrow(temps_60),30)# RCP6.0
samplemat_85<-matrix(NA,nrow=nrow(temps_85),30)# RCP8.5

# Generate Data
for(i in 1:nrow(temps_26)){
  set.seed(i)
  # RCP 2.6
  temps_26[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_26[i,2])#mu
  temps_26[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_26[i,2]))#sigma
  temps_26[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_26[i,2])#xi
  samplemat_26[i,]<-rgev(30, xi=temps_26[i,5], mu=temps_26[i,3], beta=temps_26[i,4]) # Sample Generation
  # RCP 4.5
  temps_45[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_45[i,2])#mu
  temps_45[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_45[i,2]))#sigma
  temps_45[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_45[i,2])#xi
  samplemat_45[i,]<-rgev(30, xi=temps_45[i,5], mu=temps_45[i,3], beta=temps_45[i,4]) # Sample Generation
  # RCP 6.0
  temps_60[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_60[i,2])#mu
  temps_60[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_60[i,2]))#sigma
  temps_60[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_60[i,2])#xi
  samplemat_60[i,]<-rgev(30, xi=temps_60[i,5], mu=temps_60[i,3], beta=temps_60[i,4]) # Sample Generation
  # RCP 8.5
  temps_85[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_85[i,2])#mu
  temps_85[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_85[i,2]))#sigma
  temps_85[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_85[i,2])#xi
  samplemat_85[i,]<-rgev(30, xi=temps_85[i,5], mu=temps_85[i,3], beta=temps_85[i,4]) # Sample Generation
}

# Create master data frame to hold GEV data and temperatures
temps_26<-cbind(temps_26,samplemat_26) # RCP2.6
temps_45<-cbind(temps_45,samplemat_45) # RCP4.5
temps_60<-cbind(temps_60,samplemat_60) # RCP6.0
temps_85<-cbind(temps_85,samplemat_85) # RCP8.5

# Distinguish each sample by naming each column
for (i in 1:30){
  colnames(temps_26)[i+5]<-paste("S",i,sep="") # Name each column (sample)
  colnames(temps_45)[i+5]<-paste("S",i,sep="") # Name each column (sample)
  colnames(temps_60)[i+5]<-paste("S",i,sep="") # Name each column (sample)
  colnames(temps_85)[i+5]<-paste("S",i,sep="") # Name each column (sample)
}
# Save File
save(temps_26,temps_45,temps_60,temps_85,trueparmat,
     file="Box Sync/FinalPaper/DataRepository/PositiveControlSample_Learning1923.Rdata") #Save

rm(list=ls())

#######################################################################################
# Save RCP4.5 Data in list form for Learning Study
# Each item in the list is a learning sample. The first one includes 10 years,
# the second item contains 20 years, and so on until 170 years

load("Box Sync/FinalPaper/DataRepository/PositiveControlSample_Learning1923.Rdata")
#Create a master list
tempdat<-list(tempset=temps_45$tempset1,
              year=temps_45$year,
              datset1=temps_45$S1,
              datset2=temps_45$S2,
              datset3=temps_45$S3,
              datset4=temps_45$S4,
              datset5=temps_45$S5,
              datset6=temps_45$S6,
              datset7=temps_45$S7,
              datset8=temps_45$S8,
              datset9=temps_45$S9,
              datset10=temps_45$S10)
save(tempdat,file="Box Sync/FinalPaper/DataRepository/PositiveControlSample45_1923.Rdata")
load("Box Sync/FinalPaper/DataRepository/PositiveControlSample_Learning1923.Rdata")

# Process the list for the learning study
end<-c(seq(10,170,by=10)*52,length(tempdat$tempset)) # End Years
tempdatlist<-list() # List to hold the learning study data
for(i in 1:length(end)){
  tempdatlist[[i]]<-temps_45[1:end[i],]
}
# Final list file for the learning study
finaltempdat45<-list()
for(i in 1:length(tempdatlist)){
  finaltempdat45[[i]]<-list(tempset=tempdatlist[[i]]$tempset1, # data includes temperature
                            year=tempdatlist[[i]]$year, # Year
                            datset1=tempdatlist[[i]]$S1) # GEV Data
}
# Save File
save(finaltempdat45,file="Box Sync/FinalPaper/DataRepository/FinalPositiveControlSample45_1923.Rdata")

