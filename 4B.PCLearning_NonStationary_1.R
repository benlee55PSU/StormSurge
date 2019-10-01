############################################################################################
#Generate Positive Control Sample for Learning
#Non-Stationary GEV Data
# Generated Using 10 ensemble members
##############################################################################################
# Load Data and restrict to years 1923-to-2100
load("Box Sync/FinalPaper/DataRepository/rcp_proj.RData")

  
ensembleMember<-c(54,  92,  20,  93,  95,  21,  94,  97, 106,  91)
for(i in 1:10){
  finalglobalrcp<-finalglobalrcp45[which(finalglobalrcp45[,1]%in%1923:2100),c(1,ensembleMember[i])]
  print(any(is.na(finalglobalrcp[,2])))
  
}

# Plot the Ensemble members
for(i in 1:length(ensembleMember)){
  finalglobalrcp<-finalglobalrcp45[which(finalglobalrcp45[,1]%in%1923:2100),c(1,ensembleMember[i])]
if(i==1){
  plot(x=finalglobalrcp[,1],y=finalglobalrcp[,2],typ="l",ylim=c(-1,3))
}
  else{lines(x=finalglobalrcp[,1],y=finalglobalrcp[,2],col=i)}
  
}

library(fExtremes)
## Sample Generation for Positive Control
# Initialize
for(k in 1:length(ensembleMember)){
  finalglobalrcp<-finalglobalrcp45[which(finalglobalrcp45[,1]%in%1923:2100),c(1,ensembleMember[k])]


#Generate Positive Control Sample

temps_45<-data.frame(year=rep(finalglobalrcp[,1],each=52),tempset1=rep(finalglobalrcp[,2],each=52)
                     ,mu=NA,sigma=NA,xi=NA) # RCP4.5

# Parameters from Grinsted et al. (2013)
mu0<-2.41 ;sigm0<-0.48 ;xi0<-0.54 ;amu<-0.13;asigma<-0.49; axi<-0.22 # Initialize Grinsted Parameters
trueparmat<-matrix(c(mu0,sigm0,xi0,amu,asigma,axi),nrow=1,ncol=6) # True Parameters
colnames(trueparmat)<-c("mu_0_true","sigma_0_true","xi_0_true","amu_true","asigma_true","axi_true")

# Create Matrices to hold sample data
samplemat_45<-matrix(NA,nrow=nrow(temps_45),30)# RCP4.5
# Generate Data
for(i in 1:nrow(temps_45)){
#set.seed(i*2)
  # RCP 4.5
  temps_45[i,3]<-trueparmat[1,1]*(1+trueparmat[1,4]*temps_45[i,2])#mu
  temps_45[i,4]<-exp(trueparmat[1,2]*(1+trueparmat[1,5]*temps_45[i,2]))#sigma
  temps_45[i,5]<-trueparmat[1,3]*(1+trueparmat[1,6]*temps_45[i,2])#xi
  samplemat_45[i,]<-rgev(30, xi=temps_45[i,5], mu=temps_45[i,3], beta=temps_45[i,4]) # Sample Generation
}

# Create master data frame to hold GEV data and temperatures
temps_45<-cbind(temps_45,samplemat_45) # RCP4.5

# Distinguish each sample by naming each column
for (i in 1:30){
  colnames(temps_45)[i+5]<-paste("S",i,sep="") # Name each column (sample)
}
# Save File
save(temps_45,trueparmat,
     file=paste("Box Sync/FinalPaper/DataRepository/PositiveControlSample_Learning1923_SingleRealization_",k,".Rdata",sep="")) #Save

#######################################################################################
# Save RCP4.5 Data in list form for Learning Study
# Each item in the list is a learning sample. The first one includes 10 years,
# the second item contains 20 years, and so on until 170 years

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
save(tempdat,file=paste("Box Sync/FinalPaper/DataRepository/PositiveControlSample45_1923_SingleRealization_",k,".Rdata",sep=""))

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
save(finaltempdat45,file=paste("Box Sync/FinalPaper/DataRepository/FinalPositiveControlSample45_1923_SingleRealization_",k,".Rdata",sep=""))
}
