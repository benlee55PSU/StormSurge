
#################################################################################
# Stochastic Search Variable Selection (SSVS)
#####################################################################################
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/SpikeSource_VarTime_constantp.R")

load("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/FinalPositiveControlSample45.Rdata")
maindat<-finaltempdat45[[14]]
blockmaxdat<-cbind(maindat$year ,maindat$tempset ,maindat$datset1 )
blockmaxdat<-blockmaxdat[blockmaxdat[,1]%in%1923:2010,]
datset=blockmaxdat[,3] # Set Dataset
tempset=blockmaxdat[,2] # Set temperature Set

set.seed(123)
start<-c(rep(1,9)) # Start Value
errvect<-c(0.5, 0.5, 0.5,# Proposal SDs
           0.5, 0.5 ,0.5,
           0.5, 0.5 ,0.5)
uvect<-c(0.01, 0, 0.01,        # Prior Values #1
         0.0001, 0.0001 ,0.0001,
         1, 1 ,1)
sdvect<-c(0.01, 10, 0.01,      # Prior Values #2
          10, 10 ,10,
          NA, NA ,NA)

# Adaptive Phase
burnin=1000
adaptrun1<-optimizemat(int=500,iter=15,dataset=datset,tempset1=tempset,
                       start_sim=start,errvect_sim=errvect,uvect_sim=uvect,sdvect_sim=sdvect) #OK
adaptrun1[[1]];adaptrun1[[2]];run2<-adaptrun1[[3]] # Summaries (Rejection Rate, New Proposal SD)
save(run2,file="C:/Users/Seiyon/OneDrive/FinalPaper/Processing Files/SSVS_Positive_adapt.RData")
load("C:/Users/Seiyon/OneDrive/FinalPaper/Processing Files/SSVS_Positive_adapt.RData")

MULTncrej(run2,burn=burnin) #RejectionRate
run2$errsdmat
MULTnsgevplots(run2,rm.burn=FALSE,burn=burnin) #Trace Plots With Burn-In
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots With Burn-In
CredIntervalsGEV(run2,burn=burnin) 
DensePlot(run2,burn=burnin) #Density Plots

# Better Parameters
adaptrun1[[2]]<-adaptrun1[[2]][!(is.na(adaptrun1[[2]][,1])),]
errvect<-adaptrun1[[2]][nrow(adaptrun1[[2]]),]

# FUll phase
burnin=1000
run1<-MCMC_breaks_resume(50000,5000,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,10000,res=TRUE,resdat=run2,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
MULTnsgevplots(run2,rm.burn=FALSE,burn=burnin) #Trace Plots With Burn-In
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots With Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
CredIntervalsGEV(run2,burn=burnin) 
DensePlot(run2,burn=burnin) #Density Plots
batch<-BM_MCMC(run2,int=500,burn=2000,end=nrow(run2$finmat),bline=2000) #batchmeans
plotonlyBM(batch,burn=burnin) #batchmeans plot only
meansamp<-Mean_MCMC(run2,int=500,burn=2000,end=nrow(run2$finmat),bline=2000) # Sample Means
plotonlyMean(meansamp,burn=burnin) #Sample Means Plot Only

save(run2,file="C:/Users/Seiyon/OneDrive/FinalPaper/SSVS_Positive_final.RData")
load(file="C:/Users/Seiyon/OneDrive/FinalPaper/SSVS_Positive_final.RData")


# Spike and Slab 
# Create a matrix of chosen values
parlist<-run2[[1]][-(1:burnin),1:6]
parlist[,4]<-run2[[1]][-(1:burnin),4]*run2[[1]][-c(1:(burnin-1),nrow(run2$finmat)),7]
parlist[,5]<-run2[[1]][-(1:burnin),5]*run2[[1]][-c(1:(burnin-1),nrow(run2$finmat)),8]
parlist[,6]<-run2[[1]][-(1:burnin),6]*run2[[1]][-c(1:(burnin-1),nrow(run2$finmat)),9]


#Function to put parameter sets into groups
index.par<-function(x){
  if(x[4]==0&x[5]==0&x[6]==0){return(0)}
  else if(x[4]!=0&x[5]==0&x[6]==0){return(1)}
  else if(x[4]==0&x[5]!=0&x[6]==0){return(2)}
  else if(x[4]==0&x[5]==0&x[6]!=0){return(3)}
  else if(x[4]!=0&x[5]!=0&x[6]==0){return(4)}
  else if(x[4]!=0&x[5]==0&x[6]!=0){return(5)}
  else if(x[4]==0&x[5]!=0&x[6]!=0){return(6)}
  else if(x[4]!=0&x[5]!=0&x[6]!=0){return(7)}}

group<-apply(parlist,1,index.par)# Group
#Create Summary Table
summary.tab.spike<-data.frame(
  numpar=c(0,rep(1,3),rep(2,3),3),
  desc=c("None",
         "Location","Scale","Shape",
         "Location and Scale","Location and Shape","Scale and Shape",
         "Location, Scale, and Shape"
  ),
  prop=round(as.numeric(table(group)/nrow(parlist)),4)
)
summary.tab.spike

