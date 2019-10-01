
#################################################################################
# Spike and Slab - All-At-Once Indicators
#NONSTATIONARY

#####################################################################################
source("C:/Users/Seiyon/OneDrive/FinalPaper/SourceCode/SpikeSource_VarTime_constantp.R")
load("C:\\Users\\Seiyon\\OneDrive\\FinalPaper\\Processing Files\\Grinstedmaxblockdata.RData")
blockmaxdat<-newblockmax[,c(1,4,3)]

datset=blockmaxdat[,3] # Set Dataset
tempset=blockmaxdat[,2] # Set temperature Set

start<-c(rep(1,9)) # Start Value
errvect<-c(0.5, 0.5, 0.5,# Proposal SDs
           0.5, 0.5 ,0.5,
           0.5, 0.5 ,0.5)
uvect<-c(0.01, 0, 0.01,        # Prior Values #1
         0.01, 0.01 ,0.01,
         1, 1 ,1)
sdvect<-c(0.01, 10, 0.01,      # Prior Values #2
          10, 10 ,10,
          NA, NA ,NA)

# Adaptive Phase
burnin=1000
adaptrun1<-optimizemat(int=500,iter=15,dataset=datset,tempset1=tempset,
                       start_sim=start,errvect_sim=errvect,uvect_sim=uvect,sdvect_sim=sdvect) #OK
adaptrun1[[1]];adaptrun1[[2]];run2<-adaptrun1[[3]] # Summaries (Rejection Rate, New Proposal SD)
save(run2,file="C:/Users/Seiyon/OneDrive/FinalPaper/Processing Files/SSVS_Grinsted_adapt.RData")
load("C:/Users/Seiyon/OneDrive/FinalPaper/Processing Files/SSVS_Grinsted_adapt.RData")

MULTncrej(run2,burn=burnin) #RejectionRate

MULTnsgevplots(run2,rm.burn=FALSE,burn=burnin) #Trace Plots With Burn-In
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots With Burn-In
CredIntervalsGEV(run2,burn=burnin) 
DensePlot(run2,burn=burnin) #Density Plots

# Better Parameters
adaptrun1[[2]]<-adaptrun1[[2]][!(is.na(adaptrun1[[2]][,1])),]
errvect<-adaptrun1[[2]][nrow(adaptrun1[[2]]),]

# FUll phase
burnin=5000
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

#parameter posterior_mode     Low_95_CI High_95_CI
#1       mu0     2.43800531  2.3496114295  2.5120602
#2    sigma0     0.48027024  0.4300769402  0.5250035
#3       xi0     0.49646962  0.4645720784  0.5398328
#4       amu     0.12228122  0.0000000000  0.2069132
#5    asigma     0.34199690 -0.0007202483  0.5852988
#6       axi     0.05992968 -0.0687834109  0.2852760
save(run2,file="C:/Users/Seiyon/OneDrive/FinalPaper/Processing Files/SSVS_Grinsted_final.RData")
load(file="C:/Users/Seiyon/OneDrive/FinalPaper/Processing Files/SSVS_Grinsted_final.RData")



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



numpar                       desc   prop
1      0                       None 0.8010
2      1                   Location 0.0016
3      1                      Scale 0.0733
4      1                      Shape 0.0140
5      2         Location and Scale 0.1061
6      2         Location and Shape 0.0001
7      2            Scale and Shape 0.0024
8      3 Location, Scale, and Shape 0.0015



