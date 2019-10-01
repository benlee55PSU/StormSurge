######################### Pre-Process Grinsted Surge ##############################3
#
# Create the Dataset that will be used in all draws from posterior distributions
# The raw data comes from Aslak Grinsted's website. 
# 
#
rm(list=ls())
setwd("Box Sync/FinalPaper/DataRepository/") # Set working directory
dat<-read.table("SurgeIndex.txt",header = TRUE,sep=",") # Read in data
# Add Year and week
dat$year<-substr(as.character(dat$Day),start=8,stop=12)
dat$week<-rep(1:5531,each=7,length.out =nrow(dat))

# Take weekly block maxima
blockmax<-aggregate(dat[,c(8,9)],by=list(dat$week),FUN=max,na.rm="TRUE")
blockmax<-blockmax[blockmax$SurgeIndex!=-Inf,] # Remove records with no measurements
colnames(blockmax)<-c("week","surge","year") # Add column names
### Temperature Data
load("Box Sync/FinalPaper/DataRepository/globaltemp.RData")
globtemp<-globtemp[,c(1,5,6)] # Add temperature series (actual and smoothed)
newblockmax<-merge(blockmax,globtemp,by="year") #Merge Data
save(newblockmax,file="Grinstedmaxblockdata.RData") # Save Data






