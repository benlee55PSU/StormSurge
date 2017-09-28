###############################################################################################
#
# Use historical and CMIP5 data to create tempearture time series from 1923-2100            
#
#
#################### Read in CMIP5 RCP Data ###########################################################
rm(list=ls())
library(ncdf4) # Set library for CDF files

# List CDF file names for different RCP scenarios
rcp26<-paste("Box Sync/FinalPaper/DataRepository/global_mean_tas/Annual_Mean/RCP26",
             list.files("Box Sync/FinalPaper/DataRepository/global_mean_tas/Annual_Mean/RCP26"),sep="/")

rcp45<-paste("Box Sync/FinalPaper/DataRepository/global_mean_tas/Annual_Mean/RCP45",
             list.files("Box Sync/FinalPaper/DataRepository/global_mean_tas/Annual_Mean/RCP45"),sep="/")

rcp60<-paste("Box Sync/FinalPaper/DataRepository/global_mean_tas/Annual_Mean/RCP60",
             list.files("Box Sync/FinalPaper/DataRepository/global_mean_tas/Annual_Mean/RCP60"),sep="/")

rcp85<-paste("Box Sync/FinalPaper/DataRepository/global_mean_tas/Annual_Mean/RCP85",
             list.files("Box Sync/FinalPaper/DataRepository/global_mean_tas/Annual_Mean/RCP85"),sep="/")


tim<-(1850:2100)# Set Time Span

# Update for RCP2.6
for(i in 1:length(rcp26)){
  if(i==1){
    temp_rcp26 <- ncvar_get(nc_open(rcp26[i]), "tas")[1:251]-273.15 # Convert to Celcius
  }
  else{
    temp_rcp26<-cbind(temp_rcp26,ncvar_get(nc_open(rcp26[i]), "tas")[1:251]-273.15) # Convert to Celcius
  }
}
temp_rcp26_avg<-apply(temp_rcp26,1,mean,na.rm=TRUE) # Take Emsemble Mean

# Update for RCP4.5
for(i in 1:length(rcp45)){
  if(i==1){
    temp_rcp45 <- ncvar_get(nc_open(rcp45[i]), "tas")[1:251]-273.15 # Convert to Celcius
  }
  else{
    temp_rcp45<-cbind(temp_rcp45,ncvar_get(nc_open(rcp45[i]), "tas")[1:251]-273.15) # Convert to Celcius
  }
}
temp_rcp45_avg<-apply(temp_rcp45,1,mean,na.rm=TRUE) # Take Ensemble Mean

# Update for RCP 6.0
for(i in 1:length(rcp60)){
  if(i==1){
    temp_rcp60 <- ncvar_get(nc_open(rcp60[i]), "tas")[1:251]-273.15 # Convert to Celcius
  }
  else{
    temp_rcp60<-cbind(temp_rcp60,ncvar_get(nc_open(rcp60[i]), "tas")[1:251]-273.15) # Convert to Celcius
  }
}
temp_rcp60_avg<-apply(temp_rcp60,1,mean,na.rm=TRUE) # Take Ensemble Mean

# Update for RCP 8.5
for(i in 1:length(rcp85)){
  if(i==1){
    temp_rcp85 <- ncvar_get(nc_open(rcp85[i]), "tas")[1:251]-273.15 # Convert to Celcius
  }
  else{
    temp_rcp85<-cbind(temp_rcp85,ncvar_get(nc_open(rcp85[i]), "tas")[1:251]-273.15) # Convert to Celcius
  }
}
temp_rcp85_avg<-apply(temp_rcp85,1,mean,na.rm=TRUE) # Take Ensemble Mean


#Individual
# RCP2.6
tempdat26<-cbind(tim,temp_rcp26) # Set temp Matrix
baselinetemp_rcp26<-(tempdat26[tempdat26[,1]%in% (2015),]) # Take 2015 temp
tempdat26<-(tempdat26[!tempdat26[,1]%in% (1:2014),]) # Remove records before 2015
# RCP4.5
tempdat45<-cbind(tim,temp_rcp45)
baselinetemp_rcp45<-(tempdat45[tempdat45[,1]%in% (2015),])
tempdat45<-(tempdat45[!tempdat45[,1]%in% (1:2014),])
# RCP6.0
tempdat60<-cbind(tim,temp_rcp60)
baselinetemp_rcp60<-(tempdat60[tempdat60[,1]%in% (2015),])
tempdat60<-(tempdat60[!tempdat60[,1]%in% (1:2014),])
# RCP 8.5
tempdat85<-cbind(tim,temp_rcp85)
baselinetemp_rcp85<-(tempdat85[tempdat85[,1]%in% (2015),])
tempdat85<-(tempdat85[!tempdat85[,1]%in% (1:2014),])

# Take difference between projected and 2015 
for(i in 2:ncol(tempdat26)){ # RCP 2.6
  tempdat26[,i]<-tempdat26[,i]-baselinetemp_rcp26[i]
}
for(i in 2:ncol(tempdat45)){ # RCP 4.5
  tempdat45[,i]<-tempdat45[,i]-baselinetemp_rcp45[i]
}
for(i in 2:ncol(tempdat60)){ # RCP 6.0
  tempdat60[,i]<-tempdat60[,i]-baselinetemp_rcp60[i]
}
for(i in 2:ncol(tempdat85)){ # RCP 8.5
  tempdat85[,i]<-tempdat85[,i]-baselinetemp_rcp85[i]
}

# Plot the ensemble projections
par(mfrow=c(2,2))
plot(tempdat26[,2]~tempdat26[,1],typ="l",ylim=c(min(tempdat26[,2:5]),max(tempdat26[,2:5])))
for(i in 2:ncol(tempdat26)){ # RCP 2.6
  lines(tempdat26[,i]~tempdat26[,1],col="gray") # Problem
}
plot(tempdat45[,2]~tempdat45[,1],typ="l",ylim=c(min(tempdat45[,2:5]),max(tempdat45[,2:5])))
for(i in 2:ncol(tempdat45)){  # RCP 4.5
  lines(tempdat45[,i]~tempdat45[,1],col="gray") # Problem
}
plot(tempdat60[,2]~tempdat60[,1],typ="l",ylim=c(min(tempdat60[,2:5]),max(tempdat60[,2:5])))
for(i in 2:ncol(tempdat60)){ # RCP 6.0
  lines(tempdat60[,i]~tempdat60[,1],col="gray") # Problem
}
plot(tempdat85[,2]~tempdat85[,1],typ="l",ylim=c(min(tempdat85[,2:5]),max(tempdat85[,2:5])))
for(i in 2:ncol(tempdat85)){ # RCP 8.5
  lines(tempdat85[,i]~tempdat85[,1],col="gray") # Problem
}



#Plot Ensemble Average
tempdat<-cbind(tim,temp_rcp26_avg,temp_rcp45_avg,temp_rcp60_avg,temp_rcp85_avg) # Ensemble averages
baselinetemp<-(tempdat[tempdat[,1]%in% (2015),]) # Set 2015 Temperature
tempdat<-(tempdat[!tempdat[,1]%in% (1:2014),]) # Remove records before 2015

for(i in 2:ncol(tempdat)){ # TAke anomalies fro 2015
  tempdat[,i]<-tempdat[,i]-baselinetemp[i]
}
plot(tempdat[,2]~tempdat[,1],typ="l",ylim=c(0,max(tempdat[,2:5]))) # Plot RCP 2.6
lines(tempdat[,3]~tempdat[,1],col="blue") # RCp 4.5
lines(tempdat[,4]~tempdat[,1],col="red")  # RCP 6.0
lines(tempdat[,5]~tempdat[,1],col="green") # RCP 8.5


###############ADD IN HISTORICAL
load("Box Sync/FinalPaper/DataRepository/globaltemp.RData")

newglobtemp26<-matrix(c(globtemp[,1],rep(globtemp[,5],(ncol(tempdat26)-1)))
                      ,nrow=nrow(globtemp),ncol=(ncol(tempdat26)))
newtempdat26<-cbind(tempdat26[,1],tempdat26[,-1]+globtemp[136,5])
finalglobalrcp26<-rbind(newglobtemp26,newtempdat26)

newglobtemp45<-matrix(c(globtemp[,1],rep(globtemp[,5],(ncol(tempdat45)-1)))
                      ,nrow=nrow(globtemp),ncol=(ncol(tempdat45)))
newtempdat45<-cbind(tempdat45[,1],tempdat45[,-1]+globtemp[136,5])
finalglobalrcp45<-rbind(newglobtemp45,newtempdat45)

newglobtemp60<-matrix(c(globtemp[,1],rep(globtemp[,5],(ncol(tempdat60)-1)))
                      ,nrow=nrow(globtemp),ncol=(ncol(tempdat60)))
newtempdat60<-cbind(tempdat60[,1],tempdat60[,-1]+globtemp[136,5])
finalglobalrcp60<-rbind(newglobtemp60,newtempdat60)

newglobtemp85<-matrix(c(globtemp[,1],rep(globtemp[,5],(ncol(tempdat85)-1)))
                      ,nrow=nrow(globtemp),ncol=(ncol(tempdat85)))
newtempdat85<-cbind(tempdat85[,1],tempdat85[,-1]+globtemp[136,5])
finalglobalrcp85<-rbind(newglobtemp85,newtempdat85)
save(finalglobalrcp26,finalglobalrcp45,finalglobalrcp60,finalglobalrcp85,
     file="Box Sync/FinalPaper/DataRepository/rcp_proj.RData")                        

#Ensemble
newglobtemp<-cbind(globtemp[,1],globtemp[,5],globtemp[,5],globtemp[,5],globtemp[,5])
newtempdat<-cbind(tempdat[,1],tempdat[,-1]+globtemp[136,5])
finalglobalrcp<-rbind(newglobtemp,newtempdat)
save(finalglobalrcp,file="Box Sync/FinalPaper/DataRepository/GEVPostiveControlHistoricalRCP.RData")
