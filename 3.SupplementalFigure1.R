
############################################################################################################
#
# Create Supplemental Figure 1 
# This is a temperature anomaly (from 1980-2000 average) time series of historical and projected temperature (C)
# Four plots correspond to an RCP scenario. 
#
############################################################################################################
rm(list=ls())
############################################################################################################
# Step 1:
# Load Necessary Files
load("Box Sync/FinalPaper/DataRepository/GEVPostiveControlHistoricalRCP.RData")
load("Box Sync/FinalPaper/DataRepository/rcp_proj.RData")                        
############################################################################################################
# Step 2:
# PLot figure for each RCP scenario
fontsize<-2.8
par(mfrow=c(1,1),mar=c(5, 6.5, 2, 1) + 0.1,mgp=c(4, 1.5, 0))# Plotting Parameters


# RCP 4.5
plot(finalglobalrcp45[1:136,2]~finalglobalrcp45[1:136,1],typ="l",
     xlim=c(min(finalglobalrcp45[,1],na.rm=TRUE),max(finalglobalrcp45[,1],na.rm=TRUE))
     ,ylim=c(min(finalglobalrcp45[,-1],na.rm=TRUE),max(finalglobalrcp45[,-1],na.rm=TRUE))
     ,col="black",lwd=2,
     xlab="Year",ylab=expression("Temp. Anomaly ("~degree~C~")"),
     cex.lab=fontsize,cex.main=2,cex.axis=fontsize)
for(i in 2:ncol(finalglobalrcp45)){
  lines(finalglobalrcp45[136:nrow(finalglobalrcp45),i]~finalglobalrcp45[136:nrow(finalglobalrcp45),1],col="gray")
}
lines(finalglobalrcp[136:nrow(finalglobalrcp),3]~finalglobalrcp[136:nrow(finalglobalrcp),1],col="blue",lwd=2)
abline(v=2015,lwd=2,col="red",lty=2)
legend(x=1980,y=0.5,
       legend = c('Ensemble Mean','CMIP5 Realizations'), 
       bty = "n",
       seg.len=c(0.5,0.5),
       col = c("blue", "gray"),
       lty = c(1,1),
       lwd=2,
       density=c(0,0),
       border = c(NA,NA),
       cex = 2, y.intersp = 0.7,x.intersp = 0.1
)


##################################################################################################################################
# Save PNG with dimensions 1905 vs. 982