##########################################################################################
# This file used to generate the Postive Control Data (Non-Stationary,Stationary, and Grinsted)
# Saves the three data sets into one list 
# Use this in cluster for the SSVS Positive Control Sample

#Set up Skeleton List
ssvs_sample<-list()

#Non-Stationary Data
load("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/FinalPositiveControlSample45.Rdata")
maindat<-finaltempdat45[[14]]
blockmaxdat<-cbind(maindat$year ,maindat$tempset ,maindat$datset1 )
blockmaxdat<-blockmaxdat[blockmaxdat[,1]%in%1923:2010,]
ssvs_sample[[1]]<-blockmaxdat

#Stationary Data
library(fExtremes)
maindat<-finaltempdat45[[14]]
blockmaxdat_stat<-cbind(maindat$year ,maindat$tempset ,maindat$datset1 )
blockmaxdat_stat<-blockmaxdat_stat[blockmaxdat_stat[,1]%in%1923:2010,]
blockmaxdat_stat[,3]<-rgev(length(blockmaxdat_stat[,2]),xi=0.13, mu=2.41, beta=exp(0.48)) # Set Dataset
ssvs_sample[[2]]<-blockmaxdat_stat

#Grinsted Data
load("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/Grinstedmaxblockdata.RData")
ssvs_sample[[3]]<-newblockmax[,c(1,4,3)]


save(ssvs_sample,file="/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SSVS_PositiveControl_Sample.RData")

rm(list=ls())
load("SSVS_PositiveControl_Sample.RData")

