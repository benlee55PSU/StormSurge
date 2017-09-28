####################################################################################################
# The following processes historical global average temperature from a txt file
# Downloaded from GHCN-v3 1880-08/2017 + SST: ERSST v5 1880-08/2017
# https://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.txt
################################################################################
rm(list=ls())
dat<-readLines("Box Sync/FinalPaper/DataRepository/globaltemp.txt")
dat<-dat[nchar(dat)!=0]
dat<-dat[-(grep("Year",dat)[-1])]
endw<-gregexpr("[[:blank:]]+",dat[1])[[1]][1:19]-1
startx<-c(0,endw[1:18])
finw<-c(endw-startx,6)

#Read in again as a fixedwidth format
datnew<-read.fwf("Box Sync/FinalPaper/DataRepository/globaltemp.txt",finw)
colnames(datnew)<-as.character(unlist(datnew[1,]))
colnames(datnew)<-gsub("[[:blank:]]","",colnames(datnew))
datnew<-datnew[-c(1,23,24,45,46,67,68,89,90,111,112,133,134),]
datnew<-datnew[-c(length(datnew$Year)),]


for(i in 1:ncol(datnew)){
  datnew[,i]<-as.numeric(gsub("[[:blank:]]","",datnew[,i]))
}

globtemp<-data.frame(year=datnew$Year,glob=datnew$"J-D")
globtemp$tempcel<-globtemp$glob/100
globtemp$abscelsius<-globtemp$tempcel+14 # 14 is the baseline temperature
base1980_2000<-mean(globtemp[globtemp$year%in%(1980:2000),"abscelsius"])
globtemp$relabase<-globtemp$abscelsius-base1980_2000


globtemp$smoothrelabase<-lowess(globtemp$relabase,f = 1/10)$y
plot(x=globtemp$year,y=globtemp$relabase,typ="l")
lines(x=globtemp$year,y=globtemp$smoothrelabase,typ="l")

save(globtemp,file="Box Sync/FinalPaper/DataRepository/globaltemp.RData")
rm(list=ls())

