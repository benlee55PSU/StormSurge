
rm(list=ls())
load("C:\\Users\\Seiyon\\OneDrive\\FinalPaper\\Processing Files\\SSVS_Summary_A.RData")
length(summary)

pos_summary<-matrix(NA,ncol=8,nrow=20)
for(i in 1:length(summary)){
  for(j in 0:7){
    pos_summary[i,j+1]<-round(length(summary[[i]][[1]][summary[[i]][[1]]==j])/length(summary[[i]][[1]]),3)  
  }
}


load("C:\\Users\\Seiyon\\OneDrive\\FinalPaper\\Processing Files\\SSVS_Summary_B.RData")
stat_summary<-matrix(NA,ncol=8,nrow=20)
for(i in 1:length(summary)){
  for(j in 0:7){
    stat_summary[i,j+1]<-round(length(summary[[i]][[1]][summary[[i]][[1]]==j])/length(summary[[i]][[1]]),3)  
  }
}




load("C:\\Users\\Seiyon\\OneDrive\\FinalPaper\\Processing Files\\SSVS_Summary_C.RData")
grinsted_summary<-matrix(NA,ncol=8,nrow=20)
for(i in 1:length(summary)){
  for(j in 0:7){
    grinsted_summary[i,j+1]<-round(length(summary[[i]][[1]][summary[[i]][[1]]==j])/length(summary[[i]][[1]]),3)  
  }
}
