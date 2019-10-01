#######################################################################################
# Model Comparison
# Positive Control 
# AIC, BIC, and DIC 
# Note that AIC and BIC use the MLE estimates
# (1) Stationary
# (2) Mu Only Non-stationary
# (3) Mu Sigma Non-Stationary
# (4) Fully Non-Stationary
#######################################################################################
# Initialize and Load Files
load("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/FinalPositiveControlSample45.Rdata")
#################################################################################################
#################################################################################################
# (1) Stationary
load("PositiveControl_Stat_Output.RData")
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceStat.R")
statpars<-gev.fit(xdat=datset,siglink=exp)$mle
# AIC
stataic<-statgevaic(data=datset,temps=tempset,mu0=statpars[1],sigma_0=statpars[2],xi_0=statpars[3])
stataic 
# BIC
statbic<-statgevbic(data=datset,temps=tempset,mu0=statpars[1],sigma_0=statpars[2],xi_0=statpars[3])
statbic 
#DIC
burnin=6000 # Burnin
statpars<-CredIntervalsGEV(statrun2,burn=burnin)[,2]
stat_posteriordraws<-statrun2[[1]][-(1:burnin),]
devstat_postmean<--2*newloglik(dataset=datset,tempset1=tempset,
                               mu_0=statpars[1],sigma_0=statpars[2],xi_0=statpars[3])
fulldev_stat<-apply(stat_posteriordraws,1,function(x){
  (-2*newloglik(dataset=datset,tempset1=tempset,
                mu_0=x[1],sigma_0=x[2],xi_0=x[3]))
})
avpostdev_stat<-mean(fulldev_stat)
pd_stat<-avpostdev_stat-devstat_postmean
DIC_stat<-avpostdev_stat+pd_stat
DIC_stat 
#################################################################################################
#################################################################################################
# (2) Mu Only Non-stationary
load("PositiveControl_Mu_Output.RData")
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceMu.R")
mu.mle<-gev.fit(xdat=datset,ydat=matrix(tempset,ncol=1),mul=1,siglink=exp)
mupars<-c(mu.mle$mle[c(1,3,4)],(mu.mle$mle[2]/mu.mle$mle[1]))
#AIC
muaic<-gevaic(data=datset,temps=tempset,
              mu0=mupars[1],sigma_0=mupars[2],xi_0=mupars[3],
              a_mu=mupars[4])
muaic
#BIC
mubic<-gevbic(data=datset,temps=tempset,
              mu0=mupars[1],sigma_0=mupars[2],xi_0=mupars[3],
              a_mu=mupars[4])
mubic 
# DIC
burnin=6000 # Burnin
mupars<-CredIntervalsGEV(murun2,burn=burnin)[,2]
mu_posteriordraws<-murun2[[1]][-(1:burnin),]
devmu_postmean<--2*newloglik(dataset=datset,tempset1=tempset,
                             mu_0=mupars[1],sigma_0=mupars[2],xi_0=mupars[3],
                             a_mu=mupars[4])
fulldev_mu<-apply(mu_posteriordraws,1,function(x){
  (-2*newloglik(dataset=datset,tempset1=tempset,
                mu_0=x[1],sigma_0=x[2],xi_0=x[3],
                a_mu=x[4]))})
avpostdev_mu<-mean(fulldev_mu)
pd_mu<-avpostdev_mu-devmu_postmean
DIC_mu<-avpostdev_mu+pd_mu
DIC_mu
#################################################################################################
#################################################################################################
# (3) Mu Sigma Non-Stationary
load("PositiveControl_MuSigma_Output.RData")
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSourceMuSigma.R")
musigma.mle<-gev.fit(xdat=datset,ydat=matrix(tempset,ncol=1),mul=1,sigl=1,siglink=exp,
                     method="BFGS")
musigmapars<-c(musigma.mle$mle[c(1,3,5)],(musigma.mle$mle[2]/musigma.mle$mle[1]),
               (musigma.mle$mle[4]/musigma.mle$mle[3]))

# AIC
musigmaaic<-gevaic(data=datset,temps=tempset,
                   mu0=musigmapars[1],sigma_0=musigmapars[2],xi_0=musigmapars[3],
                   a_mu=musigmapars[4],a_sigma=musigmapars[5])
musigmaaic 
# BIC
musigmabic<-gevbic(data=datset,temps=tempset,
                   mu0=musigmapars[1],sigma_0=musigmapars[2],xi_0=musigmapars[3],
                   a_mu=musigmapars[4],a_sigma=musigmapars[5])
musigmabic 
#DIC
burnin=6000 # Burnin
musigmapars<-CredIntervalsGEV(run2,burn=burnin)[,2]
musigma_posteriordraws<-run2[[1]][-(1:burnin),]
devmusigma_postmean<--2*newloglik(dataset=datset,tempset1=tempset,
                                  mu_0=musigmapars[1],sigma_0=musigmapars[2],xi_0=musigmapars[3],
                                  a_mu=musigmapars[4],a_sigma=musigmapars[5])

fulldev_musigma<-apply(musigma_posteriordraws,1,function(x){
  (-2*newloglik(dataset=datset,tempset1=tempset,
                mu_0=x[1],sigma_0=x[2],xi_0=x[3],
                a_mu=x[4],a_sigma=x[5]))
})
avpostdev_musigma<-mean(fulldev_musigma)

pd_musigma<-avpostdev_musigma-devmusigma_postmean
DIC_musigma<-avpostdev_musigma+pd_musigma
DIC_musigma#20966.91
#################################################################################################
#################################################################################################
# (4) Fully Non-Stationary
load("PositiveControl_FullNS_Output.RData")
source("/Users/su_skl5261/Box Sync/FinalPaper/DataRepository/SourceCode/UniformSource.R")
ns.mle<-gev.fit(xdat=datset,ydat=matrix(tempset,ncol=1),
                mul=1,sigl=1,shl=1,siglink=exp,method="BFGS")
nstatpars<-c(ns.mle$mle[c(1,3,5)],
             (ns.mle$mle[2]/ns.mle$mle[1]),
             (ns.mle$mle[4]/ns.mle$mle[3]),
             (ns.mle$mle[6]/ns.mle$mle[5]))

#AIC
nstataic<-gevaic(data=datset,temps=tempset,
                 mu0=nstatpars[1],sigma_0=nstatpars[2],xi_0=nstatpars[3],
                 a_mu=nstatpars[4],a_sigma=nstatpars[5],a_xi=nstatpars[6])
nstataic
#BIC
nsstatbic<-gevbic(data=datset,temps=tempset, 
                  mu0=nstatpars[1],sigma_0=nstatpars[2],xi_0=nstatpars[3],
                  a_mu=nstatpars[4],a_sigma=nstatpars[5],a_xi=nstatpars[6])
nsstatbic 
#DIC
burnin=6000 # Burnin
nstatpars<-CredIntervalsGEV(run2,burn=burnin)[,2] # Posterior Means
ns_posteriordraws<-run2[[1]][-(1:burnin),]
devns_postmean<--2*newloglik(dataset=datset,tempset1=tempset,
                             mu_0=nstatpars[1],sigma_0=nstatpars[2],xi_0=nstatpars[3],
                             a_mu=nstatpars[4],a_sigma=nstatpars[5],a_xi=nstatpars[6])
fulldev_ns<-apply(ns_posteriordraws,1,function(x){
  (-2*newloglik(dataset=datset,tempset1=tempset,
                mu_0=x[1],sigma_0=x[2],xi_0=x[3],
                a_mu=x[4],a_sigma=x[5],a_xi=x[6]))})
avpostdev_ns<-mean(fulldev_ns)
pd_ns<-avpostdev_ns-devns_postmean
DIC_ns<-avpostdev_ns+pd_ns
DIC_ns
