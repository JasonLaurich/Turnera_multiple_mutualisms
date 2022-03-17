# This file contains the code for generating tables partitioning variation in mutualistic phenotypes for all 10,000 posterior estimates of MCMC objects estimate on randomized data sets for each population.
# It will generate summary tables that can then be extracted from a server and used to estimate mean plus 95% HPD intervals of null trait heritability.
# This will be run in parallel on a server due to computational demands.  
library(parallel)
library(foreach)
library(doParallel)


#Use 25 cores. At 1000 replicates, that's 40/core.
#Browntown
registerDoParallel(25)
foreach (i=1:25) %dopar%{
  
  IDscale <- function(inputdf,IDcol=1) {
    tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
    tmp[,IDcol]<- inputdf[,IDcol]#restore id
    return(tmp)
  }
  
  library(MCMCglmm)

  ped <- read.table("pedmat.txt",stringsAsFactors=F,header=F) #pedigree
  names(ped)<-c('id','sire','dam','dam.pop','sire.pop')

# Sample blank vectors to record variance associated with block, breeding values, and units (error) 
  ELAblk<-vector()
  ELAani<-vector()
  ELAuni<-vector()

  EFNblk<-vector()
  EFNani<-vector()
  EFNuni<-vector()

  FNblk<-vector()
  FNani<-vector()
  FNuni<-vector()

  HRKblk<-vector()
  HRKani<-vector()
  HRKuni<-vector()

  GRTblk<-vector()
  GRTani<-vector()
  GRTuni<-vector()
  
  x<-i-1
  a<-x*40+1
  z<-a+39

  for (i in a:z){
    data<-read.table(paste('randatBT',i,'.txt',sep=""), sep = "\t")
    names(data)<-c('id','site','mat','elave','lnefn','lnfn','herk','growth','block')
    
    data<-data[,-c(2,3)]
    
    master<-merge(x=data,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    droplevels(master)
    head(master)
    ped<-master[,-c(2:7)]
    
    cln.data<-as.data.frame(IDscale(data))
    
    
    #Add number of sires to the pedigree at the top (19 maternal lines)
    glm.ped <- data.frame(
      id=c(unique(ped$dam),ped$id),
      sire = c(rep(NA,length.out=19), ped$sire), 
      dam=c(rep(NA,length.out=19),ped$sire))
    
    colnames(cln.data)[1]<-"animal"
    colnames(glm.ped)[1]<-"animal"
    
    
    #Need to change block to a factor
    cln.data$block<-as.numeric(as.factor(cln.data$block))
    
    prior<-list(G=list(G1=list(V=0.01, nu=3),G2=list(V=diag(1)*0.5,nu=3)), R=list(V=diag(1)*0.5,nu=3))
    
    MCela<-MCMCglmm(elave ~ 1 , 
                  random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                  family="gaussian",
                  pedigree=glm.ped, data=cln.data, prior=prior,
                  verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
  
    sum<-as.data.frame(summary(MCela$VCV)$statistics)
    ELAblk<-c(ELAblk,sum[1,1])  
    ELAani<-c(ELAani,sum[2,1])
    ELAuni<-c(ELAuni,sum[3,1])
    
    MCefn<-MCMCglmm(lnefn ~ 1 , 
                  random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                  family="gaussian",
                  pedigree=glm.ped, data=cln.data, prior=prior,
                  verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
  
    sum<-as.data.frame(summary(MCefn$VCV)$statistics)
    EFNblk<-c(EFNblk,sum[1,1])  
    EFNani<-c(EFNani,sum[2,1])
    EFNuni<-c(EFNuni,sum[3,1])
    
    MCfn<-MCMCglmm(lnfn ~ 1 , 
                  random=~block + animal,
                  family="gaussian",
                  pedigree=glm.ped, data=cln.data, prior=prior,
                  verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
  
    sum<-as.data.frame(summary(MCfn$VCV)$statistics)
    FNblk<-c(FNblk,sum[1,1])  
    FNani<-c(FNani,sum[2,1])
    FNuni<-c(FNuni,sum[3,1])
    
    MCherk<-MCMCglmm(herk ~ 1 , 
                  random=~block + animal,
                  family="gaussian",
                  pedigree=glm.ped, data=cln.data, prior=prior,
                  verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCherk$VCV)$statistics)
    HRKblk<-c(HRKblk,sum[1,1])  
    HRKani<-c(HRKani,sum[2,1])
    HRKuni<-c(HRKuni,sum[3,1])
    
    MCgrt<-MCMCglmm(growth ~ 1 , 
                  random=~block + animal,
                  family="gaussian",
                  pedigree=glm.ped, data=cln.data, prior=prior,
                  verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
  
    sum<-as.data.frame(summary(MCgrt$VCV)$statistics)
    GRTblk<-c(GRTblk,sum[1,1])  
    GRTani<-c(GRTani,sum[2,1])
    GRTuni<-c(GRTuni,sum[3,1])
  }

# I have found it preferable to save individual estimates as the code runs, rather than just compiling and saving a single table.
# This will provide a buffer against server failure or disconnection.
    
  ELA<-cbind(ELAblk,ELAani,ELAuni)
  write.table(ELA, file=paste("BT_ELA",a,"her_ran.txt"), sep="\t")

  EFN<-cbind(EFNblk,EFNani,EFNuni)
  write.table(EFN, file=paste("BT_EFN",a,"her_ran.txt"), sep="\t")

  FN<-cbind(FNblk,FNani,FNuni)
  write.table(FN, file=paste("BT_FN",a,"her_ran.txt"), sep="\t")

  HRK<-cbind(HRKblk,HRKani,HRKuni)
  write.table(HRK, file=paste("BT_HRK",a,"her_ran.txt"), sep="\t")

  GRT<-cbind(GRTblk,GRTani,GRTuni)
  write.table(GRT, file=paste("BT_GRT",a,"her_ran.txt"), sep="\t")

}

stopImplicitCluster()

#Assemble the tables
ELA_her<-vector()
EFN_her<-vector()
FN_her<-vector()
HRK_her<-vector()
GRT_her<-vector()

for (i in 1:25){
  a <- (i-1)*40 +1 
  
  ELA<-read.table(paste('BT_ELA',a,'her_ran.txt'),sep="\t",header=T)
  ELA_her<- rbind(ELA_her,ELA)
}

for (j in 1:1000){
  ELA_her$H[j]<-ELA_her[j,2]/sum(ELA_her[j,1:3])
} 
write.table(ELA_her, file="BT_her_ELA.txt", sep="\t")

#Now EFN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  EFN<-read.table(paste('BT_EFN',a,'her_ran.txt'),sep="\t",header=T)
  EFN_her<- rbind(EFN_her,EFN)
}

for (j in 1:1000){
  EFN_her$H[j]<-EFN_her[j,2]/sum(EFN_her[j,1:3])
} 
write.table(EFN_her, file="BT_her_EFN.txt", sep="\t")

#Now FN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  FN<-read.table(paste('BT_FN',a,'her_ran.txt'),sep="\t",header=T)
  FN_her<- rbind(FN_her,FN)
}

for (j in 1:1000){
  FN_her$H[j]<-FN_her[j,2]/sum(FN_her[j,1:3])
} 
write.table(FN_her, file="BT_her_FN.txt", sep="\t")

#Now HRK
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  HRK<-read.table(paste('BT_HRK',a,'her_ran.txt'),sep="\t",header=T)
  HRK_her<- rbind(HRK_her,HRK)
}

for (j in 1:1000){
  HRK_her$H[j]<-HRK_her[j,2]/sum(HRK_her[j,1:3])
} 
write.table(HRK_her, file="BT_her_HRK.txt", sep="\t")

#Now GRT
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  GRT<-read.table(paste('BT_GRT',a,'her_ran.txt'),sep="\t",header=T)
  GRT_her<- rbind(GRT_her,GRT)
}

for (j in 1:1000){
  GRT_her$H[j]<-GRT_her[j,2]/sum(GRT_her[j,1:3])
} 
write.table(GRT_her, file="BT_her_GRT.txt", sep="\t")


#Onto Cave now
#Use 25 cores. At 1000 replicates, that's 40/core.
registerDoParallel(25)
foreach (i=1:25) %dopar%{
  
  IDscale <- function(inputdf,IDcol=1) {
    tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
    tmp[,IDcol]<- inputdf[,IDcol]#restore id
    return(tmp)
  }
  
  library(MCMCglmm)
  
  ped <- read.table("pedmat.txt",stringsAsFactors=F,header=F) #pedigree
  names(ped)<-c('id','sire','dam','dam.pop','sire.pop')
  
  ELAblk<-vector()
  ELAani<-vector()
  ELAuni<-vector()
  
  EFNblk<-vector()
  EFNani<-vector()
  EFNuni<-vector()
  
  FNblk<-vector()
  FNani<-vector()
  FNuni<-vector()
  
  HRKblk<-vector()
  HRKani<-vector()
  HRKuni<-vector()
  
  GRTblk<-vector()
  GRTani<-vector()
  GRTuni<-vector()
  
  x<-i-1
  a<-x*40+1
  z<-a+39
  
  for (i in a:z){
    data<-read.table(paste('randatCa',i,'.txt',sep=""), sep = "\t")
    names(data)<-c('id','site','mat','elave','lnefn','lnfn','herk','growth','block')

    data<-data[,-c(2,3)]
    
    master<-merge(x=data,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    droplevels(master)
    head(master)
    ped<-master[,-c(2:7)]
    
    cln.data<-as.data.frame(IDscale(data))
    
    
    #Add number of sires to the pedigree at the top (19 maternal lines)
    glm.ped <- data.frame(
      id=c(unique(ped$dam),ped$id),
      sire = c(rep(NA,length.out=20), ped$sire), 
      dam=c(rep(NA,length.out=20),ped$sire))
    
    colnames(cln.data)[1]<-"animal"
    colnames(glm.ped)[1]<-"animal"
    
    
    #Need to change block to a factor
    cln.data$block<-as.numeric(as.factor(cln.data$block))
    
    prior<-list(G=list(G1=list(V=0.01, nu=3),G2=list(V=diag(1)*0.5,nu=3)), R=list(V=diag(1)*0.5,nu=3))
    
    MCela<-MCMCglmm(elave ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCela$VCV)$statistics)
    ELAblk<-c(ELAblk,sum[1,1])  
    ELAani<-c(ELAani,sum[2,1])
    ELAuni<-c(ELAuni,sum[3,1])
    
    MCefn<-MCMCglmm(lnefn ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCefn$VCV)$statistics)
    EFNblk<-c(EFNblk,sum[1,1])  
    EFNani<-c(EFNani,sum[2,1])
    EFNuni<-c(EFNuni,sum[3,1])
    
    MCfn<-MCMCglmm(lnfn ~ 1 , 
                   random=~block + animal,
                   family="gaussian",
                   pedigree=glm.ped, data=cln.data, prior=prior,
                   verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCfn$VCV)$statistics)
    FNblk<-c(FNblk,sum[1,1])  
    FNani<-c(FNani,sum[2,1])
    FNuni<-c(FNuni,sum[3,1])
    
    MCherk<-MCMCglmm(herk ~ 1 , 
                     random=~block + animal,
                     family="gaussian",
                     pedigree=glm.ped, data=cln.data, prior=prior,
                     verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCherk$VCV)$statistics)
    HRKblk<-c(HRKblk,sum[1,1])  
    HRKani<-c(HRKani,sum[2,1])
    HRKuni<-c(HRKuni,sum[3,1])
    
    MCgrt<-MCMCglmm(growth ~ 1 , 
                    random=~block + animal,
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCgrt$VCV)$statistics)
    GRTblk<-c(GRTblk,sum[1,1])  
    GRTani<-c(GRTani,sum[2,1])
    GRTuni<-c(GRTuni,sum[3,1])
  }
  
  ELA<-cbind(ELAblk,ELAani,ELAuni)
  write.table(ELA, file=paste("Ca_ELA",a,"her_ran.txt"), sep="\t")
  EFN<-cbind(EFNblk,EFNani,EFNuni)
  write.table(EFN, file=paste("Ca_EFN",a,"her_ran.txt"), sep="\t")
  FN<-cbind(FNblk,FNani,FNuni)
  write.table(FN, file=paste("Ca_FN",a,"her_ran.txt"), sep="\t")
  HRK<-cbind(HRKblk,HRKani,HRKuni)
  write.table(HRK, file=paste("Ca_HRK",a,"her_ran.txt"), sep="\t")
  GRT<-cbind(GRTblk,GRTani,GRTuni)
  write.table(GRT, file=paste("Ca_GRT",a,"her_ran.txt"), sep="\t")
  
}

stopImplicitCluster()

#Assemble the tables
ELA_her<-vector()
EFN_her<-vector()
FN_her<-vector()
HRK_her<-vector()
GRT_her<-vector()

for (i in 1:25){
  a <- (i-1)*40 +1 
  
  ELA<-read.table(paste('Ca_ELA',a,'her_ran.txt'),sep="\t",header=T)
  ELA_her<- rbind(ELA_her,ELA)
}

for (j in 1:1000){
  ELA_her$H[j]<-ELA_her[j,2]/sum(ELA_her[j,1:3])
} 
write.table(ELA_her, file="Ca_her_ELA.txt", sep="\t")

#Now EFN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  EFN<-read.table(paste('Ca_EFN',a,'her_ran.txt'),sep="\t",header=T)
  EFN_her<- rbind(EFN_her,EFN)
}

for (j in 1:1000){
  EFN_her$H[j]<-EFN_her[j,2]/sum(EFN_her[j,1:3])
} 
write.table(EFN_her, file="Ca_her_EFN.txt", sep="\t")

#Now FN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  FN<-read.table(paste('Ca_FN',a,'her_ran.txt'),sep="\t",header=T)
  FN_her<- rbind(FN_her,FN)
}

for (j in 1:1000){
  FN_her$H[j]<-FN_her[j,2]/sum(FN_her[j,1:3])
} 
write.table(FN_her, file="Ca_her_FN.txt", sep="\t")

#Now HRK
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  HRK<-read.table(paste('Ca_HRK',a,'her_ran.txt'),sep="\t",header=T)
  HRK_her<- rbind(HRK_her,HRK)
}

for (j in 1:1000){
  HRK_her$H[j]<-HRK_her[j,2]/sum(HRK_her[j,1:3])
} 
write.table(HRK_her, file="Ca_her_HRK.txt", sep="\t")

#Now GRT
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  GRT<-read.table(paste('Ca_GRT',a,'her_ran.txt'),sep="\t",header=T)
  GRT_her<- rbind(GRT_her,GRT)
}

for (j in 1:1000){
  GRT_her$H[j]<-GRT_her[j,2]/sum(GRT_her[j,1:3])
} 
write.table(GRT_her, file="Ca_her_GRT.txt", sep="\t")

#Onto MosquitoCove now
#Use 25 cores. At 1000 replicates, that's 40/core.
registerDoParallel(25)
foreach (i=1:25) %dopar%{
  
  IDscale <- function(inputdf,IDcol=1) {
    tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
    tmp[,IDcol]<- inputdf[,IDcol]#restore id
    return(tmp)
  }
  
  library(MCMCglmm)
  
  ped <- read.table("pedmat.txt",stringsAsFactors=F,header=F) #pedigree
  names(ped)<-c('id','sire','dam','dam.pop','sire.pop')
  
  ELAblk<-vector()
  ELAani<-vector()
  ELAuni<-vector()
  
  EFNblk<-vector()
  EFNani<-vector()
  EFNuni<-vector()
  
  FNblk<-vector()
  FNani<-vector()
  FNuni<-vector()
  
  HRKblk<-vector()
  HRKani<-vector()
  HRKuni<-vector()
  
  GRTblk<-vector()
  GRTani<-vector()
  GRTuni<-vector()
  
  x<-i-1
  a<-x*40+1
  z<-a+39
  
  for (i in a:z){
    data<-read.table(paste('randatMC',i,'.txt',sep=""), sep = "\t")
    names(data)<-c('id','site','mat','elave','lnefn','lnfn','herk','growth','block')
    
    data<-data[,-c(2,3)]
    
    master<-merge(x=data,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    droplevels(master)
    head(master)
    ped<-master[,-c(2:7)]
    
    cln.data<-as.data.frame(IDscale(data))
    
    
    #Add number of sires to the pedigree at the top (19 maternal lines)
    glm.ped <- data.frame(
      id=c(unique(ped$dam),ped$id),
      sire = c(rep(NA,length.out=20), ped$sire), 
      dam=c(rep(NA,length.out=20),ped$sire))
    
    colnames(cln.data)[1]<-"animal"
    colnames(glm.ped)[1]<-"animal"
    
    
    #Need to change block to a factor
    cln.data$block<-as.numeric(as.factor(cln.data$block))
    
    prior<-list(G=list(G1=list(V=0.01, nu=3),G2=list(V=diag(1)*0.5,nu=3)), R=list(V=diag(1)*0.5,nu=3))
    
    MCela<-MCMCglmm(elave ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCela$VCV)$statistics)
    ELAblk<-c(ELAblk,sum[1,1])  
    ELAani<-c(ELAani,sum[2,1])
    ELAuni<-c(ELAuni,sum[3,1])
    
    MCefn<-MCMCglmm(lnefn ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCefn$VCV)$statistics)
    EFNblk<-c(EFNblk,sum[1,1])  
    EFNani<-c(EFNani,sum[2,1])
    EFNuni<-c(EFNuni,sum[3,1])
    
    MCfn<-MCMCglmm(lnfn ~ 1 , 
                   random=~block + animal,
                   family="gaussian",
                   pedigree=glm.ped, data=cln.data, prior=prior,
                   verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCfn$VCV)$statistics)
    FNblk<-c(FNblk,sum[1,1])  
    FNani<-c(FNani,sum[2,1])
    FNuni<-c(FNuni,sum[3,1])
    
    MCherk<-MCMCglmm(herk ~ 1 , 
                     random=~block + animal,
                     family="gaussian",
                     pedigree=glm.ped, data=cln.data, prior=prior,
                     verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCherk$VCV)$statistics)
    HRKblk<-c(HRKblk,sum[1,1])  
    HRKani<-c(HRKani,sum[2,1])
    HRKuni<-c(HRKuni,sum[3,1])
    
    MCgrt<-MCMCglmm(growth ~ 1 , 
                    random=~block + animal,
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCgrt$VCV)$statistics)
    GRTblk<-c(GRTblk,sum[1,1])  
    GRTani<-c(GRTani,sum[2,1])
    GRTuni<-c(GRTuni,sum[3,1])
  }
  
  ELA<-cbind(ELAblk,ELAani,ELAuni)
  write.table(ELA, file=paste("MC_ELA",a,"her_ran.txt"), sep="\t")
  EFN<-cbind(EFNblk,EFNani,EFNuni)
  write.table(EFN, file=paste("MC_EFN",a,"her_ran.txt"), sep="\t")
  FN<-cbind(FNblk,FNani,FNuni)
  write.table(FN, file=paste("MC_FN",a,"her_ran.txt"), sep="\t")
  HRK<-cbind(HRKblk,HRKani,HRKuni)
  write.table(HRK, file=paste("MC_HRK",a,"her_ran.txt"), sep="\t")
  GRT<-cbind(GRTblk,GRTani,GRTuni)
  write.table(GRT, file=paste("MC_GRT",a,"her_ran.txt"), sep="\t")
  
}

stopImplicitCluster()

#Assemble the tables
ELA_her<-vector()
EFN_her<-vector()
FN_her<-vector()
HRK_her<-vector()
GRT_her<-vector()

for (i in 1:25){
  a <- (i-1)*40 +1 
  
  ELA<-read.table(paste('MC_ELA',a,'her_ran.txt'),sep="\t",header=T)
  ELA_her<- rbind(ELA_her,ELA)
}

for (j in 1:1000){
  ELA_her$H[j]<-ELA_her[j,2]/sum(ELA_her[j,1:3])
} 
write.table(ELA_her, file="MC_her_ELA.txt", sep="\t")

#Now EFN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  EFN<-read.table(paste('MC_EFN',a,'her_ran.txt'),sep="\t",header=T)
  EFN_her<- rbind(EFN_her,EFN)
}

for (j in 1:1000){
  EFN_her$H[j]<-EFN_her[j,2]/sum(EFN_her[j,1:3])
} 
write.table(EFN_her, file="MC_her_EFN.txt", sep="\t")

#Now FN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  FN<-read.table(paste('MC_FN',a,'her_ran.txt'),sep="\t",header=T)
  FN_her<- rbind(FN_her,FN)
}

for (j in 1:1000){
  FN_her$H[j]<-FN_her[j,2]/sum(FN_her[j,1:3])
} 
write.table(FN_her, file="MC_her_FN.txt", sep="\t")

#Now HRK
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  HRK<-read.table(paste('MC_HRK',a,'her_ran.txt'),sep="\t",header=T)
  HRK_her<- rbind(HRK_her,HRK)
}

for (j in 1:1000){
  HRK_her$H[j]<-HRK_her[j,2]/sum(HRK_her[j,1:3])
} 
write.table(HRK_her, file="MC_her_HRK.txt", sep="\t")

#Now GRT
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  GRT<-read.table(paste('MC_GRT',a,'her_ran.txt'),sep="\t",header=T)
  GRT_her<- rbind(GRT_her,GRT)
}

for (j in 1:1000){
  GRT_her$H[j]<-GRT_her[j,2]/sum(GRT_her[j,1:3])
} 
write.table(GRT_her, file="MC_her_GRT.txt", sep="\t")

#Onto Murdock now
#Use 25 cores. At 1000 replicates, that's 40/core.
registerDoParallel(25)
foreach (i=1:25) %dopar%{
  
  IDscale <- function(inputdf,IDcol=1) {
    tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
    tmp[,IDcol]<- inputdf[,IDcol]#restore id
    return(tmp)
  }
  
  library(MCMCglmm)
  
  ped <- read.table("pedmat.txt",stringsAsFactors=F,header=F) #pedigree
  names(ped)<-c('id','sire','dam','dam.pop','sire.pop')
  
  ELAblk<-vector()
  ELAani<-vector()
  ELAuni<-vector()
  
  EFNblk<-vector()
  EFNani<-vector()
  EFNuni<-vector()
  
  FNblk<-vector()
  FNani<-vector()
  FNuni<-vector()
  
  HRKblk<-vector()
  HRKani<-vector()
  HRKuni<-vector()
  
  GRTblk<-vector()
  GRTani<-vector()
  GRTuni<-vector()
  
  x<-i-1
  a<-x*40+1
  z<-a+39
  
  for (i in a:z){
    data<-read.table(paste('randatMu',i,'.txt',sep=""), sep = "\t")
    names(data)<-c('id','site','mat','elave','lnefn','lnfn','herk','growth','block')
    
    data<-data[,-c(2,3)]
    
    master<-merge(x=data,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    droplevels(master)
    head(master)
    ped<-master[,-c(2:7)]
    
    cln.data<-as.data.frame(IDscale(data))
    
    
    #Add number of sires to the pedigree at the top (19 maternal lines)
    glm.ped <- data.frame(
      id=c(unique(ped$dam),ped$id),
      sire = c(rep(NA,length.out=19),ped$sire), 
      dam=c(rep(NA,length.out=19),ped$sire))
    
    colnames(cln.data)[1]<-"animal"
    colnames(glm.ped)[1]<-"animal"
    
    
    #Need to change block to a factor
    cln.data$block<-as.numeric(as.factor(cln.data$block))
    
    prior<-list(G=list(G1=list(V=0.01, nu=3),G2=list(V=diag(1)*0.5,nu=3)), R=list(V=diag(1)*0.5,nu=3))
    
    MCela<-MCMCglmm(elave ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCela$VCV)$statistics)
    ELAblk<-c(ELAblk,sum[1,1])  
    ELAani<-c(ELAani,sum[2,1])
    ELAuni<-c(ELAuni,sum[3,1])
    
    MCefn<-MCMCglmm(lnefn ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCefn$VCV)$statistics)
    EFNblk<-c(EFNblk,sum[1,1])  
    EFNani<-c(EFNani,sum[2,1])
    EFNuni<-c(EFNuni,sum[3,1])
    
    MCfn<-MCMCglmm(lnfn ~ 1 , 
                   random=~block + animal,
                   family="gaussian",
                   pedigree=glm.ped, data=cln.data, prior=prior,
                   verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCfn$VCV)$statistics)
    FNblk<-c(FNblk,sum[1,1])  
    FNani<-c(FNani,sum[2,1])
    FNuni<-c(FNuni,sum[3,1])
    
    MCherk<-MCMCglmm(herk ~ 1 , 
                     random=~block + animal,
                     family="gaussian",
                     pedigree=glm.ped, data=cln.data, prior=prior,
                     verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCherk$VCV)$statistics)
    HRKblk<-c(HRKblk,sum[1,1])  
    HRKani<-c(HRKani,sum[2,1])
    HRKuni<-c(HRKuni,sum[3,1])
    
    MCgrt<-MCMCglmm(growth ~ 1 , 
                    random=~block + animal,
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCgrt$VCV)$statistics)
    GRTblk<-c(GRTblk,sum[1,1])  
    GRTani<-c(GRTani,sum[2,1])
    GRTuni<-c(GRTuni,sum[3,1])
  }
  
  ELA<-cbind(ELAblk,ELAani,ELAuni)
  write.table(ELA, file=paste("Mu_ELA",a,"her_ran.txt"), sep="\t")
  EFN<-cbind(EFNblk,EFNani,EFNuni)
  write.table(EFN, file=paste("Mu_EFN",a,"her_ran.txt"), sep="\t")
  FN<-cbind(FNblk,FNani,FNuni)
  write.table(FN, file=paste("Mu_FN",a,"her_ran.txt"), sep="\t")
  HRK<-cbind(HRKblk,HRKani,HRKuni)
  write.table(HRK, file=paste("Mu_HRK",a,"her_ran.txt"), sep="\t")
  GRT<-cbind(GRTblk,GRTani,GRTuni)
  write.table(GRT, file=paste("Mu_GRT",a,"her_ran.txt"), sep="\t")
  
}

stopImplicitCluster()

#Assemble the tables
ELA_her<-vector()
EFN_her<-vector()
FN_her<-vector()
HRK_her<-vector()
GRT_her<-vector()

for (i in 1:25){
  a <- (i-1)*40 +1 
  
  ELA<-read.table(paste('Mu_ELA',a,'her_ran.txt'),sep="\t",header=T)
  ELA_her<- rbind(ELA_her,ELA)
}

for (j in 1:1000){
  ELA_her$H[j]<-ELA_her[j,2]/sum(ELA_her[j,1:3])
} 
write.table(ELA_her, file="Mu_her_ELA.txt", sep="\t")

#Now EFN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  EFN<-read.table(paste('Mu_EFN',a,'her_ran.txt'),sep="\t",header=T)
  EFN_her<- rbind(EFN_her,EFN)
}

for (j in 1:1000){
  EFN_her$H[j]<-EFN_her[j,2]/sum(EFN_her[j,1:3])
} 
write.table(EFN_her, file="Mu_her_EFN.txt", sep="\t")

#Now FN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  FN<-read.table(paste('Mu_FN',a,'her_ran.txt'),sep="\t",header=T)
  FN_her<- rbind(FN_her,FN)
}

for (j in 1:1000){
  FN_her$H[j]<-FN_her[j,2]/sum(FN_her[j,1:3])
} 
write.table(FN_her, file="Mu_her_FN.txt", sep="\t")

#Now HRK
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  HRK<-read.table(paste('Mu_HRK',a,'her_ran.txt'),sep="\t",header=T)
  HRK_her<- rbind(HRK_her,HRK)
}

for (j in 1:1000){
  HRK_her$H[j]<-HRK_her[j,2]/sum(HRK_her[j,1:3])
} 
write.table(HRK_her, file="Mu_her_HRK.txt", sep="\t")

#Now GRT
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  GRT<-read.table(paste('Mu_GRT',a,'her_ran.txt'),sep="\t",header=T)
  GRT_her<- rbind(GRT_her,GRT)
}

for (j in 1:1000){
  GRT_her$H[j]<-GRT_her[j,2]/sum(GRT_her[j,1:3])
} 
write.table(GRT_her, file="Mu_her_GRT.txt", sep="\t")

#Onto St Ann now
#Use 25 cores. At 1000 replicates, that's 40/core.
registerDoParallel(25)
foreach (i=1:25) %dopar%{
  
  IDscale <- function(inputdf,IDcol=1) {
    tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
    tmp[,IDcol]<- inputdf[,IDcol]#restore id
    return(tmp)
  }
  
  library(MCMCglmm)
  
  ped <- read.table("pedmat.txt",stringsAsFactors=F,header=F) #pedigree
  names(ped)<-c('id','sire','dam','dam.pop','sire.pop')
  
  ELAblk<-vector()
  ELAani<-vector()
  ELAuni<-vector()
  
  EFNblk<-vector()
  EFNani<-vector()
  EFNuni<-vector()
  
  FNblk<-vector()
  FNani<-vector()
  FNuni<-vector()
  
  HRKblk<-vector()
  HRKani<-vector()
  HRKuni<-vector()
  
  GRTblk<-vector()
  GRTani<-vector()
  GRTuni<-vector()
  
  x<-i-1
  a<-x*40+1
  z<-a+39
  
  for (i in a:z){
    data<-read.table(paste('randatSA',i,'.txt',sep=""), sep = "\t")
    names(data)<-c('id','site','mat','elave','lnefn','lnfn','herk','growth','block')
    
    data<-data[,-c(2,3)]
    
    master<-merge(x=data,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    droplevels(master)
    head(master)
    ped<-master[,-c(2:7)]
    
    cln.data<-as.data.frame(IDscale(data))
    
    
    #Add number of sires to the pedigree at the top (19 maternal lines)
    glm.ped <- data.frame(
      id=c(unique(ped$dam),ped$id),
      sire = c(rep(NA,length.out=19), ped$sire), 
      dam=c(rep(NA,length.out=19),ped$sire))
    
    colnames(cln.data)[1]<-"animal"
    colnames(glm.ped)[1]<-"animal"
    
    
    #Need to change block to a factor
    cln.data$block<-as.numeric(as.factor(cln.data$block))
    
    prior<-list(G=list(G1=list(V=0.01, nu=3),G2=list(V=diag(1)*0.5,nu=3)), R=list(V=diag(1)*0.5,nu=3))
    
    MCela<-MCMCglmm(elave ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCela$VCV)$statistics)
    ELAblk<-c(ELAblk,sum[1,1])  
    ELAani<-c(ELAani,sum[2,1])
    ELAuni<-c(ELAuni,sum[3,1])
    
    MCefn<-MCMCglmm(lnefn ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCefn$VCV)$statistics)
    EFNblk<-c(EFNblk,sum[1,1])  
    EFNani<-c(EFNani,sum[2,1])
    EFNuni<-c(EFNuni,sum[3,1])
    
    MCfn<-MCMCglmm(lnfn ~ 1 , 
                   random=~block + animal,
                   family="gaussian",
                   pedigree=glm.ped, data=cln.data, prior=prior,
                   verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCfn$VCV)$statistics)
    FNblk<-c(FNblk,sum[1,1])  
    FNani<-c(FNani,sum[2,1])
    FNuni<-c(FNuni,sum[3,1])
    
    MCherk<-MCMCglmm(herk ~ 1 , 
                     random=~block + animal,
                     family="gaussian",
                     pedigree=glm.ped, data=cln.data, prior=prior,
                     verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCherk$VCV)$statistics)
    HRKblk<-c(HRKblk,sum[1,1])  
    HRKani<-c(HRKani,sum[2,1])
    HRKuni<-c(HRKuni,sum[3,1])
    
    MCgrt<-MCMCglmm(growth ~ 1 , 
                    random=~block + animal,
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCgrt$VCV)$statistics)
    GRTblk<-c(GRTblk,sum[1,1])  
    GRTani<-c(GRTani,sum[2,1])
    GRTuni<-c(GRTuni,sum[3,1])
  }
  
  ELA<-cbind(ELAblk,ELAani,ELAuni)
  write.table(ELA, file=paste("SA_ELA",a,"her_ran.txt"), sep="\t")
  EFN<-cbind(EFNblk,EFNani,EFNuni)
  write.table(EFN, file=paste("SA_EFN",a,"her_ran.txt"), sep="\t")
  FN<-cbind(FNblk,FNani,FNuni)
  write.table(FN, file=paste("SA_FN",a,"her_ran.txt"), sep="\t")
  HRK<-cbind(HRKblk,HRKani,HRKuni)
  write.table(HRK, file=paste("SA_HRK",a,"her_ran.txt"), sep="\t")
  GRT<-cbind(GRTblk,GRTani,GRTuni)
  write.table(GRT, file=paste("SA_GRT",a,"her_ran.txt"), sep="\t")
  
}

stopImplicitCluster()

#Assemble the tables
ELA_her<-vector()
EFN_her<-vector()
FN_her<-vector()
HRK_her<-vector()
GRT_her<-vector()

for (i in 1:25){
  a <- (i-1)*40 +1 
  
  ELA<-read.table(paste('SA_ELA',a,'her_ran.txt'),sep="\t",header=T)
  ELA_her<- rbind(ELA_her,ELA)
}

for (j in 1:1000){
  ELA_her$H[j]<-ELA_her[j,2]/sum(ELA_her[j,1:3])
} 
write.table(ELA_her, file="SA_her_ELA.txt", sep="\t")

#Now EFN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  EFN<-read.table(paste('SA_EFN',a,'her_ran.txt'),sep="\t",header=T)
  EFN_her<- rbind(EFN_her,EFN)
}

for (j in 1:1000){
  EFN_her$H[j]<-EFN_her[j,2]/sum(EFN_her[j,1:3])
} 
write.table(EFN_her, file="SA_her_EFN.txt", sep="\t")

#Now FN
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  FN<-read.table(paste('SA_FN',a,'her_ran.txt'),sep="\t",header=T)
  FN_her<- rbind(FN_her,FN)
}

for (j in 1:1000){
  FN_her$H[j]<-FN_her[j,2]/sum(FN_her[j,1:3])
} 
write.table(FN_her, file="SA_her_FN.txt", sep="\t")

#Now HRK
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  HRK<-read.table(paste('SA_HRK',a,'her_ran.txt'),sep="\t",header=T)
  HRK_her<- rbind(HRK_her,HRK)
}

for (j in 1:1000){
  HRK_her$H[j]<-HRK_her[j,2]/sum(HRK_her[j,1:3])
} 
write.table(HRK_her, file="SA_her_HRK.txt", sep="\t")

#Now GRT
for (i in 1:25){
  a <- (i-1)*40 +1 
  
  GRT<-read.table(paste('SA_GRT',a,'her_ran.txt'),sep="\t",header=T)
  GRT_her<- rbind(GRT_her,GRT)
}

for (j in 1:1000){
  GRT_her$H[j]<-GRT_her[j,2]/sum(GRT_her[j,1:3])
} 
write.table(GRT_her, file="SA_her_GRT.txt", sep="\t")

####################################################################################################################

# Now we will move on to the island metapopulation.

#Use 25 cores. At 1000 replicates, that's 40/core.
registerDoParallel(25)
foreach (i=1:25) %dopar%{
  
  IDscale <- function(inputdf,IDcol=1) {
    tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
    tmp[,IDcol]<- inputdf[,IDcol]#restore id
    return(tmp)
  }
  
  library(MCMCglmm)
  
  IDscale <- function(inputdf,IDcol=1) {
    tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
    tmp[,IDcol]<- inputdf[,IDcol]#restore id
    return(tmp)
  }
  
  ped <- read.table("pedmat.txt",stringsAsFactors=F,header=F) #pedigree
  names(ped)<-c('id','sire','dam','dam.pop','sire.pop')
  
  ELAblk<-vector()
  ELAani<-vector()
  ELAuni<-vector()
  
  EFNblk<-vector()
  EFNani<-vector()
  EFNuni<-vector()
  
  FNblk<-vector()
  FNani<-vector()
  FNuni<-vector()
  
  HRKblk<-vector()
  HRKani<-vector()
  HRKuni<-vector()
  
  GRTblk<-vector()
  GRTani<-vector()
  GRTuni<-vector()
  
  x<-i-1
  a<-x*40+1
  z<-a+39
  
  for (i in a:z){
    data<-read.table(paste('randat',i,'.txt',sep=""), sep = "\t")
    names(data)<-c('id','site','mat','elave','lnefn','lnfn','herk','growth','block')
    
    dataG<-data[,-c(2,3)]
    cln.data <- as.data.frame(IDscale(dataG))
    
    master<-merge(x=cln.data,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    ped<-master[,-c(2:7)]
    ped<-droplevels(ped)
    
    cln.data$id<-as.integer(cln.data$id)
    ped$id<-as.integer(ped$id)
    
    allped<-ped[,-c(4,5)]
    
    glm.ped <- data.frame(
      id=c(unique(allped$dam),allped$id),
      sire = c(rep(NA,length.out=190), allped$sire), 
      dam=c(rep(NA,length.out=190),allped$sire))
    
    colnames(cln.data)[1]<-"animal"
    colnames(glm.ped)[1]<-"animal"
    
    #Need to change block to a character
    cln.data$block<-as.character(cln.data$block)
    
    prior<-list(G=list(G1=list(V=0.01, nu=3),G2=list(V=diag(1)*0.5,nu=3)), R=list(V=diag(1)*0.5,nu=3))
    
    MCela<-MCMCglmm(elave ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCela$VCV)$statistics)
    ELAblk<-c(ELAblk,sum[1,1])  
    ELAani<-c(ELAani,sum[2,1])
    ELAuni<-c(ELAuni,sum[3,1])
    
    MCefn<-MCMCglmm(lnefn ~ 1 , 
                    random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCefn$VCV)$statistics)
    EFNblk<-c(EFNblk,sum[1,1])  
    EFNani<-c(EFNani,sum[2,1])
    EFNuni<-c(EFNuni,sum[3,1])
    
    MCfn<-MCMCglmm(lnfn ~ 1 , 
                   random=~block + animal,
                   family="gaussian",
                   pedigree=glm.ped, data=cln.data, prior=prior,
                   verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCfn$VCV)$statistics)
    FNblk<-c(FNblk,sum[1,1])  
    FNani<-c(FNani,sum[2,1])
    FNuni<-c(FNuni,sum[3,1])
    
    MCherk<-MCMCglmm(herk ~ 1 , 
                     random=~block + animal,
                     family="gaussian",
                     pedigree=glm.ped, data=cln.data, prior=prior,
                     verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCherk$VCV)$statistics)
    HRKblk<-c(HRKblk,sum[1,1])  
    HRKani<-c(HRKani,sum[2,1])
    HRKuni<-c(HRKuni,sum[3,1])
    
    MCgrt<-MCMCglmm(growth ~ 1 , 
                    random=~block + animal,
                    family="gaussian",
                    pedigree=glm.ped, data=cln.data, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=10000)
    
    sum<-as.data.frame(summary(MCgrt$VCV)$statistics)
    GRTblk<-c(GRTblk,sum[1,1])  
    GRTani<-c(GRTani,sum[2,1])
    GRTuni<-c(GRTuni,sum[3,1])
  }
  
  ELA<-cbind(ELAblk,ELAani,ELAuni)
  write.table(ELA, file=paste("Jam_ELA",a,"her_ran.txt"), sep="\t")
  
  EFN<-cbind(EFNblk,EFNani,EFNuni)
  write.table(EFN, file=paste("Jam_EFN",a,"her_ran.txt"), sep="\t")
  
  FN<-cbind(FNblk,FNani,FNuni)
  write.table(FN, file=paste("Jam_FN",a,"her_ran.txt"), sep="\t")
  
  HRK<-cbind(HRKblk,HRKani,HRKuni)
  write.table(HRK, file=paste("Jam_HRK",a,"her_ran.txt"), sep="\t")
  
  GRT<-cbind(GRTblk,GRTani,GRTuni)
  write.table(GRT, file=paste("Jam_GRT",a,"her_ran.txt"), sep="\t")
  
}

stopImplicitCluster()