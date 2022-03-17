#This code will run a huge number of random G matrices on Stephen's server. Here we will run G's for the 5 populations with higher replication (10000) but store only the relevant output
# Also serves as sample code for fitting null Gw estimates.

library(parallel)
library(foreach)
library(doParallel)


#Use 25 cores. At 1000 replicates, that's 40/core. OK this is a third version, designed to fill out the tables that were interupted due to internet failure.
#25-40, 65-80, etc. 
registerDoParallel(25)
foreach (i=1:25) %dopar%{
  
  library(MCMCglmm)
  
  IDscale <- function(inputdf,IDcol=1) {
    tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
    tmp[,IDcol]<- inputdf[,IDcol]#restore id
    return(tmp)
  }
  
  ###ONLY FOR MY TESTING, not on the server 
  
  #setwd("C:/Users/Jason/Desktop/ExperimentalData/TurnChap1/Final Analysis/Curated files")
  
  prior<-list(G=list(G1=list(V=0.01, nu=2),G2=list(V=diag(4)*0.5,nu=2)), R=list(V=diag(4)*0.5,nu=2)) #prior
  
  ped <- read.table("pedmat.txt",stringsAsFactors=F,header=F) #pedigree
  names(ped)<-c('id','sire','dam','dam.pop','sire.pop')
  
  #To set the starting seed for each core
  x<-i-1
  a<-x*40+1
  z<-a+39
  
  #Create the table that will store the data for these GMS
  
  df_B<-data.frame(matrix(ncol = 10, nrow = 0))
  names(df)<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
  
  df_C<-data.frame(matrix(ncol = 10, nrow = 0))
  names(df)<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
  
  df_MC<-data.frame(matrix(ncol = 10, nrow = 0))
  names(df)<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
  
  df_Mu<-data.frame(matrix(ncol = 10, nrow = 0))
  names(df)<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
  
  df_S<-data.frame(matrix(ncol = 10, nrow = 0))
  names(df)<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
  
  #Now loop through loading datasets, doing transformations and running and saving GMs
  for (i in a:z){
    dataB<-read.table(paste('randatBT',i,'.txt',sep=""), sep = "\t")
    names(dataB)<-c('id','site','mat','elave','efn_nec','lnfn','herk','growth','block')
    
    dataC<-read.table(paste('randatCa',i,'.txt',sep=""), sep = "\t")
    names(dataC)<-c('id','site','mat','elave','efn_nec','lnfn','herk','growth','block')
    
    dataMC<-read.table(paste('randatMC',i,'.txt',sep=""), sep = "\t")
    names(dataMC)<-c('id','site','mat','elave','efn_nec','lnfn','herk','growth','block')
    
    dataMu<-read.table(paste('randatMu',i,'.txt',sep=""), sep = "\t")
    names(dataMu)<-c('id','site','mat','elave','efn_nec','lnfn','herk','growth','block')
    
    dataS<-read.table(paste('randatSA',i,'.txt',sep=""), sep = "\t")
    names(dataS)<-c('id','site','mat','elave','efn_nec','lnfn','herk','growth','block')
    
    dataB<-dataB[,-c(2,3)]
    dataC<-dataC[,-c(2,3)]
    dataMC<-dataMC[,-c(2,3)]
    dataMu<-dataMu[,-c(2,3)]
    dataS<-dataS[,-c(2,3)]
    
    ped$id<-as.integer(ped$id)
    
    #Get our pedigree ligned up with our data points
    masterB<-merge(x=dataB,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    pedB<-masterB[,-c(2:7)]
    pedB<-droplevels(pedB)
    
    #Trim pedigree file
    pedB<-pedB[,-c(4,5)]
    
    #Add number of sires to the pedigree at the top (190 maternal lines)
    glm.pedB <- data.frame(
      id=c(unique(pedB$dam),pedB$id),
      sire = c(rep(NA,length.out=19), pedB$sire), 
      dam=c(rep(NA,length.out=19),pedB$sire))
    
    #Get our pedigree ligned up with our data points
    masterC<-merge(x=dataC,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    pedC<-masterC[,-c(2:7)]
    pedC<-droplevels(pedC)
    
    #Trim pedigree file
    pedC<-pedC[,-c(4,5)]
    
    #Add number of sires to the pedigree at the top (190 maternal lines)
    glm.pedC <- data.frame(
      id=c(unique(pedC$dam),pedC$id),
      sire = c(rep(NA,length.out=20), pedC$sire), 
      dam=c(rep(NA,length.out=20),pedC$sire))
    
    #Get our pedigree ligned up with our data points
    masterMC<-merge(x=dataMC,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    pedMC<-masterMC[,-c(2:7)]
    pedMC<-droplevels(pedMC)
    
    #Trim pedigree file
    pedMC<-pedMC[,-c(4,5)]
    
    #Add number of sires to the pedigree at the top (190 maternal lines)
    glm.pedMC <- data.frame(
      id=c(unique(pedMC$dam),pedMC$id),
      sire = c(rep(NA,length.out=20), pedMC$sire), 
      dam=c(rep(NA,length.out=20),pedMC$sire))
    
    #Get our pedigree ligned up with our data points
    masterMu<-merge(x=dataMu,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    pedMu<-masterMu[,-c(2:7)]
    pedMu<-droplevels(pedMu)
    
    #Trim pedigree file
    pedMu<-pedMu[,-c(4,5)]
    
    #Add number of sires to the pedigree at the top (190 maternal lines)
    glm.pedMu <- data.frame(
      id=c(unique(pedMu$dam),pedMu$id),
      sire = c(rep(NA,length.out=19), pedMu$sire), 
      dam=c(rep(NA,length.out=19),pedMu$sire))
    
    #Get our pedigree ligned up with our data points
    masterS<-merge(x=dataS,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
    pedS<-masterS[,-c(2:7)]
    pedS<-droplevels(pedS)
    
    #Trim pedigree file
    pedS<-pedS[,-c(4,5)]
    
    #Add number of sires to the pedigree at the top (190 maternal lines)
    glm.pedS <- data.frame(
      id=c(unique(pedS$dam),pedS$id),
      sire = c(rep(NA,length.out=19), pedS$sire), 
      dam=c(rep(NA,length.out=19),pedS$sire))
    
    colnames(dataB)[1]<-"animal"
    colnames(glm.pedB)[1]<-"animal"
    
    colnames(dataC)[1]<-"animal"
    colnames(glm.pedC)[1]<-"animal"
    
    colnames(dataMC)[1]<-"animal"
    colnames(glm.pedMC)[1]<-"animal"
    
    colnames(dataMu)[1]<-"animal"
    colnames(glm.pedMu)[1]<-"animal"
    
    colnames(dataS)[1]<-"animal"
    colnames(glm.pedS)[1]<-"animal"
    
    #Need to change block to a character
    dataB$block<-as.numeric(dataB$block)
    dataB$block<-as.factor(dataB$block)
    
    dataC$block<-as.numeric(dataC$block)
    dataC$block<-as.factor(dataC$block)
    
    dataMC$block<-as.numeric(dataMC$block)
    dataMC$block<-as.factor(dataMC$block)
    
    dataMu$block<-as.numeric(dataMu$block)
    dataMu$block<-as.factor(dataMu$block)
    
    dataS$block<-as.numeric(dataS$block)
    dataS$block<-as.factor(dataS$block)
    
    GMran<-MCMCglmm(cbind(elave,efn_nec,lnfn,herk) ~ 1 , 
                    random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                    family=c("gaussian","gaussian","gaussian","gaussian"),
                    pedigree=glm.pedB, data=dataB, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=1010000, thin=100, burnin=10000)
    
    sum<-summary(GMran$VCV)
    val<-as.data.frame(sum$statistics)
    
    var<-val[,1]
    
    vars<-var[2:5]
    vars<-c(vars,var[7:9])
    vars<-c(vars,var[12:13])
    vars<-c(vars,var[17])
    
    df_B[i,]<-vars
    
    write.table(df_B, file=paste("BT_GVCV_est",a,"_ran_final.txt"), sep="\t")
    
    GMran<-MCMCglmm(cbind(elave,efn_nec,lnfn,herk) ~ 1 , 
                    random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                    family=c("gaussian","gaussian","gaussian","gaussian"),
                    pedigree=glm.pedC, data=dataC, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=1010000, thin=100, burnin=10000)
    
    sum<-summary(GMran$VCV)
    val<-as.data.frame(sum$statistics)
    
    var<-val[,1]
    
    vars<-var[2:5]
    vars<-c(vars,var[7:9])
    vars<-c(vars,var[12:13])
    vars<-c(vars,var[17])
    
    df_C[i,]<-vars
    
    write.table(df_C, file=paste("Ca_GVCV_est",a,"_ran_final.txt"), sep="\t")
    
    GMran<-MCMCglmm(cbind(elave,efn_nec,lnfn,herk) ~ 1 , 
                    random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                    family=c("gaussian","gaussian","gaussian","gaussian"),
                    pedigree=glm.pedMC, data=dataMC, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=1010000, thin=100, burnin=10000)
    
    sum<-summary(GMran$VCV)
    val<-as.data.frame(sum$statistics)
    
    var<-val[,1]
    
    vars<-var[2:5]
    vars<-c(vars,var[7:9])
    vars<-c(vars,var[12:13])
    vars<-c(vars,var[17])
    
    df_MC[i,]<-vars
    
    write.table(df_MC, file=paste("MC_GVCV_est",a,"_ran_final.txt"), sep="\t")
    
    GMran<-MCMCglmm(cbind(elave,efn_nec,lnfn,herk) ~ 1 , 
                    random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                    family=c("gaussian","gaussian","gaussian","gaussian"),
                    pedigree=glm.pedMu, data=dataMu, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=1010000, thin=100, burnin=10000)
    
    sum<-summary(GMran$VCV)
    val<-as.data.frame(sum$statistics)
    
    var<-val[,1]
    
    vars<-var[2:5]
    vars<-c(vars,var[7:9])
    vars<-c(vars,var[12:13])
    vars<-c(vars,var[17])
    
    df_Mu[i,]<-vars
    
    write.table(df_Mu, file=paste("Mu_GVCV_est",a,"_ran_final.txt"), sep="\t")
    
    GMran<-MCMCglmm(cbind(elave,efn_nec,lnfn,herk) ~ 1 , 
                    random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                    family=c("gaussian","gaussian","gaussian","gaussian"),
                    pedigree=glm.pedS, data=dataS, prior=prior,
                    verbose=FALSE,pr=TRUE,nitt=1010000, thin=100, burnin=10000)
    
    sum<-summary(GMran$VCV)
    val<-as.data.frame(sum$statistics)
    
    var<-val[,1]
    
    vars<-var[2:5]
    vars<-c(vars,var[7:9])
    vars<-c(vars,var[12:13])
    vars<-c(vars,var[17])
    
    df_S[i,]<-vars
    
    write.table(df_S, file=paste("SA_GVCV_est",a,"_ran_final.txt"), sep="\t")
    
  }
  
  
}

stopImplicitCluster()
