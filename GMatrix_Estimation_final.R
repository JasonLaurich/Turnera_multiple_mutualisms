# This file estimates the G matrices for the Jamaican metapopulation and the 5 individual populations of T. ulmifolia
# It also generates univariate MCMC objects for Bayesian heriability estimates. 

# This will also run the GM for the subset of the Jamaican metapopulation data for which we collected seed set data from pollination trials. 

# Note that much of this code is designed to run in parallel, as these models take a long time to run. This should be considered sample code, and in reality was run on a server. 

# Finally, it contains code for the generation of randomized datasets and running the random GMs, along with sample code for assessing the significance of genetic variance terms in the G matrices.

require(MCMCglmm)

#Will require changing directory
setwd("C:/Users/Jason/Desktop/ExperimentalData/TurnChap1/Final Analysis/Curated files/Dryad files")

# We're going to start making the MCMC objects for the Jamaican metapopulation

# Upload data files and reorganize them 

df<- read.table("pheno_efnpernec.txt",stringsAsFactors=F,header=F)#phenotypes
names(df)<-c('id','site','mat','elave','efn','fn','herk','growth','efn_nec')

ped <- read.table("pedmat.txt",stringsAsFactors=F,header=F)#pedigree
names(ped)<-c('id','sire','dam','dam.pop','sire.pop')

df$mn_elave<-df$elave/mean(df$elave, na.rm=T)
df$mn_efn_nec<-df$efn_nec/mean(df$efn_nec, na.rm=T)
df$mn_fn<-df$fn/mean(df$fn, na.rm=T)
df$mn_herk<-df$herk/mean(df$herk, na.rm=T)

id<-c(1:2000)
block<-c(rep(1,1000),rep(2,1000))
blk<-as.data.frame(cbind(id,block))
df<-merge(x=df,y=blk,by.x="id",by.y="id",all.x=F,all.y=F)

master<-merge(x=df,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
allped<-master[,c(1,15:16)]
allped<-droplevels(allped)

jam.ped <- data.frame(
  id=c(unique(allped$dam),allped$id),
  sire = c(rep(NA,length.out=198), allped$sire), 
  dam=c(rep(NA,length.out=198),allped$sire))


dataJam<-dataG[,-c(2,3)]
jam.data <- as.data.frame(IDscale(dataJam))
id<-c(1:2000)
block<-c(rep(1,1000),rep(2,1000))
blk<-as.data.frame(cbind(id,block))
jam.pheno<-merge(x=jam.data,y=blk,by.x="id",by.y="id",all.x=F,all.y=F)

master<-merge(x=jam.pheno,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
allped<-master[,c(1,10:11)]
allped<-droplevels(allped)

jam.ped <- data.frame(
  id=c(unique(allped$dam),allped$id),
  sire = c(rep(NA,length.out=198), allped$sire), 
  dam=c(rep(NA,length.out=198),allped$sire))

colnames(df)[1]<-"animal"
colnames(jam.ped)[1]<-"animal"

#Need to change block to a character
df$block<-as.character(df$block)

prior<-list(G=list(G1=list(V=0.01, nu=3),G2=list(V=diag(4)*0.5,nu=3)), R=list(V=diag(4)*0.5,nu=3))

GM_fin_lrg<-MCMCglmm(cbind(mn_elave,mn_efn_nec,mn_fn,mn_herk) ~ 1 , 
                 random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                 family=c("gaussian","gaussian","gaussian","gaussian"),
                 pedigree=jam.ped, data=jam.pheno, prior=prior,
                 verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)

save(GM_fin_lrg, file='GM_final_lrg.Rdata')

# Now we run the G matrix (prior4) and univariate MCMC objects (prior1)
# Note that this is where one can insert different priors (ie. to test and select one from many)

# Univariate prior
prior1<-list(G=list(G1=list(V=0.01, nu=2),G2=list(V=diag(1)*0.5,nu=2)), R=list(V=diag(1)*0.5,nu=2))

# 4x4 (G matrix) prior
prior4<-list(G=list(G1=list(V=0.01, nu=2),G2=list(V=diag(4)*0.5,nu=2)), R=list(V=diag(4)*0.5,nu=2))

# Univariate objects
MCela<-MCMCglmm(elave ~ 1 , 
                random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                family="gaussian",
                pedigree=jam.ped, data=jam.pheno, prior=prior1,
                verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
save(MCela, file=paste('Jam_ELA_Uni.Rdata', sep=""))

MCefn<-MCMCglmm(efn_nec ~ 1 , 
                random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                family="gaussian",
                pedigree=jam.ped, data=jam.pheno, prior=prior1,
                verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
save(MCefn, file=paste('Jam_EFN_Uni.Rdata', sep=""))

MCfn<-MCMCglmm(lnfn ~ 1 , 
               random=~block + animal,
               family="gaussian",
               pedigree=jam.ped, data=jam.pheno, prior=prior1,
               verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
save(MCfn, file=paste('Jam_FN_Uni.Rdata', sep=""))

MCherk<-MCMCglmm(herk ~ 1 , 
                 random=~block + animal,
                 family="gaussian",
                 pedigree=jam.ped, data=jam.pheno, prior=prior1,
                 verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
save(MCherk, file=paste('Jam_Hrk_Uni.Rdata', sep=""))

#G Matrix

GM_Jam<-MCMCglmm(cbind(elave,efn_nec,lnfn,herk) ~ 1 , 
                     random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                     family=c("gaussian","gaussian","gaussian","gaussian"),
                     pedigree=jam.ped, data=phenoG, prior=prior4,
                     verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)

save(GM_Jam, file='GM_Jam.Rdata')

#####################################################################################################################################

# Now we will run independent G matrices and univariate MCMC objects for the 5 largest populations of T. ulmifolia we sampled: 
# Browntown, Cave, Mosquito Cove, Murdock, and St. Ann's Bay.

df$site<-as.factor(df$site)
df_pop<-df[(df$site=='BT'|df$site=='Ca'|df$site=='MC'|df$site=='Mu'|df$site=='SA'),]
df_pop<-droplevels(df_pop)
data_pop<-data_pop[,-c(2,3)]

pop.pheno<-merge(x=df_pop,y=blk,by.x="id",by.y="id",all.x=F,all.y=F)

master<-merge(x=pop.pheno,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
ped<-master[,c(1,9:12)]
ped<-droplevels(ped)

ped$dam.pop<-as.numeric(as.factor(ped$dam.pop))
ped$sire.pop<- as.numeric(as.factor(ped$sire.pop))

psi2<-as.data.frame(1:5)
pop.ped <- lapply(1:nrow(psi2), function(z) ped[which(ped$dam.pop==psi2[z,1]),])

pop.pheno <- lapply(1:nrow(psi2), function(z) pop.pheno[which(ped$dam.pop==psi2[z,1]),])

pop.phenoC <- lapply(1:5,function(z) as.data.frame(IDscale(pop.pheno[[z]])))
for (i in 1:5){
  pop.phenoC[[i]]$block<-as.numeric(as.factor(pop.phenoC[[i]]$block))
}

# Function to get the number of unique sires and mothers. Eg. population size, in this breeding design they are equivalent.
popsire=function(z){
  sires<-z[,3]
  uniq<-unique(sires)
  return(length(uniq))
}

glm.pop.ped <- lapply(1:5, function(z) as.matrix(data.frame(
  id=c(unique(pop.ped[[z]]$dam),pop.ped[[z]]$id),
  sire = c(rep(NA,length.out=popsire(pop.ped[[z]])), pop.ped[[z]]$sire), 
  dam=c(rep(NA,length.out=popsire(pop.ped[[z]])),pop.ped[[z]]$dam)))
)

for(i in 1:5){colnames(glm.pop.ped[[i]])[1] <- "animal"}
glm.pop.pheno  <- pop.phenoC
for(i in 1:5){
  colnames(glm.pop.pheno[[i]])[1] <- "animal" }

#The following chunk of code will run a G matrix in addition to each univariate MCMC object for each population on a different core.

registerDoParallel(5)
foreach (i=1:5) %dopar%{
  library(MCMCglmm)
  
  GM<-MCMCglmm(cbind(elave,efn_nec,lnfn,herk) ~ 1 , 
                  random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                  family=c("gaussian","gaussian","gaussian","gaussian"),
                  pedigree=glm.pop.ped[[i]], data=glm.pop.pheno[[i]], prior=prior4,
                  verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
  save(GM, file=paste('GMpop_',i,'.Rdata', sep=""))
  
  MCela<-MCMCglmm(elave ~ 1 , 
                  random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                  family="gaussian",
                  pedigree=glm.pop.ped[[i]], data=glm.pop.pheno[[i]], prior=prior1,
                  verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
  save(MCela, file=paste('MCela_pop',i,'.Rdata', sep=""))
  
  MCefn<-MCMCglmm(efn_nec ~ 1 , 
                  random=~block + animal, #no trait main effect, should be 0 since data is centered (above)
                  family="gaussian",
                  pedigree=glm.pop.ped[[i]], data=glm.pop.pheno[[i]], prior=prior1,
                  verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
  save(MCefn, file=paste('MCefn_pop',i,'.Rdata', sep=""))
  
  MCfn<-MCMCglmm(lnfn ~ 1 , 
                 random=~block + animal,
                 family="gaussian",
                 pedigree=glm.pop.ped[[i]], data=glm.pop.pheno[[i]], prior=prior1,
                 verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
  save(MCfn, file=paste('MCfn_pop',i,'.Rdata', sep=""))
  
  MCherk<-MCMCglmm(herk ~ 1 , 
                   random=~block + animal,
                   family="gaussian",
                   pedigree=glm.pop.ped[[i]], data=glm.pop.pheno[[i]], prior=prior1,
                   verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)
  save(MCherk, file=paste('MCherk_pop',i,'.Rdata', sep=""))
  
}

stopImplicitCluster()

#####################################################################################################################

# Now we will run the G matrix with the seedset data as well 
# This will require uploading the seedset data and re-organizing our files

sd<-read.table('seedset.txt',header=F, na.strings="",sep="\t")
names(sd)<-c('id','mat','pop','trt','seed')

sd$lnseed<-log(sd$seed+1)

# So we want to use estimates of selfed seed set here. 
self<-sd[sd$trt=="self",]

cln.pheno<-merge(x=dataG,y=blk,by.x="id",by.y="id",all.x=F,all.y=F)

self.pheno<-merge(x=cln.pheno,y=self,by.x="id",by.y="id",all.x=F,all.y=F)

phenoSD<-self.pheno[(!is.na(self.pheno$efn_nec)&!is.na(self.pheno$elave)&!is.na(self.pheno$fn)&!is.na(self.pheno$herk)),]
phenoSD<-phenoSD[,c(1,4:11,17)]
phenoSD <- as.data.frame(IDscale(phenoSD))

master<-merge(x=phenoSD,y=ped,by.x="id",by.y="id",all.x=F,all.y=F)
allped<-master[,c(1,11:12)]

#79 unique sires. 

phenoSD$id<-as.integer(phenoG$id)
allped$id<-as.integer(allped$id)

sd.ped <- data.frame(
  id=c(unique(allped$dam),allped$id),
  sire = c(rep(NA,length.out=79), allped$sire), 
  dam=c(rep(NA,length.out=79),allped$sire))

colnames(phenoSD)[1]<-"animal"
colnames(sd.ped)[1]<-"animal"

#Need to change block to a character
phenoSD$block<-as.character(phenoG$block)

# 5x5 (G matrix) prior
prior5<-list(G=list(G1=list(V=0.01, nu=3),G2=list(V=diag(5)*0.5,nu=3)), R=list(V=diag(5)*0.5,nu=3))

GM_seed<-MCMCglmm(cbind(elave,efn_nec,lnfn,herk,lnseed) ~ 1 , 
                     random=~block + us(trait):animal, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                     family=c("gaussian","gaussian","gaussian","gaussian","gaussian"),
                     pedigree=jam.ped, data=phenoG, prior=prior5,
                     verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)

save(GM_seed, file='GM_seed.Rdata')


#############################Generation of randomized data sets################################

#For the island-wide metapopulation
randat<- as.data.frame(dataG[,1:3])

for (i in 1:1000){
    randat_i<-cbind(randat,sample(dataG$elave),sample(dataG$efn_nec),sample(dataG$lnfn),sample(dataG$herk),sample(dataG$block))
    write.table(randat_i, file=paste('randat',i,'.txt', sep=""), sep="\t")
  }

#For the 5 individual, large populations.

dataB<-dataG[dataG$site=='BT',]
dataC<-dataG[dataG$site=='Ca',]
dataMC<-dataG[dataG$site=='MC',]
dataMu<-dataG[dataG$site=='Mu',]
dataS<-dataG[dataG$site=='SA',]

randatB<- as.data.frame(dataB[,1:3])
randatC<- as.data.frame(dataC[,1:3])
randatMC<- as.data.frame(dataMC[,1:3])
randatMu<- as.data.frame(dataMu[,1:3])
randatS<- as.data.frame(dataS[,1:3])

#BT
for (i in 1:1000){
  randat_i<-cbind(randatB, sample(dataB$elave),sample(dataB$lnefn),sample(dataB$efn_nec),sample(dataB$lnfn),sample(dataB$herk),sample(dataB$growth),sample(dataB$block))
  write.table(randat_i, file=paste('randatBT',i,'.txt', sep=""), sep="\t")
}
#Cave
for (i in 1:1000){
  randat_i<-cbind(randatC, sample(dataC$elave),sample(dataC$lnefn),sample(dataC$efn_nec),sample(dataC$lnfn),sample(dataC$herk),sample(dataC$growth),sample(dataC$block))
  write.table(randat_i, file=paste('randatCa',i,'.txt', sep=""), sep="\t")
}
#MC
for (i in 1:1000){
  randat_i<-cbind(randatMC, sample(dataMC$elave),sample(dataMC$lnefn),sample(dataMC$efn_nec),sample(dataMC$lnfn),sample(dataMC$herk),sample(dataMC$growth),sample(dataMC$block))
  write.table(randat_i, file=paste('randatMC',i,'.txt', sep=""), sep="\t")
}
#Mu
for (i in 1:1000){
  randat_i<-cbind(randatMu, sample(dataMu$elave),sample(dataMu$lnefn),sample(dataMu$efn_nec),sample(dataMu$lnfn),sample(dataMu$herk),sample(dataMu$growth),sample(dataMu$block))
  write.table(randat_i, file=paste('randatMu',i,'.txt', sep=""), sep="\t")
}
#SA
for (i in 1:1000){
  randat_i<-cbind(randatS, sample(dataS$elave),sample(dataS$lnefn),sample(dataS$efn_nec),sample(dataS$lnfn),sample(dataS$herk),sample(dataS$growth),sample(dataS$block))
  write.table(randat_i, file=paste('randatSA',i,'.txt', sep=""), sep="\t")
}

#GMs and univariate MCMC objects can then be run on the randomized data sets. 

################################################################################################################

# Assessing the significance of genetic variance estimates (diagonals in the G matrices)
# Bayesian heritability estimates can be found in the heritability code. 

#Browntown
load('GM1.Rdata')
sum<-summary(GM$VCV)
val<-as.data.frame(sum$statistics)

data<-read.table('BT_ran_GVCV.txt', header=T)

table<-data.frame(matrix(NA, nrow=10,ncol=4))
var<-val[,1]
table[1:10,1]<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
table[1:4,2]<-var[2:5]
table[5:7,2]<-var[7:9]
table[8:9,2]<-var[12:13]
table[10,2]<-var[17]

for (i in c(1,5,8,10)){
     quant<-quantile(data[,i],seq(0,1,by=0.025))
     LQ<-quant[2]
     UQ<-quant[40]
     table[i,3]<-LQ
     table[i,4]<-UQ
}

#Cave
load('GM2.Rdata')
sum<-summary(GM$VCV)
val<-as.data.frame(sum$statistics)

data<-read.table('Ca_ran_GVCV.txt', header=T)

table<-data.frame(matrix(NA, nrow=10,ncol=4))
var<-val[,1]
table[1:10,1]<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
table[1:4,2]<-var[2:5]
table[5:7,2]<-var[7:9]
table[8:9,2]<-var[12:13]
table[10,2]<-var[17]

for (i in c(1,5,8,10)){
  quant<-quantile(data[,i],seq(0,1,by=0.025),na.rm =T)
  LQ<-quant[2]
  UQ<-quant[40]
  table[i,3]<-LQ
  table[i,4]<-UQ
}

#Mosquito Cove
load('GM3.Rdata')
sum<-summary(GM$VCV)
val<-as.data.frame(sum$statistics)

data<-read.table('MC_ran_GVCV.txt', header=T)

table<-data.frame(matrix(NA, nrow=10,ncol=4))
var<-val[,1]
table[1:10,1]<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
table[1:4,2]<-var[2:5]
table[5:7,2]<-var[7:9]
table[8:9,2]<-var[12:13]
table[10,2]<-var[17]

for (i in c(1,5,8,10)){
  quant<-quantile(data[,i],seq(0,1,by=0.025),na.rm =T)
  LQ<-quant[2]
  UQ<-quant[40]
  table[i,3]<-LQ
  table[i,4]<-UQ
}

#Murdock
load('GM4.Rdata')
sum<-summary(GM$VCV)
val<-as.data.frame(sum$statistics)

data<-read.table('Mu_ran_GVCV.txt', header=T)

table<-data.frame(matrix(NA, nrow=10,ncol=4))
var<-val[,1]
table[1:10,1]<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
table[1:4,2]<-var[2:5]
table[5:7,2]<-var[7:9]
table[8:9,2]<-var[12:13]
table[10,2]<-var[17]

for (i in c(1,5,8,10)){
  quant<-quantile(data[,i],seq(0,1,by=0.025),na.rm =T)
  LQ<-quant[2]
  UQ<-quant[40]
  table[i,3]<-LQ
  table[i,4]<-UQ
}

#St. Ann's Bay
load('GM5.Rdata')
sum<-summary(GM$VCV)
val<-as.data.frame(sum$statistics)

data<-read.table('SA_ran_GVCV.txt', header=T)

table<-data.frame(matrix(NA, nrow=10,ncol=4))
var<-val[,1]
table[1:10,1]<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
table[1:4,2]<-var[2:5]
table[5:7,2]<-var[7:9]
table[8:9,2]<-var[12:13]
table[10,2]<-var[17]

for (i in c(1,5,8,10)){
  quant<-quantile(data[,i],seq(0,1,by=0.025),na.rm =T)
  LQ<-quant[2]
  UQ<-quant[40]
  table[i,3]<-LQ
  table[i,4]<-UQ
}

#Island Metapopulation
load('Jam_GM.Rdata')
sum<-summary(GM_efn_nec$VCV)
val<-as.data.frame(sum$statistics)

data<-read.table('Jam_ran_GVCV.txt', header=T)

table<-data.frame(matrix(NA, nrow=10,ncol=4))
var<-val[,1]
table[1:10,1]<-c('elV','el-efn','el-fn','el-hrk','efnV','efn-fn','efn-hrk','fnV','fn-hrk','hrkV')
table[1:4,2]<-var[2:5]
table[5:7,2]<-var[7:9]
table[8:9,2]<-var[12:13]
table[10,2]<-var[17]

for (i in c(1,5,8,10)){
  quant<-quantile(data[,i],seq(0,1,by=0.025),na.rm =T)
  LQ<-quant[2]
  UQ<-quant[40]
  table[i,3]<-LQ
  table[i,4]<-UQ
}

################################################################################################################
               
# Fitting average G within populations and D at the same time
                     
GM_Jam_pop_mean<-MCMCglmm(cbind(mn_elave,mn_efn_nec,mn_fn,mn_herk) ~ 1 , 
                          random=~block + us(trait):animal + us(trait):site, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                          family=c("gaussian","gaussian","gaussian","gaussian"),
                          pedigree=jam.ped, data=jam.pheno, prior=prior,
                          verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)

save(GM_Jam_pop_mean, file='JamGwD_final.Rdata')                     
                     
