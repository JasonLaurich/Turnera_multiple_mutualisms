# This file estimates the D Matrix (divergence among population means), and then compares the eigenstructure of D to Gw and the 5 population Gs using various methods.

library(car)
library(psych)
library(MCMCglmm)

setwd("C:/Users/Jason/Desktop/Data/TurnChap1/Dryad Files")

pheno<- read.table("pheno.txt",stringsAsFactors=F,header=F)#phenotypes
names(pheno)<-c('id','site','mat','elave','efn','fn','herk','efn_nec','growth')

data<-pheno[!(pheno$site=='Extra'),]
data<-data[(data$efn>0),]
data<-droplevels(data)

dataG<-data[(!is.na(data$efn_nec)&!is.na(data$elave)&!is.na(data$fn)&!is.na(data$herk)),]

dataG$lnfn<-log(dataG$fn+1)

dataG$site<-as.factor(dataG$site)
str(dataG)

list_df <- split(dataG, dataG$site)

#Now we want to assess D_max the axis of greatest variation among population trait means.
#Generate vectors for population mean phenotypes. 

M_list<- matrix(nrow=17, ncol=4) 
for (i in 1:17){
  data<-list_df[i]
  data<-as.data.frame(data)
  data2<-data[c(4,8,10,7)]
  for (j in 1:4){
    M_list[i,j]<-mean(data2[,j])
  }
}

M_list<-as.data.frame(M_list)
names(M_list)<-c('ela','efn','fn','herk')
MCor<-corr.test(M_list, method="pearson")
var(M_list)
cor(M_list)

MCor

# What if I scale the data first?

IDscale <- function(inputdf,IDcol=1) {
  tmp <- scale(inputdf, center = TRUE, scale = TRUE)#
  tmp[,IDcol]<- inputdf[,IDcol]#restore id
  return(tmp)
}

# Remove site information so we can scale our trait data.
dataJam<-dataG[,-c(2,3)]
jam.data <- as.data.frame(IDscale(dataJam))

# Add it back in 
jam.data$site<-dataG$site

list_df_scaled <- split(jam.data, jam.data$site)

#Generate vectors for population mean (scaled) phenotypes. 

M_list_scaled<- matrix(nrow=17, ncol=4) 
for (i in 1:17){
  data<-list_df_scaled[i]
  data<-as.data.frame(data)
  data2<-data[c(2,6,8,5)]
  for (j in 1:4){
    M_list_scaled[i,j]<-mean(data2[,j])
  }
}

M_list_scaled<-as.data.frame(M_list_scaled)
names(M_list_scaled)<-c('ela','efn','fn','herk')
MCor<-corr.test(M_list, method="pearson")
var(M_list_scaled)
cor(M_list_scaled)

# The var() looks very similar to our observed island-wide G. Let's run a MCMCglmm object for confirmation.

prior<-list(G=NULL, R=list(V=diag(4)*0.5,nu=2))

D_Jam<-MCMCglmm(cbind(ela,efn,fn,herk) ~ 1 , 
                random=NULL, rcov=~us(trait):units, #no trait main effect, should be 0 since data is centered (above)
                family=c("gaussian","gaussian","gaussian","gaussian"),
                pedigree=NULL, data=M_list_scaled, prior=prior,
                #verbose=FALSE,pr=TRUE,nitt=1010, thin=10, burnin=10)
                verbose=FALSE,pr=TRUE,nitt=10100000, thin=1000, burnin=100000)

summary(D_Jam$VCV)

# Save it.
save(D_Jam, file='D_Jam.Rdata')

# Load it.
#save(D_Jam, file='D_Jam.Rdata')

# Ok let's summarize all this
n<-4
traitnames <- c("ela","efn","fn","herk") 
MCMCsamp <- 10000

MCMC_Darray <- array(,c(MCMCsamp,(n^2)))
MCMC_Darray  <- D_Jam$VCV

Darray <- array(,c(n,n,MCMCsamp))
dimnames(Darray) <- list(traitnames,traitnames)

for (j in 1:MCMCsamp){
  #block- 1, animal data 2-26, error 27-51
  D <- matrix(MCMC_Darray[j,1:16],ncol= n)
  Darray[,,j] <- D
}


Dsum <- array(,c(n,n,3))
stats<- c("Mean","L_HPD","U_HPD")
dimnames(Dsum) <- list(traitnames,traitnames,stats)

for (j in 1:4){
  for (k in 1:4){
    Dsum[j,k,1]<-mean(Darray[j,k,])
    Dsum[j,k,2]<-quantile(Darray[j,k,],probs=c(.025))
    Dsum[j,k,3]<-quantile(Darray[j,k,],probs=c(.975))
  }
}

Dsum

round(Dsum, digits=2)

# OK, so we've got our D matrix, now let's compare it to the island metapopulation

load('GM_Jam.Rdata')
summary(GMpr$VCV)

# Create an empty array to store our posterior estimates of G - the + 1 is for the block entry

m <- 2
r <- 2

MCMCarray <- array(,c(MCMCsamp,(n^2)*r+1))
MCMCarray <- GMpr$VCV

Garray <- array(,c(n,n,MCMCsamp))
dimnames(Garray) <- list(traitnames,traitnames)
Parray <- array(,c(n,n,MCMCsamp))
dimnames(Parray) <- list(traitnames,traitnames)

for (j in 1:MCMCsamp){
  #block- 1, animal data 2-26, error 27-51
  G <- matrix(MCMCarray[j,2:17],ncol= n)
  CE <- matrix(MCMCarray[j,1],ncol= n) 
  R <- matrix(MCMCarray[j,18:33],ncol= n)
  Garray[,,j] <- G
  Parray[,,j] <- G + R
}

# So now we have isolated and formatted our estimates of D and G. Assemble them into one array so we can do a comparison

Comp_array<- array(,c(n,n,m,MCMCsamp))
names <- c("G","D") 
dimnames(Comp_array) <- list(traitnames,traitnames,names)

Comp_array[,,1,]<-Garray
Comp_array[,,2,]<-Darray

#OK, let's start with Krzanowski's

#############################################################################################################################################################

# Krzanowski's common subspace

#START
kr.subspace <- function(Gs, vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]] 
  if(length(vec) != m){stop("vec must have length = m")}
  h <- function (g, v){
    AA <- array(, c(n, n, m))  
    for (k in 1:m){
      g.vec <- eigen(g[,,k])$vectors[,1:(v[k])] 
      AA[,,k] <- g.vec %*% t(g.vec)
    }
    H <- apply(AA, 1:2, sum)
    list(H = H, AA = AA)
  }
  #internal function to calculate AA and H
  MCMC.H <- array(, c(n, n, MCMCsamp))
  dimnames(MCMC.H) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[4]])      
  MCMC.AA <- array(, c(n, n, m, MCMCsamp))
  dimnames(MCMC.AA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]], dimnames(Gs)[[4]])
  for (i in 1:MCMCsamp){
    kr <- h(Gs[,,,i], v = vec)
    MCMC.H[,,i] <- kr$H
    MCMC.AA[,,,i] <- kr$AA
  }	
  #calculate AA and H for the ith MCMC sample of the G array		
  avH <- apply(MCMC.H, 1:2, mean)
  rownames(avH) <- dimnames(Gs)[[1]]
  colnames(avH) <- dimnames(Gs)[[1]]
  #calculate the posterior mean H
  avAA <- apply(MCMC.AA, 1:3, mean)
  dimnames(avAA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]])
  #calculate the posterior mean AA
  avH.vec <- eigen(avH)$vectors
  #eigenanalysis of posterior mean H	
  proj<- function(a, b) t(b) %*% a %*% b
  #internal function to do projection
  avH.theta <- matrix(, n, m)
  for (i in 1:n){
    for (i in 1:n){
      avH.theta[i,] <- acos(sqrt(apply(avAA, 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  #angles between the eigenvectors posterior mean H and the posterior mean subspaces of each population
  MCMC.H.val <- matrix(, MCMCsamp, n)
  colnames(MCMC.H.val) <- paste("h", 1:n, sep="")
  for (i in 1:n){
    MCMC.H.val[,i] <- apply(MCMC.H, 3, proj, b = avH.vec[,i])
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean H 
  MCMC.H.theta <- array(, c(n, m, MCMCsamp))
  rownames(MCMC.H.theta) <- paste("h", 1:n, sep="")
  colnames(MCMC.H.theta) <- dimnames(Gs)[[3]]
  for(i in 1:n){
    for(j in 1:MCMCsamp){
      MCMC.H.theta[i,,j] <- acos(sqrt(apply(MCMC.AA[,,,j], 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  #posterior distribution of the angles between the eigenvectors of posterior mean H and the MCMC samples of the subspaces of each population
  list(avAA = avAA, avH = avH, MCMC.AA = MCMC.AA, MCMC.H = MCMC.H, MCMC.H.val = MCMC.H.val, MCMC.H.theta = MCMC.H.theta)
}
#END
m <- dim(Comp_array)[[3]]

MCMCG.kr.DG <- kr.subspace(Comp_array, vec = rep(2,m))

# Save the object so that it can be looked at later.
save(MCMCG.kr.DG, file='MCMCG.kr.DG.Rdata')
#load('MCMCG.kr.DG.Rdata')
MCMCG.kr.DG$avH
eigen(MCMCG.kr.DG$avH)

# Load up the table containing the results of Krzanowski's analysis fitted for the Morrissey-adjusted null matrices. 
kr_ran<-read.table('krzan.txt', na.strings="", header=T, sep='\t')

# Load up the MCMC object containing the results of Krzanowski's analysis on the G matrices we generated by randomly assigning phenotypes
load('MCMCG.kr.randG.Rdata')

HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr.DG$MCMC.H.val)), colMeans(MCMCG.kr.DG$MCMC.H.val))
HPD.H.val

HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr$MCMC.H.val)), colMeans(MCMCG.kr$MCMC.H.val), HPDinterval(as.mcmc(kr_ran)), colMeans(kr_ran), HPDinterval(as.mcmc(MCMCG.kr.randG$MCMC.H.val)), colMeans(MCMCG.kr.randG$MCMC.H.val))
HPD.H.val

vec<-c('h1','h2','h3','h4')
HPD_K<-cbind(HPD.H.val,vec)
HPD_K<-as.data.frame(HPD_K)
names(HPD_K)<-c('min','max','mean','H')
HPD_K$mean<-as.numeric(as.character(HPD_K$mean))
HPD_K$min<-as.numeric(as.character(HPD_K$min))
HPD_K$max<-as.numeric(as.character(HPD_K$max))

#Generate supplemental figure, showing shared subspace between Gw and D
kGD<-ggplot(data=HPD_K,aes(y=mean,x=H))+geom_point(size=3)
kGD<-kGD + geom_errorbar(aes(ymin=min, ymax=max), size=1,width=0.1)
kGD<-kGD +theme(axis.title.y=element_blank(), axis.title.x=element_blank())
kGD<-kGD + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
kGD<-kGD + theme_bw()
kGD<-kGD + scale_x_discrete(labels=c(expression(h[1]),expression(h[2]), expression(h[3]), expression(h[4])))
kGD<-kGD  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
kGD<-kGD  + ylab("Eigenvalue")
kGD<-kGD  + xlab("Eigenvector of H")
kGD<-kGD  + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
kGD<-kGD  + theme(panel.border = element_blank())
kGD<-kGD  + theme(axis.line = element_line(color = 'black'))
kGD<-kGD  + theme(legend.text = element_text(size=12, face="bold"))
kGD<-kGD  + theme(axis.text.x = element_text(size=12, face="bold", colour = "black"))
kGD<-kGD  + theme(axis.text.y = element_text(size=12, face="bold", colour = "black"))
kGD

########################################################################################################################
#Eigenstructure comparison.

Gsum_Jam <- array(,c(n,n,3))
dimnames(Gsum_Jam) <- list(traitnames,traitnames,stats)

for (j in 1:4){
  for (k in 1:4){
    Gsum_Jam[j,k,1]<-mean(Garray[j,k,])
    Gsum_Jam[j,k,2]<-quantile(Garray[j,k,],probs=c(.025))
    Gsum_Jam[j,k,3]<-quantile(Garray[j,k,],probs=c(.975))
  }
}

Jam.eigen<-eigen(Gsum_Jam[,,1]) 
Jam.eigen

D.eigen<-eigen(Dsum[,,1])
D.eigen

gmax<-Jam.eigen$vectors[,1]
gmax

dmax<-D.eigen$vectors[,1]
dmax

acos(t(gmax)%*%dmax)
# 2.79 degrees

#Let's get variation in this theta.

gmax0<-Jam.eigen$vectors[,2]

dmax0<-D.eigen$vectors[,2]

acos(t(gmax)%*%dmax)


#so dmax is our z. Can we calculate e(z) and c(z)

install.packages('matlib')
library(matlib)

e_z<- t(dmax)%*%Gsum_Jam[,,1]%*%dmax
e_z

c_z1<- t(dmax)%*%inv(Gsum_Jam[,,1])%*%dmax
c_z<-c_z1^-1

#Compare D to H?

MCMCG.kr.DG$avH
H.eig<-eigen(MCMCG.kr.DG$avH)

D.eigen
H.eig

Hmax<-H.eig$vectors[,1]
H2<-H.eig$vectors[,2]

d2<-D.eigen$vectors[,2]

acos(t(Hmax)%*%dmax)

acos(t(H2)%*%d2)

#########################################################################################################
# Fourth order?

#START
covtensor <- function(Gs){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  neigten <- n*(n+1)/2 
  #Number of eigentensors
  MCMC.S <- array(,c(neigten, neigten, MCMCsamp))
  dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep=""), paste("e", 1:neigten, sep=""))
  for (k in 1:MCMCsamp){
    MCMCG <- Gs[,,,k] 
    MCMCvarmat <- t(apply(MCMCG, 3, diag)) 
    #find the variances of the kth G and store them 
    MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle)) 
    #find the covariances of the kth G and store them
    MCMC.S[1:n,1:n, k] <- cov(MCMCvarmat, MCMCvarmat) 
    #fill the upper left quadrant of the kth S
    MCMC.S[(n+1):neigten,(n+1):neigten, k] <- 2*cov(MCMCcovmat, MCMCcovmat)
    #fill the lower right quadrant of the kth S
    MCMC.S[1:n,(n+1):neigten, k] <- sqrt(2)*cov(MCMCvarmat, MCMCcovmat)
    #fill the upper right quadrant of the kth S
    MCMC.S[(n+1):neigten,1:n, k] <- sqrt(2)*cov(MCMCcovmat, MCMCvarmat)
    #fill the lower left quadrant of the kthS
  }  
  av.S <- apply(MCMC.S, 1:2, mean)
  #posterior mean S
  av.S.val <- eigen(av.S)$values
  #eigenvalues of posterior mean S 
  av.S.vec <- eigen(av.S)$vectors
  #eigenvalues of posterior mean S
  eTmat <- array(, c(n, n, neigten))
  dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep=""))  
  for (i in 1:neigten){
    emat <- matrix(0, n, n) 
    lowerTriangle(emat) <- 1/sqrt(2)*av.S.vec[(n+1):neigten,i]
    emat <- emat + t(emat)
    diag(emat) <- av.S.vec[1:n,i]
    eTmat[,,i] <- emat 
  }
  #construct the second-order eigentensors of posterior mean S
  eT.eigen <- array(, c(n+1, n, neigten))
  for (i in 1:neigten){
    eT.eigen[1,,i] <- t(eigen(eTmat[,,i])$values) 
    #Eigenvalues of the ith eigentensor
    eT.eigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors 
    #Eigenvectors of the ith eigentensor
    eT.eigen[,,i] <- eT.eigen[,order(abs(eT.eigen[1,,i]), decreasing = T), i]
  }
  MCMC.S.val <- matrix(, MCMCsamp, neigten)
  colnames(MCMC.S.val) <- paste("E", 1:neigten, sep="")
  for (i in 1:MCMCsamp){
    for(j in 1:neigten){
      MCMC.S.val[i,j] <- t(av.S.vec[,j]) %*% MCMC.S[,,i] %*% av.S.vec[,j]
    }
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean S
  av.G.coord <- array(, c(m, neigten, 1))
  dimnames(av.G.coord) <- list(names, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the jth avG for the eigentensors of posterior mean S
  MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
  dimnames(MCMC.G.coord) <- list(names, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    MCMC.G.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
  tensor.summary <- data.frame(rep(av.S.val,each=n), t(data.frame(eT.eigen)))
  colnames(tensor.summary) <- c("S.eigval", "eT.val", traitnames)
  rownames(tensor.summary)<- paste(paste("e", rep(1:neigten, each=n), sep=""), rep(1:n,neigten), sep=".")
  list(tensor.summary = tensor.summary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
}
#END

MCMC.covtensor.DG <- covtensor(Comp_array)

nnonzero <- min(n*(n+1)/2,m-1)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor.DG$MCMC.S.val[,1]), prob=0.95), mean(MCMC.covtensor.DG$MCMC.S.val[,1]))
HPD.eT.2<- cbind(HPDinterval(as.mcmc(MCMC.covtensor.DG$MCMC.S.val[,1:nnonzero]), prob=0.95), mean(MCMC.covtensor.DG$MCMC.S.val[,1]))
round(HPD.eT.val, 3)
round(HPD.eT.2, 3)

#OK, so D and Garray give us our 10,000 replicates. Let's use that to calculate HPD intervals for theta, size and eccentricty.

# Vector to stor theta, size & eccentriticty of D, and for G
stat_array <- array(,c(MCMCsamp,5))
stat_array<-as.data.frame(stat_array)
names(stat_array)<-c('theta','Dsize','Decc','Gsize','Gecc')


for (i in 1:MCMCsamp){
  #First, fit the eigen for each.
  D.eig<- eigen(Darray[,,i])
  D.max<-D.eig$vectors[,1]
  
  G.eig<-eigen(Garray[,,i])
  G.max<-G.eig$vectors[,1]
  
  stat_array[i,1]<- acos(t(G.max)%*%D.max)
  
  stat_array[i,2]<-sum(D.eig$values)
  stat_array[i,3]<-D.eig$values[1]/sum(D.eig$values)
  
  stat_array[i,4]<-sum(G.eig$values)
  stat_array[i,5]<-G.eig$values[1]/sum(G.eig$values)

}

HPDinterval(as.mcmc(stat_array[,1:5]), prob=0.95)
colMeans(stat_array[,1:5])


##################################################################################################
#OK, so now we want to compare the D to the 5 largest Gs.

for (i in 1:5){
  load(paste('GMpop_final',i,".Rdata",sep=""))
  assign(paste("GM",i,sep=""),GM)
}

# specify number of MCMC samples
MCMCsamp <- 10000
# number of traits 
n <- 4
# number of matrices to compare
m <- 5 
# number of random effects specified in the model. In mine animal, block and residuals. But we only have two block levels. 
r <- 2
# trait names
traitnames <- c("ela","efn","fn","herk") 
# matrix labels 
Gnames <- c("BT","Ca","MC","Mu","SA")

# Create an empty array to store our posterior estimates of G - the + 1 is for the block entry
MCMCarray <- array(,c(MCMCsamp,(n^2)*r+1,m))

# Store G matrices as mth element of dim [3]
MCMCarray[,,1] <- GM1$VCV
MCMCarray[,,2] <- GM2$VCV
MCMCarray[,,3] <- GM3$VCV
MCMCarray[,,4] <- GM4$VCV
MCMCarray[,,5] <- GM5$VCV

Garray <- array(,c(n,n,m,MCMCsamp))
dimnames(Garray) <- list(traitnames,traitnames,Gnames)
Parray <- array(,c(n,n,m,MCMCsamp))
dimnames(Parray) <- list(traitnames,traitnames,Gnames)

for (i in 1:m){
  for (j in 1:MCMCsamp){
    #block- 1, animal data 2-26, error 27-51
    G <- matrix(MCMCarray[j,2:17,i],ncol= n)
    CE <- matrix(MCMCarray[j,1,i],ncol= n) 
    R <- matrix(MCMCarray[j,18:33,i],ncol= n)
    Garray[,,i,j] <- G
    Parray[,,i,j] <- G + R
  }
}

# summarize all posterior estimates so we can perform some basic eigenanalysis on our G matrices
Gsum <- array(,c(n,n,3,m))
stats<- c("Mean","L_HPD","U_HPD")
dimnames(Gsum) <- list(traitnames,traitnames,stats,Gnames)

for (i in 1:5){
  for (j in 1:4){
    for (k in 1:4){
      Gsum[j,k,1,i]<-mean(Garray[j,k,i,])
      Gsum[j,k,2,i]<-quantile(Garray[j,k,i,],probs=c(.025))
      Gsum[j,k,3,i]<-quantile(Garray[j,k,i,],probs=c(.975))
    }
  }
}

# That gives us the mean values for genetic variance and covariance for each population, and the 95% HPD intervals.
# Evaluate the eigenstructure of our G matrices. 

BT.eigen<-eigen(Gsum[,,1,1]) 
Ca.eigen<-eigen(Gsum[,,1,2]) 
MC.eigen<-eigen(Gsum[,,1,3]) 
Mu.eigen<-eigen(Gsum[,,1,4]) 
SA.eigen<-eigen(Gsum[,,1,5]) 

#Eigenstructure comparison to each G matrix.


Gsum_Jam <- array(,c(n,n,3))
dimnames(Gsum_Jam) <- list(traitnames,traitnames,stats)

for (j in 1:4){
  for (k in 1:4){
    Gsum_Jam[j,k,1]<-mean(Garray[j,k,])
    Gsum_Jam[j,k,2]<-quantile(Garray[j,k,],probs=c(.025))
    Gsum_Jam[j,k,3]<-quantile(Garray[j,k,],probs=c(.975))
  }
}

Jam.eigen<-eigen(Gsum_Jam[,,1]) 
Jam.eigen

D.eigen<-eigen(Dsum[,,1])
D.eigen

gmax<-Jam.eigen$vectors[,1]
gmax

dmax<-D.eigen$vectors[,1]
dmax

acos(t(gmax)%*%dmax)
# 2.79 degrees

#Let's get variation in this theta.

gmax0<-Jam.eigen$vectors[,2]

dmax0<-D.eigen$vectors[,2]

acos(t(gmax)%*%dmax)


#so dmax is our z. Can we calculate e(z) and c(z)

install.packages('matlib')
library(matlib)

e_z<- t(dmax)%*%Gsum_Jam[,,1]%*%dmax
e_z

c_z1<- t(dmax)%*%inv(Gsum_Jam[,,1])%*%dmax
c_z<-c_z1^-1

#Compare D to H?

MCMCG.kr.DG$avH
H.eig<-eigen(MCMCG.kr.DG$avH)

D.eigen
H.eig

Hmax<-H.eig$vectors[,1]
H2<-H.eig$vectors[,2]

d2<-D.eigen$vectors[,2]

acos(t(Hmax)%*%dmax)

acos(t(H2)%*%d2)

#Let's run Krzanowski's on D v BT

Comp_D_BT<- array(,c(n,n,2,MCMCsamp))
names <- c("BT","D") 
dimnames(Comp_D_BT) <- list(traitnames,traitnames,names)

Comp_D_BT[,,1,]<-Garray[,,1,]
Comp_D_BT[,,2,]<-Darray

m <- dim(Comp_D_BT)[[3]]

MCMCG.kr.D_BT <- kr.subspace(Comp_D_BT, vec = rep(2,m))

cbind(HPDinterval(as.mcmc(MCMCG.kr.D_BT$MCMC.H.val)), colMeans(MCMCG.kr.D_BT$MCMC.H.val))

#lower     upper          
#h1 1.182763874 1.9710667 1.6288964
#h2 0.487358326 1.9597331 1.2159231
#h3 0.122464035 1.4451261 0.7734569
#h4 0.005956595 0.9526952 0.3817236

Comp_D_Ca<- array(,c(n,n,2,MCMCsamp))
names <- c("Ca","D") 
dimnames(Comp_D_Ca) <- list(traitnames,traitnames,names)

Comp_D_Ca[,,1,]<-Garray[,,2,]
Comp_D_Ca[,,2,]<-Darray

MCMCG.kr.D_Ca <- kr.subspace(Comp_D_Ca, vec = rep(2,m))

cbind(HPDinterval(as.mcmc(MCMCG.kr.D_Ca$MCMC.H.val)), colMeans(MCMCG.kr.D_Ca$MCMC.H.val))

#lower     upper          
#h1 1.180287469 1.9880503 1.6363989
#h2 0.709792227 1.9200508 1.3393246
#h3 0.151838019 1.2517604 0.7317247
#h4 0.005513309 0.8495489 0.2925519

Comp_D_MC<- array(,c(n,n,2,MCMCsamp))
names <- c("MC","D") 
dimnames(Comp_D_MC) <- list(traitnames,traitnames,names)

Comp_D_MC[,,1,]<-Garray[,,3,]
Comp_D_MC[,,2,]<-Darray

MCMCG.kr.D_MC <- kr.subspace(Comp_D_MC, vec = rep(2,m))

cbind(HPDinterval(as.mcmc(MCMCG.kr.D_MC$MCMC.H.val)), colMeans(MCMCG.kr.D_MC$MCMC.H.val))

#lower    upper          
#h1 1.450799802 1.997691 1.8149705
#h2 0.126361060 1.680487 0.9546399
#h3 0.076474249 1.651128 0.8827661
#h4 0.008102288 0.908571 0.3476235

Comp_D_Mu<- array(,c(n,n,2,MCMCsamp))
names <- c("Mu","D") 
dimnames(Comp_D_Mu) <- list(traitnames,traitnames,names)

Comp_D_Mu[,,1,]<-Garray[,,4,]
Comp_D_Mu[,,2,]<-Darray

MCMCG.kr.D_Mu <- kr.subspace(Comp_D_Mu, vec = rep(2,m))

cbind(HPDinterval(as.mcmc(MCMCG.kr.D_Mu$MCMC.H.val)), colMeans(MCMCG.kr.D_Mu$MCMC.H.val))

#lower    upper          
#h1 0.93145711 1.953771 1.4871865
#h2 0.65458828 1.898826 1.1982963
#h3 0.07585202 1.533802 0.7788562
#h4 0.03240870 1.049211 0.5356610

Comp_D_SA<- array(,c(n,n,2,MCMCsamp))
names <- c("SA","D") 
dimnames(Comp_D_SA) <- list(traitnames,traitnames,names)

Comp_D_SA[,,1,]<-Garray[,,5,]
Comp_D_SA[,,2,]<-Darray

MCMCG.kr.D_SA <- kr.subspace(Comp_D_SA, vec = rep(2,m))

cbind(HPDinterval(as.mcmc(MCMCG.kr.D_SA$MCMC.H.val)), colMeans(MCMCG.kr.D_SA$MCMC.H.val))

#lower     upper          
#h1 1.3063586433 1.9981314 1.7588526
#h2 0.6130048994 1.9564645 1.2834537
#h3 0.0144113457 1.3237660 0.6253285
#h4 0.0006991012 0.8964551 0.3323652

# So.. MC and Mu look a little different? Still wide HPD bars.

#Calculate the thetas.

BTmax<-BT.eigen$vectors[,1]
Camax<-Ca.eigen$vectors[,1]
MCmax<-MC.eigen$vectors[,1]
Mumax<-Mu.eigen$vectors[,1]
SAmax<-SA.eigen$vectors[,1]

acos(t(BTmax)%*%dmax)
# 0.839

acos(t(Camax)%*%dmax)
# 2.22

acos(t(MCmax)%*%dmax)
# 2.986

acos(t(Mumax)%*%dmax)
# 1.483

acos(t(SAmax)%*%dmax)
# 2.251


#Let's get variation in this theta.

gmax0<-Jam.eigen$vectors[,2]

dmax0<-D.eigen$vectors[,2]

acos(t(gmax)%*%dmax)
