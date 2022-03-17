# This R script uploads G matrices from the five largest populations of Turnera ulmifolia we sampled and compares them using a variety of methods.
# Most of the code is sourced from Aguirre et. al. 2014, with modifications of the null matrix fitting as per Morrissey et al. 2019.

# This code will run all 4 comparison techniques we employed. (1) Random Skewers, (2) Krzanowski's common subspace, (3) 4th order covariance tensor and (4) Flury's hierarchical analysis.

library(MCMCglmm)
library(gdata)
library(matrixcalc)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(parallel)
library(foreach)
library(doParallel)
library(Rmisc)

#Set directory
setwd("C:/Users/Jason/Desktop/Data/TurnChap1/Dryad Files")

#Load 5 G matrices
for (i in 1:5){
  load(paste('GMpop_final',i,".Rdata",sep=""))
  assign(paste("GM",i,sep=""),GM)
}

#Examine the estimates of genetic variance and correlation.
summary(GM1$VCV)
summary(GM2$VCV)
summary(GM3$VCV)
summary(GM4$VCV)
summary(GM5$VCV)

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

#Now, let's load the randomly generated G matrices for the purposes of significance testing
BT<-read.table('BT_ranG_final.txt', header= T)
Ca<-read.table('Ca_ranG_final.txt', header= T)
MC<-read.table('MC_ranG_final.txt', header= T)
Mu<-read.table('Mu_ranG_final.txt', header= T)
SA<-read.table('SA_ranG_final.txt', header= T)

BT_ran_sum<-sapply(BT, CI)
Ca_ran_sum<-sapply(Ca, CI)
MC_ran_sum<-sapply(MC, CI)
Mu_ran_sum<-sapply(Mu, CI)
SA_ran_sum<-sapply(SA, CI)

ranG_sum<-as.data.frame(rbind(BT_ran_sum, Ca_ran_sum, MC_ran_sum, Mu_ran_sum, SA_ran_sum))
names(ranG_sum)<-c('ela','ela-efn','ela-fn','ela-hrk','efn','efn-fn','efn-hrk','fn','fn-hrk','hrk')
ranG_sum$pop<-c(rep('BT',3),rep('Ca',3),rep('MC',3),rep('Mu',3),rep('SA',3))

write.table(ranG_sum, file="ranGs_summary.txt", sep="\t")

ranG_gvar<-ranG_sum[,c(1,5,8,10,11)]
write.table(ranG_gvar, file="ranGs_genetic_variance.txt", sep="\t")

##########################################################################################################################################################

# The random skewer method (Aguirre 2014)
# Assesses the amount of genetic variance in each matrix in the direction of each random skewer.
# Skewers that are different assembled into R matrix, eigenvectors of R projected back onto Gs to determine what axes of differentiation really matter.

R.proj <- function(Gs,p,vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]]
  rand.vec <-matrix(,vec,n)
  for (i in 1:vec){
    b <- rnorm(n,0,1)
    rand.vec[i,] <- b/(sqrt(sum(b^2)))
  }
  #generate unit length random vectors  
  proj<- function(G,b) t(b) %*% G %*% (b)
  #internal function to do projection
  G.proj <- array(,c(MCMCsamp, m, vec))
  colnames(G.proj) <- dimnames(Gs)[[3]]
  for (i in 1:vec){
    G.proj[,,i]<- t(apply(Gs, 3:4, proj, b = rand.vec[i,]))
  }
  #project each random vector through each MCMC sample of each G
  prs <- cbind(rep(1:m, each = m), 1:m) 
  prs.comp <- prs[prs[,1] < prs[,2], , drop = FALSE] 
  #setting up an index for HPD comparisons
  proj.score <-matrix(,vec,((m^2 - m)/2))
  for (k in 1:vec){
    HPD.int <- HPDinterval(as.mcmc(G.proj[,,k]), prob = p)
    proj.score[k,] <- ifelse(HPD.int[prs.comp[,1],1] > HPD.int[prs.comp[,2],2] | HPD.int[prs.comp[,2],1] > HPD.int[prs.comp[,1],2],1,0) 
  }
  #for a given random vector, examine if the HPD intervals of any pair of G matrices overlap
  vec.score <-cbind(rand.vec, proj.score)
  colnames(vec.score) <- c(1:n, paste(dimnames(Gs)[[3]][prs.comp[, 1]], ".vs.", dimnames(Gs)[[3]][prs.comp[, 2]], sep = ""))
  #collate the random vectors and the outcome of their projection on the G matrices
  sig.vec <- subset(vec.score, rowSums(vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0) 
  #collate just the random vectors that resulted in significant differences in variance
  if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger vec or lower p"); eig.R <- "Na"}
  else{
    eig.R <- eigen(cov(sig.vec[,1:n]))
    rownames(eig.R$vectors) <- dimnames(Gs)[[1]]
    colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
  }  
  #eigen analysis of the R matrix
  list(G.proj = G.proj, vec.score = vec.score, eig.R = eig.R)
}

MCMC.R.proj <- R.proj(Garray, p = 0.95, vec = 10000)

table(rowSums(MCMC.R.proj$vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0 )

#Something like 4424 of 10,000 vectors detected differences among our populations. Will vary slightly every time code is run.
#Save for future analyses
save(MCMC.R.proj, file='MCMC.R.proj.Rdata')
#load('MCMC.R.proj.Rdata')

#Check eigenstructure of R
lapply(MCMC.R.proj$eig.R, round, digits = 3)

#Function to do projection back onto G
proj<- function(G, b) t(b) %*% G %*% (b)

#Genetic variance in each population in the direction of the eigenvectors of R
R.vec.proj <- array(, c(MCMCsamp, m, n))
for (i in 1:n){
  R.vec.proj[,,i] <- t(apply(Garray, 3:4, proj, b = MCMC.R.proj$eig.R$vectors[,i]))
}

#HPD intervals for the genetic variance in each population in the direction of the eigenvectors of R
HPD.R.vec.proj <- array(, c(m, 2, n))
for (i in 1:n){
  HPD.R.vec.proj[,,i] <- HPDinterval(as.mcmc(R.vec.proj[,,i]), prob = 0.95)    
}

#Save HPD intervals and the R.vec.proh
write.table(R.vec.proj, file="R_vector_projection_final.txt", sep="\t")
write.table(HPD.R.vec.proj, file="HPD_int_R_final.txt", sep="\t")

#Generate supplemental figure detailing the amount of genetic variance in the direction of each eigenvector of R for all populations
#First eigenvector
mean<-colMeans(R.vec.proj[,,1])
pop<-Gnames
min<-HPD.R.vec.proj[,1,1]
max<-HPD.R.vec.proj[,2,1]
R1<-data.frame(mean,pop,min,max)

r1<-ggplot(data=R1,aes(y=mean,x=pop))+geom_point()
r1<-r1+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
r1<-r1 + geom_errorbar(aes(ymin=min, ymax=max))
r1<-r1 +theme(axis.title.y=element_blank(), axis.title.x=element_blank())
r1<-r1 + ggtitle(expression(r[1]))

#Second eigenvector
mean<-colMeans(R.vec.proj[,,2])
min<-HPD.R.vec.proj[,1,2]
max<-HPD.R.vec.proj[,2,2]
R2<-data.frame(mean,pop,min,max)

r2<-ggplot(data=R2,aes(y=mean,x=pop))+geom_point()
r2<-r2+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
r2<-r2 + geom_errorbar(aes(ymin=min, ymax=max))
r2<-r2 +theme(axis.title.y=element_blank(), axis.title.x=element_blank())
r2<-r2 + ggtitle(expression(r[2]))

#Third eigenvector
mean<-colMeans(R.vec.proj[,,3])
min<-HPD.R.vec.proj[,1,3]
max<-HPD.R.vec.proj[,2,3]
R3<-data.frame(mean,pop,min,max)

r3<-ggplot(data=R3,aes(y=mean,x=pop))+geom_point()
r3<-r3+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
r3<-r3 + geom_errorbar(aes(ymin=min, ymax=max))
r3<-r3 +theme(axis.title.y=element_blank(), axis.title.x=element_blank())
r3<-r3 + ggtitle(expression(r[3]))

#Fourth eigenvector
mean<-colMeans(R.vec.proj[,,4])
min<-HPD.R.vec.proj[,1,4]
max<-HPD.R.vec.proj[,2,4]
R4<-data.frame(mean,pop,min,max)

r4<-ggplot(data=R4,aes(y=mean,x=pop))+geom_point()
r4<-r4+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
r4<-r4 + geom_errorbar(aes(ymin=min, ymax=max))
r4<-r4 +theme(axis.title.y=element_blank(), axis.title.x=element_blank())
r4<-r4 + ggtitle(expression(r[4]))

y.grob <- textGrob("Genetic variance in direction of eigenvector", gp=gpar(fontsize=15), rot=90)
x.grob <- textGrob("Population", gp=gpar(fontsize=15))

plot<-plot_grid(r1,r2,r3,r4,ncol=2,nrow=2)
grid.arrange(arrangeGrob(plot,left=y.grob, bottom=x.grob))

##############################################################################################################################################################

# Random skewer method (Roff 2012)
# Now moving on to the implementation of Roff's (2012) random skewer method. Briefly, projects same number of vectors onto posterior estimates of G,
# but now calculates the vector correlations of the responses between populations to see if their evolutionary responses to randomly generated skewers differ
# and compare angles to randomly generated G matrices (null hypothesis: no differences)

# This first part is a modification of the function included in Aguirre code. I deconstructed the function to make it more manageable
# Set the number of random vectors. n, m and MCMCsamp already specified at top of code.
vec<-10000
n <- dim(Garray)[[1]]
m <- dim(Garray)[[3]]
MCMCsamp <- dim(Garray)[[4]]
rand.vec <-matrix(,vec,n)
Vec_rsp<-array(rep(1,vec*n*m), dim=c(vec,n,m))
for (i in 1:vec){
  b <- rnorm(n,0,1)
  rand.vec[i,] <- b/(sqrt(sum(b^2)))
}
#generate unit length random vectors  
  
#Below are my modifications/additions.
#I project vectors through each MCMC sample of our Gs

# Matrix to hold vector responses for each of the 10000 estimates for the 5 pops
Vec_rsp<-array(rep(1,MCMCsamp*n*m), dim=c(MCMCsamp,n,m))
# Matrix to hold 10000 estimates of vector corrs. 10 is the number of population comparisons.
vect_cor<-array(rep(1,MCMCsamp,10), dim=c(MCMCsamp,10))
# Matrix to store mean vector correlations (for all the MCMC samples of each G) for each vector. The means and HPD intervals of these is what we report.
mean_vect_cor<-array(rep(1,vec*10), dim=c(vec,10))
   
#for each vector
for (i in 1:vec){
  #Record vector response for each MCMC estimate
  for(j in 1:MCMCsamp){
    for (p in 1:m){
      Vec_rsp[j,,p]<-Garray[,,p,j]%*%rand.vec[i,]
    }
    
    #Calculate the vector correlation, but don't store.
    vect_cor[j,1]<-t(Vec_rsp[j,,1])%*%Vec_rsp[j,,2]/sqrt(t(Vec_rsp[j,,1])%*%Vec_rsp[j,,1]*t(Vec_rsp[j,,2])%*%Vec_rsp[j,,2])
    vect_cor[j,2]<-t(Vec_rsp[j,,1])%*%Vec_rsp[j,,3]/sqrt(t(Vec_rsp[j,,1])%*%Vec_rsp[j,,1]*t(Vec_rsp[j,,3])%*%Vec_rsp[j,,3])
    vect_cor[j,3]<-t(Vec_rsp[j,,1])%*%Vec_rsp[j,,4]/sqrt(t(Vec_rsp[j,,1])%*%Vec_rsp[j,,1]*t(Vec_rsp[j,,4])%*%Vec_rsp[j,,4])
    vect_cor[j,4]<-t(Vec_rsp[j,,1])%*%Vec_rsp[j,,5]/sqrt(t(Vec_rsp[j,,1])%*%Vec_rsp[j,,1]*t(Vec_rsp[j,,5])%*%Vec_rsp[j,,5])
    vect_cor[j,5]<-t(Vec_rsp[j,,2])%*%Vec_rsp[j,,3]/sqrt(t(Vec_rsp[j,,2])%*%Vec_rsp[j,,2]*t(Vec_rsp[j,,3])%*%Vec_rsp[j,,3])
    vect_cor[j,6]<-t(Vec_rsp[j,,2])%*%Vec_rsp[j,,4]/sqrt(t(Vec_rsp[j,,2])%*%Vec_rsp[j,,2]*t(Vec_rsp[j,,4])%*%Vec_rsp[j,,4])
    vect_cor[j,7]<-t(Vec_rsp[j,,2])%*%Vec_rsp[j,,5]/sqrt(t(Vec_rsp[j,,2])%*%Vec_rsp[j,,2]*t(Vec_rsp[j,,5])%*%Vec_rsp[j,,5])
    vect_cor[j,8]<-t(Vec_rsp[j,,3])%*%Vec_rsp[j,,4]/sqrt(t(Vec_rsp[j,,3])%*%Vec_rsp[j,,3]*t(Vec_rsp[j,,4])%*%Vec_rsp[j,,4])
    vect_cor[j,9]<-t(Vec_rsp[j,,3])%*%Vec_rsp[j,,5]/sqrt(t(Vec_rsp[j,,3])%*%Vec_rsp[j,,3]*t(Vec_rsp[j,,5])%*%Vec_rsp[j,,5])
    vect_cor[j,10]<-t(Vec_rsp[j,,4])%*%Vec_rsp[j,,5]/sqrt(t(Vec_rsp[j,,4])%*%Vec_rsp[j,,4]*t(Vec_rsp[j,,5])%*%Vec_rsp[j,,5])
 
    for (a in 1:10){
      mean_vect_cor[i,a]<-mean(vect_cor[,a])}
  
  }
}  

write.table(mean_vect_cor, file="mean_vec_corr.txt", sep="\t")
mean_vect_cor<-read.table(file="mean_vec_corr_observedG_final.txt", sep="\t")

#Get the hpds
HPD<-array(rep(1,10*3),dim=c(10,3))
for (i in 1:10){
  HPD[i,1]<-mean(mean_vect_cor[,i])
}

LQ<-vector()
UQ<-vector()
for (i in 1:10){
  Q<-quantile(mean_vect_cor[,i], seq(0,1,by=0.025))
  LQ<-c(LQ, Q[2])
  UQ<-c(UQ, Q[40])
 }
HPD[,2]<-LQ
HPD[,3]<-UQ

# We generated some randomized G matrices for comparing with observed G matrices for significance testing 2 ways

# First, we randomly assigned phenotypes to individuals in each population and fit 1,000 new G matrices on these null data. This was done in parallel on a server.
# We then estimated differences between these randomized differences using the same analyses conducted on our observed Gs, treating the mean values of 
# genetic variance and covariance of each null matrix as a single estimate for a total of 1000 null comparisons. 

# Secondly, we used the modifications suggested by Morrisey et. al (2019), and added estimates of environmental variance to randomly assigned breeding values
# in each population before fitting G matrices to these data and carrying out analyses.

# The results of both these methods are imported here for comparison. 

Vec_rsp_ran<-array(rep(1,MCMCsamp*n*m), dim=c(MCMCsamp,n,m))
# Matrix to hold 10000 estimates of vector corrs. 10 is the number of population comparisons.
vect_cor_ran<-array(rep(1,MCMCsamp,10), dim=c(MCMCsamp,10))
# Matrix to store mean vector correlations (for all the MCMC samples of each G) for each vector. The means and HPD intervals of these is what we report.
mean_vect_cor_ran<-array(rep(1,vec*10), dim=c(vec,10))

#for each vector
for (i in 1:vec){
  
  #Record vector response for each MCMC estimate
  for(j in 1:MCMCsamp){
    for (p in 1:m){
      Vec_rsp_ran[j,,p]<-rand.Garray[,,p,j]%*%rand.vec[i,]
    }
    
    #Calculate the vector correlation, but don't store.
    vect_cor_ran[j,1]<-t(Vec_rsp_ran[j,,1])%*%Vec_rsp_ran[j,,2]/sqrt(t(Vec_rsp_ran[j,,1])%*%Vec_rsp_ran[j,,1]*t(Vec_rsp_ran[j,,2])%*%Vec_rsp_ran[j,,2])
    vect_cor_ran[j,2]<-t(Vec_rsp_ran[j,,1])%*%Vec_rsp_ran[j,,3]/sqrt(t(Vec_rsp_ran[j,,1])%*%Vec_rsp_ran[j,,1]*t(Vec_rsp_ran[j,,3])%*%Vec_rsp_ran[j,,3])
    vect_cor_ran[j,3]<-t(Vec_rsp_ran[j,,1])%*%Vec_rsp_ran[j,,4]/sqrt(t(Vec_rsp_ran[j,,1])%*%Vec_rsp_ran[j,,1]*t(Vec_rsp_ran[j,,4])%*%Vec_rsp_ran[j,,4])
    vect_cor_ran[j,4]<-t(Vec_rsp_ran[j,,1])%*%Vec_rsp_ran[j,,5]/sqrt(t(Vec_rsp_ran[j,,1])%*%Vec_rsp_ran[j,,1]*t(Vec_rsp_ran[j,,5])%*%Vec_rsp_ran[j,,5])
    vect_cor_ran[j,5]<-t(Vec_rsp_ran[j,,2])%*%Vec_rsp_ran[j,,3]/sqrt(t(Vec_rsp_ran[j,,2])%*%Vec_rsp_ran[j,,2]*t(Vec_rsp_ran[j,,3])%*%Vec_rsp_ran[j,,3])
    vect_cor_ran[j,6]<-t(Vec_rsp_ran[j,,2])%*%Vec_rsp_ran[j,,4]/sqrt(t(Vec_rsp_ran[j,,2])%*%Vec_rsp_ran[j,,2]*t(Vec_rsp_ran[j,,4])%*%Vec_rsp_ran[j,,4])
    vect_cor_ran[j,7]<-t(Vec_rsp_ran[j,,2])%*%Vec_rsp_ran[j,,5]/sqrt(t(Vec_rsp_ran[j,,2])%*%Vec_rsp_ran[j,,2]*t(Vec_rsp_ran[j,,5])%*%Vec_rsp_ran[j,,5])
    vect_cor_ran[j,8]<-t(Vec_rsp_ran[j,,3])%*%Vec_rsp_ran[j,,4]/sqrt(t(Vec_rsp_ran[j,,3])%*%Vec_rsp_ran[j,,3]*t(Vec_rsp_ran[j,,4])%*%Vec_rsp_ran[j,,4])
    vect_cor_ran[j,9]<-t(Vec_rsp_ran[j,,3])%*%Vec_rsp_ran[j,,5]/sqrt(t(Vec_rsp_ran[j,,3])%*%Vec_rsp_ran[j,,3]*t(Vec_rsp_ran[j,,5])%*%Vec_rsp_ran[j,,5])
    vect_cor_ran[j,10]<-t(Vec_rsp_ran[j,,4])%*%Vec_rsp_ran[j,,5]/sqrt(t(Vec_rsp_ran[j,,4])%*%Vec_rsp_ran[j,,4]*t(Vec_rsp_ran[j,,5])%*%Vec_rsp_ran[j,,5])
    
    #Store the mean vector correlation for each vector
    for (a in 1:10){
      mean_vect_cor_ran[i,a]<-mean(vect_cor_ran[,a])
    }
  }  
}


#write.table(mean_vect_cor_ran, file="ran_vec_corr_1-10000.txt", sep="\t")
#Comparison against estimates obtained with Morrissey (2019)'s approach, where G matrices are fit on reconstructed phenotypes reconstituted from randomized breeding values and environmental variance

morrissey_vect_cor<-read.table(file="ranskew.txt", sep="\t")

#Get the hpds
HPD_ran<-array(rep(1,10*3),dim=c(10,3))
for (i in 1:10){
  HPD_ran[i,1]<-mean(morrissey_vect_cor[,i])
}

LQ_ran<-vector()
UQ_ran<-vector()
for (i in 1:10){
  Q_ran<-quantile(morrissey_vect_cor[,i], seq(0,1,by=0.025))
  LQ_ran<-c(LQ_ran, Q_ran[2])
  UQ_ran<-c(UQ_ran, Q_ran[40])
}
HPD_ran[,2]<-LQ_ran
HPD_ran[,3]<-UQ_ran

HPD_ran


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
m <- dim(Garray)[[3]]

MCMCG.kr <- kr.subspace(Garray, vec = rep(2,m))

# Save the object so that it can be looked at later.
save(MCMCG.kr, file='MCMCG.kr.Rdata')
#load('MCMCG.kr.Rdata')

# Load up the table containing the results of Krzanowski's analysis fitted for the Morrissey-adjusted null matrices. 
kr_ran<-read.table('krzan.txt', na.strings="", header=T, sep='\t')

# Load up the MCMC object containing the results of Krzanowski's analysis on the G matrices we generated by randomly assigning phenotypes
load('MCMCG.kr.randG.Rdata')

HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr$MCMC.H.val)), colMeans(MCMCG.kr$MCMC.H.val), HPDinterval(as.mcmc(kr_ran)), colMeans(kr_ran), HPDinterval(as.mcmc(MCMCG.kr.randG$MCMC.H.val)), colMeans(MCMCG.kr.randG$MCMC.H.val))
HPD.H.val

HPD_k<-HPD.H.val[,1:3]
HPD_kr<-HPD.H.val[,4:6]
HPD_kr_ran<-HPD.H.val[,7:9]
vec<-rep(c('h1','h2','h3','h4'),2)
G<-c(rep('obs',4),rep('ran',4))
HPD_K<-rbind(HPD_k,HPD_kr)
HPD_K<-cbind(G,vec,HPD_K)
HPD_K<-as.data.frame(HPD_K)
names(HPD_K)<-c('G','H','min','max','mean')
HPD_K$mean<-as.numeric(as.character(HPD_K$mean))
HPD_K$min<-as.numeric(as.character(HPD_K$min))
HPD_K$max<-as.numeric(as.character(HPD_K$max))

#Generate supplemental figure, comparing the results of the observed matrices to Morrisey's null
k2<-ggplot(data=HPD_K,aes(y=mean,x=H, group=G, colour=G))+geom_point(position=position_dodge(width=0.4), size=3)
k2<-k2 + geom_errorbar(aes(ymin=min, ymax=max), position=position_dodge(width=0.4), size=1,width=0.1)
k2<-k2 +theme(axis.title.y=element_blank(), axis.title.x=element_blank())
k2<-k2 + ylim(0,5)
k2<-k2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k2<-k2 + theme_bw()
#k2<-k2 + scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5))
k2<-k2 + theme(legend.position = c(0.6,0.9))
k2<-k2 + theme(legend.title=element_blank())
k2<-k2 + scale_color_manual(labels=c("Observed","Randomized"),values = c("dodgerblue3","red3"))
k2<-k2+ scale_x_discrete(labels=c(expression(h[1]),expression(h[2]), expression(h[3]), expression(h[4])))
k2<-k2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k2<-k2 + ylab("Eigenvalue")
k2<-k2 + xlab("Eigenvector of H")
k2<-k2 + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
k2<-k2 + theme(panel.border = element_blank())
k2<-k2 + theme(axis.line = element_line(color = 'black'))
k2<-k2 + theme(legend.text = element_text(size=12, face="bold"))
k2<-k2 + theme(axis.text.x = element_text(size=12, face="bold", colour = "black"))
k2<-k2 + theme(axis.text.y = element_text(size=12, face="bold", colour = "black"))

#Generate supplemental figure, comparing the results of the observed matrices to the g matrices we fit from randomized phenotypes

HPD_KR<-rbind(HPD_k,HPD_kr_ran)
HPD_KR<-cbind(G,vec,HPD_KR)
HPD_KR<-as.data.frame(HPD_KR)
names(HPD_KR)<-c('G','H','min','max','mean')
HPD_KR$mean<-as.numeric(as.character(HPD_KR$mean))
HPD_KR$min<-as.numeric(as.character(HPD_KR$min))
HPD_KR$max<-as.numeric(as.character(HPD_KR$max))

k2.1<-ggplot(data=HPD_KR,aes(y=mean,x=H, group=G, colour=G))+geom_point(position=position_dodge(width=0.4), size=3)
k2.1<-k2.1 + geom_errorbar(aes(ymin=min, ymax=max), position=position_dodge(width=0.4), size=1,width=0.1)
k2.1<-k2.1 +theme(axis.title.y=element_blank(), axis.title.x=element_blank())
k2.1<-k2.1 + ylim(0,5)
k2.1<-k2.1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k2.1<-k2.1 + theme_bw()
#k2.1<-k2.1 + scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5))
k2.1<-k2.1 + theme(legend.position = c(0.6,0.9))
k2.1<-k2.1 + theme(legend.title=element_blank())
k2.1<-k2.1 + scale_color_manual(labels=c("Observed","Randomized"),values = c("dodgerblue3","red3"))
k2.1<-k2.1+ scale_x_discrete(labels=c(expression(h[1]),expression(h[2]), expression(h[3]), expression(h[4])))
k2.1<-k2.1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k2.1<-k2.1 + ylab("Eigenvalue")
k2.1<-k2.1 + xlab("Eigenvector of H")
k2.1<-k2.1 + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
k2.1<-k2.1 + theme(panel.border = element_blank())
k2.1<-k2.1 + theme(axis.line = element_line(color = 'black'))

k2.1

round(eigen(MCMCG.kr$avH)$vectors, 3)
round(apply(MCMCG.kr$MCMC.H.theta, 1:2, mean), 1)

#Take a look at how much variance each eignevector is explainging with more detail
val <- matrix(, n, m)
for (i in 1:m){
  avG <- apply(Garray, 1:3, mean)
  val[,i] <- round((cumsum(t(eigen(avG[,,i])$values))/sum(eigen(avG[,,i])$values)*100),3)
}
val

#Let's get the HPD intervals of the angles here.
MCMCG.kr$MCMC.H.theta
angmat<-array(rep(1,10000*5*4), dim=c(10000,5,4))

for (i in 1:4){
  for (a in 1:10000){
    angmat[a,1:5,i]<-MCMCG.kr$MCMC.H.theta[i,1:5,a]
  }
}

anghpd<-array(rep(1,4*2*5), dim=c(4,2,5))

for (i in 1:5){
  LQ<-vector()
  UQ<-vector()
  for (a in 1:4){
    quant<-quantile(angmat[,i,a],seq(0,1,by=0.025))
    LQ<-c(LQ,quant[2])
    UQ<-c(UQ,quant[40]) 
  }
  anghpd[,1,i]<-LQ
  anghpd[,2,i]<-UQ
}

#This line updates the number of eigenvectors to the number required to contain 90% of the variation rather than setting k at 1/2 the number of traits. In or case, 3
MCMCG.kr1 <- kr.subspace(Garray, vec = rep(3,m))

MCMCG.kr.rand1 <- kr.subspace(rand.Garray, vec = rep(3,m))

HPD.H.val1 <- cbind(HPDinterval(as.mcmc(MCMCG.kr1$MCMC.H.val)), colMeans(MCMCG.kr1$MCMC.H.val), HPDinterval(as.mcmc(MCMCG.kr.rand1$MCMC.H.val)), colMeans(MCMCG.kr.rand1$MCMC.H.val))
HPD.H.val1

HPD_k1<-HPD.H.val1[,1:3]
HPD_kr1<-HPD.H.val1[,4:6]
HPD_K1<-rbind(HPD_k1,HPD_kr1)
HPD_K1<-cbind(G,vec,HPD_K1)
HPD_K1<-as.data.frame(HPD_K1)
names(HPD_K1)<-c('G','H','min','max','mean')
HPD_K1$mean<-as.numeric(as.character(HPD_K1$mean))
HPD_K1$min<-as.numeric(as.character(HPD_K1$min))
HPD_K1$max<-as.numeric(as.character(HPD_K1$max))

#Generate supplemental figure
k3<-ggplot(data=HPD_K1,aes(y=mean,x=H, group=G, colour=G))+geom_point(position=position_dodge(width=0.4), size=3)
k3<-k3 + geom_errorbar(aes(ymin=min, ymax=max), position=position_dodge(width=0.4), size=1,width=0.1)
k3<-k3 +theme(axis.title.y=element_blank(), axis.title.x=element_blank())
k3<-k3 + ggtitle('B') +theme(plot.title=element_text(hjust=0.01))
k3<-k3 + scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5))
k3<-k3 + theme(legend.position = "none")
k3<-k3 + theme(legend.title=element_blank())
k3<-k3 + scale_color_manual(labels=c("Observed","Randomized"),values = c("dodgerblue3","red3"))
k3<-k3 + scale_x_discrete(labels=c(expression(h[1]),expression(h[2]), expression(h[3]), expression(h[4])))

y.grob <- textGrob("Eigenvalue", gp=gpar(fontsize=15), rot=90)
x.grob <- textGrob("Eigenvector of H", gp=gpar(fontsize=15))

plotk<-plot_grid(k2,k3,ncol=2,nrow=1)
grid.arrange(arrangeGrob(plotk,left=y.grob, bottom=x.grob))

round(eigen(MCMCG.kr$avH)$vectors, 3)
round(apply(MCMCG.kr$MCMC.H.theta, 1:2, mean), 1)

#########################################################################################################################

#Genetic covariance tensor

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
  dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the jth avG for the eigentensors of posterior mean S
  MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
  dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
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

MCMC.covtensor <- covtensor(Garray)

# Save the object so that it can be looked at later.
save(MCMC.covtensor, file='MCMC.covtensor.Rdata')
#load('MCMC.covtensor.Rdata')

# Load up the table containing the results of the covariance tensor analysis fitted for the Morrissey-adjusted null matrices. 
cov_ran<-read.table('fourth_order.txt', na.strings="", header=T, sep='\t')

# Load up the MCMC object containing the results of the covariance tensor analysis on the G matrices we generated by randomly assigning phenotypes
load('MCMC.covtensor.randG.Rdata')

nnonzero <- min(n*(n+1)/2,m-1)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), prob=0.95), colMeans(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), HPDinterval(as.mcmc(cov_ran[,1:nnonzero]), prob=0.95), colMeans(cov_ran[,1:nnonzero]))
HPD.eT.2<- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), prob=0.95), colMeans(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), HPDinterval(as.mcmc(MCMC.covtensor.randG$MCMC.S.val[,1:nnonzero]), prob=0.95), colMeans(MCMC.covtensor.randG$MCMC.S.val[,1:nnonzero]))

round(HPD.eT.val, 3)
round(HPD.eT.2, 3)

HPD_t<-HPD.eT.val[,1:3]
HPD_tr<-HPD.eT.val[,4:6]
HPD_T<-rbind(HPD_t,HPD_tr)
vec<-rep(c('E1','E2','E3','E4'),2)
G<-c(rep('obs',4),rep('ran',4))
HPD_T<-cbind(G,vec,HPD_T)
HPD_T<-as.data.frame(HPD_T)
names(HPD_T)<-c('G','H','min','max','mean')
HPD_T$mean<-as.numeric(as.character(HPD_T$mean))
HPD_T$min<-as.numeric(as.character(HPD_T$min))
HPD_T$max<-as.numeric(as.character(HPD_T$max))

#Figure 4 Panel A

tens<-ggplot(data=HPD_T,aes(y=mean,x=H, group=G, colour=G))+geom_point(position=position_dodge(width=0.4), size=3)
tens<-tens + geom_errorbar(aes(ymin=min, ymax=max), position=position_dodge(width=0.4), size=1,width=0.2)
tens<-tens + theme_bw()
tens<-tens + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
tens<-tens + labs(x= expression(bold(paste("Eigentensor of ", Sigma))))
tens<-tens + theme(axis.title.x = element_text(size = 12))
tens<-tens + labs(y= expression(bold(paste("Variance ( ", alpha, ") in the direction of E"))))
tens<-tens + theme(axis.title.y = element_text(size = 12))
tens<-tens + scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3))
tens<-tens + theme(legend.position = c(0.4,0.9))
tens<-tens + theme(legend.title=element_blank())
tens<-tens + scale_color_manual(labels=c("Observed","Randomized"),values = c("dodgerblue3","red3"))
tens<-tens + ggtitle(expression(bold("A"))) +theme(plot.title=element_text(hjust=-0.01))
tens<-tens + scale_x_discrete(labels = c(expression(bold(E[1])),expression(bold(E[2])), expression(bold(E[3])), expression(bold(E[4]))))
tens<-tens + theme(axis.text.x = element_text(size=12, face="bold", colour = "black"))
tens<-tens + theme(axis.text.y = element_text(size=12, face="bold", colour = "black"))
tens<-tens + theme(panel.border = element_blank())
tens<-tens + theme(axis.line = element_line())
tens<-tens + theme(legend.text = element_text(size=12, face="bold"))
tens

round(MCMC.covtensor$tensor.summary[1:(n*3),2:dim(MCMC.covtensor$tensor.summary)[2]], 3)

proj<- function(G, b) t(b) %*% G %*% (b)
#Function to do projection to get the amount of variation each populations exhibits

e11 <- c(as.numeric(MCMC.covtensor$tensor.summary[1,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e11.proj <- apply(Garray, 3:4, proj, b = e11)

HPD.e11 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.95)
HPD.e11<-as.data.frame(HPD.e11)
HPD.e11$mean<-rowMeans(e11.proj)
HPD.e11$pop<-Gnames
HPD.e11$pop<-as.factor(HPD.e11$pop)

#Figure 4 Panel B

e11<-ggplot(data=HPD.e11,aes(y=mean,x=pop))+geom_point(size=3)
e11<-e11 + geom_errorbar(aes(ymin=lower, ymax=upper), size=1,width=0.2)
e11<-e11 + theme_bw()
e11<-e11 + theme(axis.title.y=element_blank(), axis.title.x=element_blank())
e11<-e11 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
e11<-e11 + ggtitle(expression(bold(paste("B   ", e[11])))) +theme(plot.title=element_text(hjust=-0.01))
e11<-e11 + theme(axis.text.x = element_text(size=12, face="bold", colour = "black"))
e11<-e11 + theme(axis.text.y = element_text(size=12, face="bold", colour = "black"))
e11<-e11 + theme(panel.border = element_blank())
e11<-e11 + theme(axis.line = element_line())
e11

e21 <- c(as.numeric(MCMC.covtensor$tensor.summary[(n+1),3:dim(MCMC.covtensor$tensor.summary)[2]]))
e21.proj <- apply(Garray, 3:4, proj, b = e21)

HPD.e21 <- HPDinterval(t(as.mcmc(e21.proj)),prob = 0.95)

HPD.e21<-as.data.frame(HPD.e21)
HPD.e21$mean<-rowMeans(e21.proj)
HPD.e21$pop<-Gnames
HPD.e21$pop<-as.factor(HPD.e21$pop)

#Figure 4 Panel C

e21<-ggplot(data=HPD.e21,aes(y=mean,x=pop))+geom_point(size=3)
e21<-e21 + geom_errorbar(aes(ymin=lower, ymax=upper), size=1,width=0.2)
e21<-e21 + theme_bw()
e21<-e21 + theme(axis.title.y=element_blank(), axis.title.x=element_blank())
e21<-e21 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
e21<-e21 + ggtitle(expression(bold(paste("C   ", e[21])))) +theme(plot.title=element_text(hjust=-0.01))
e21<-e21 + theme(axis.text.x = element_text(size=12, face="bold", colour = "black"))
e21<-e21 + theme(axis.text.y = element_text(size=12, face="bold", colour = "black"))
e21<-e21 + theme(panel.border = element_blank())
e21<-e21 + theme(axis.line = element_line())
e21


y.grob <- textGrob(expression(bold(paste("Genetic variance ( ", lambda, ") in the direction of e"))), gp=gpar(fontsize=12), rot=90)
x.grob <- textGrob(expression(bold("Population", gp=gpar(fontsize=12))))

plot_e12<-plot_grid(e11,e21,ncol=1,nrow=2)

plot_tens<-plot_grid(tens,grid.arrange(arrangeGrob(plot_e12,left=y.grob, bottom=x.grob)),ncol=2,nrow=1)

lapply(MCMC.covtensor$eig.R, round, digits = 3)
MCMC.covtensor$tensor.summary
MCMC.covtensor$eTmat
MCMC.covtensor$av.S
MCMC.covtensor$MCMC.S.val

##############################################################################################################################################################

# Flury's Hierarchical Method 

library(cpc)

# First we need the sample sizes of each group, which are:
n_sam<-c(167,162,153,156,162)

Gmats<- array(,c(n,n,m))
dimnames(Gmats) <- list(traitnames,traitnames,Gnames)
for (i in 1:5){
  Gmats[,,i]<-Gsum[,,1,i]
}

# Estimate stepwise cpc and FG cpc to diagonalize several covariance matrices
step_cpc<-stepwisecpc(Gmats, n_sam)
FG_cpc<-FG(Gmats, n_sam)

# Flury's hieracarchy test.
# Will look at both approaches to matrix diagonalization
Flury_res_step<-flury.test(Gmats, n_sam, B=step_cpc$B, p = 4, qmax= 2)
Flury_res_fg<-flury.test(Gmats, n_sam, B=FG_cpc$B, p = 4, qmax= 2)

#Run for all posterior estimates of G.
#Just work with a dataframe, so 7 columns, the first one will be the label for MCMCsample

Fl_stp_HPD <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(Fl_stp_HPD) <- c('Model', 'Chi.square', 'DF','Chi2.div.df','AIC','No.of.CPCs','ID')
for (i in 1:10000){
  step.cpc<-stepwisecpc(Garray[,,,i], n_sam)
  res<-flury.test(Garray[,,,i], n_sam, B=step.cpc$B, p = 4, qmax= 2)
  res<-as.data.frame(res)
  res$ID<-rep(i,6)
  Fl_stp_HPD<-rbind(Fl_stp_HPD,res)
}

Fl_stp_HPD$Model<-as.factor(Fl_stp_HPD$Model)
str(Fl_stp_HPD)

Fl_fg_HPD <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(Fl_fg_HPD) <- c('Model', 'Chi.square', 'DF','Chi2.div.df','AIC','No.of.CPCs','ID')
for (i in 1:10000){
  fg.cpc<-FG(Garray[,,,i], n_sam)
  res<-flury.test(Garray[,,,i], n_sam, B=fg.cpc$B, p = 4, qmax= 2)
  res<-as.data.frame(res)
  res$ID<-rep(i,6)
  Fl_fg_HPD<-rbind(Fl_fg_HPD,res)
}

Fl_fg_HPD$Model<-as.factor(Fl_fg_HPD$Model)
str(Fl_fg_HPD)

aggregate(x=Fl_stp_HPD$AIC,
          by = list(Fl_stp_HPD$Model),
          FUN = mean)

aggregate(x=Fl_fg_HPD$AIC,
          by = list(Fl_fg_HPD$Model),
          FUN = mean)
