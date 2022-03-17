library(MCMCglmm)
library(foreach)
library(doParallel)
library(gdata)
library(matrixcalc)

setwd("C:/Users/Jason/Desktop/Data/TurnChap1/Dryad Files")

# Load randomized G matrices
BT<-read.table('BT_ranG_final.txt', header= T)
Ca<-read.table('Ca_ranG_final.txt', header= T)
MC<-read.table('MC_ranG_final.txt', header= T)
Mu<-read.table('Mu_ranG_final.txt', header= T)
SA<-read.table('SA_ranG_final.txt', header= T)
  
#number of MCMC samples
MCMCsamp <- 1000
#number of traits 
n <- 4
#number of matrices to compare
m <- 5 
#trait names
traitnames <- c("ela","efn","fn","herk") 
#matrix labels 
Gnames <- c("BT","Ca","MC","Mu","SA")
  
Garray <- array(,c(n,n,m,MCMCsamp))
dimnames(Garray) <- list(traitnames,traitnames,Gnames)

# Load the randomized G matrices into Garray
  
  for (i in 1:1000){
    Garray[1,1,1,i]<-BT[i,1] 
    Garray[1,2,1,i]<-BT[i,2]
    Garray[2,1,1,i]<-BT[i,2]
    Garray[1,3,1,i]<-BT[i,3]
    Garray[3,1,1,i]<-BT[i,3]
    Garray[1,4,1,i]<-BT[i,4]
    Garray[4,1,1,i]<-BT[i,4]
    Garray[2,2,1,i]<-BT[i,5]
    Garray[2,3,1,i]<-BT[i,6]
    Garray[3,2,1,i]<-BT[i,6]
    Garray[2,4,1,i]<-BT[i,7]
    Garray[4,2,1,i]<-BT[i,7]
    Garray[3,3,1,i]<-BT[i,8]
    Garray[3,4,1,i]<-BT[i,9]
    Garray[4,3,1,i]<-BT[i,9]
    Garray[4,4,1,i]<-BT[i,10]
    
    Garray[1,1,2,i]<-Ca[i,1] 
    Garray[1,2,2,i]<-Ca[i,2]
    Garray[2,1,2,i]<-Ca[i,2]
    Garray[1,3,2,i]<-Ca[i,3]
    Garray[3,1,2,i]<-Ca[i,3]
    Garray[1,4,2,i]<-Ca[i,4]
    Garray[4,1,2,i]<-Ca[i,4]
    Garray[2,2,2,i]<-Ca[i,5]
    Garray[2,3,2,i]<-Ca[i,6]
    Garray[3,2,2,i]<-Ca[i,6]
    Garray[2,4,2,i]<-Ca[i,7]
    Garray[4,2,2,i]<-Ca[i,7]
    Garray[3,3,2,i]<-Ca[i,8]
    Garray[3,4,2,i]<-Ca[i,9]
    Garray[4,3,2,i]<-Ca[i,9]
    Garray[4,4,2,i]<-Ca[i,10]
    
    Garray[1,1,3,i]<-MC[i,1] 
    Garray[1,2,3,i]<-MC[i,2]
    Garray[2,1,3,i]<-MC[i,2]
    Garray[1,3,3,i]<-MC[i,3]
    Garray[3,1,3,i]<-MC[i,3]
    Garray[1,4,3,i]<-MC[i,4]
    Garray[4,1,3,i]<-MC[i,4]
    Garray[2,2,3,i]<-MC[i,5]
    Garray[2,3,3,i]<-MC[i,6]
    Garray[3,2,3,i]<-MC[i,6]
    Garray[2,4,3,i]<-MC[i,7]
    Garray[4,2,3,i]<-MC[i,7]
    Garray[3,3,3,i]<-MC[i,8]
    Garray[3,4,3,i]<-MC[i,9]
    Garray[4,3,3,i]<-MC[i,9]
    Garray[4,4,3,i]<-MC[i,10]
    
    Garray[1,1,4,i]<-Mu[i,1] 
    Garray[1,2,4,i]<-Mu[i,2]
    Garray[2,1,4,i]<-Mu[i,2]
    Garray[1,3,4,i]<-Mu[i,3]
    Garray[3,1,4,i]<-Mu[i,3]
    Garray[1,4,4,i]<-Mu[i,4]
    Garray[4,1,4,i]<-Mu[i,4]
    Garray[2,2,4,i]<-Mu[i,5]
    Garray[2,3,4,i]<-Mu[i,6]
    Garray[3,2,4,i]<-Mu[i,6]
    Garray[2,4,4,i]<-Mu[i,7]
    Garray[4,2,4,i]<-Mu[i,7]
    Garray[3,3,4,i]<-Mu[i,8]
    Garray[3,4,4,i]<-Mu[i,9]
    Garray[4,3,4,i]<-Mu[i,9]
    Garray[4,4,4,i]<-Mu[i,10]
    
    Garray[1,1,5,i]<-SA[i,1] 
    Garray[1,2,5,i]<-SA[i,2]
    Garray[2,1,5,i]<-SA[i,2]
    Garray[1,3,5,i]<-SA[i,3]
    Garray[3,1,5,i]<-SA[i,3]
    Garray[1,4,5,i]<-SA[i,4]
    Garray[4,1,5,i]<-SA[i,4]
    Garray[2,2,5,i]<-SA[i,5]
    Garray[2,3,5,i]<-SA[i,6]
    Garray[3,2,5,i]<-SA[i,6]
    Garray[2,4,5,i]<-SA[i,7]
    Garray[4,2,5,i]<-SA[i,7]
    Garray[3,3,5,i]<-SA[i,8]
    Garray[3,4,5,i]<-SA[i,9]
    Garray[4,3,5,i]<-SA[i,9]
    Garray[4,4,5,i]<-SA[i,10]
  }
  
# Ok so 1 to 10 are the vector correlations for each population pair, 11 to 14 are the H values for the 4 eigenvectors (Krzanowski's).
# 4 more for the 4th order covariance tensor
ranskew_stats_ranG<-array(,c(1000,10))
k_stats_ranG<-array(,c(1000,4))
four_stats_ranG<-array(,c(1000,4))
  
Vec_rsp_ran<-array(rep(1,1000*n*m), dim=c(1000,n,m))
# Matrix to hold 1000 estimates of vector corrs. 10 is the number of population comparisons.
vect_cor_ran<-array(rep(1,1000,10), dim=c(1000,10))
# Matrix to store mean vector correlations (for all the MCMC samples of each G) for each vector. The means and HPD intervals of these is what we report.
mean_vect_cor_ran<-array(rep(1,1000*10), dim=c(1000,10))
  
 
# Now we do Krzanowski's analysis
    
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
    AA <- array(, c(4, 4, 5))  
    for (k in 1:5){
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
  for (z in 1:MCMCsamp){
    kr <- h(Gs[,,,z], v = vec)
    MCMC.H[,,z] <- kr$H
    MCMC.AA[,,,z] <- kr$AA
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
  
MCMCG.kr.rand <- kr.subspace(Garray, vec = rep(2,m))

K.ran<-colMeans(MCMCG.kr.rand$MCMC.H.val)
K.ran<-as.vector(K.ran)
    
    print("Krzan summed")
    
    k_stats[x,1:4]<-K.ran
    
    write.table(k_stats_ranG, file=paste("k_stats_ranG",y,".txt"), sep="\t")
    
    
    #Generate multivariate normal selection vectors. 
    vec<-1000
    n <- dim(Garray)[[1]]
    m <- dim(Garray)[[3]]
    MCMCsamp <- dim(Garray)[[4]]
    rand.vec <-matrix(,vec,n)
    Vec_rsp<-array(rep(1,vec*n*m), dim=c(vec,n,m))
    for (i in 1:vec){
      b <- rnorm(n,0,1)
      rand.vec[i,] <- b/(sqrt(sum(b^2)))
    }
    
    print("vectors generated")
    
    #for each vector
    for (j in 1:1000){
      
      #Record vector response for each MCMC estimate
      for(k in 1:1000){
        for (p in 1:m){
          Vec_rsp_ran[k,,p]<-Garray[,,p,k]%*%rand.vec[k,]
        }
        
        #Calculate the vector correlation, but don't store.
        vect_cor_ran[k,1]<-t(Vec_rsp_ran[k,,1])%*%Vec_rsp_ran[k,,2]/sqrt(t(Vec_rsp_ran[k,,1])%*%Vec_rsp_ran[k,,1]*t(Vec_rsp_ran[k,,2])%*%Vec_rsp_ran[k,,2])
        vect_cor_ran[k,2]<-t(Vec_rsp_ran[k,,1])%*%Vec_rsp_ran[k,,3]/sqrt(t(Vec_rsp_ran[k,,1])%*%Vec_rsp_ran[k,,1]*t(Vec_rsp_ran[k,,3])%*%Vec_rsp_ran[k,,3])
        vect_cor_ran[k,3]<-t(Vec_rsp_ran[k,,1])%*%Vec_rsp_ran[k,,4]/sqrt(t(Vec_rsp_ran[k,,1])%*%Vec_rsp_ran[k,,1]*t(Vec_rsp_ran[k,,4])%*%Vec_rsp_ran[k,,4])
        vect_cor_ran[k,4]<-t(Vec_rsp_ran[k,,1])%*%Vec_rsp_ran[k,,5]/sqrt(t(Vec_rsp_ran[k,,1])%*%Vec_rsp_ran[k,,1]*t(Vec_rsp_ran[k,,5])%*%Vec_rsp_ran[k,,5])
        vect_cor_ran[k,5]<-t(Vec_rsp_ran[k,,2])%*%Vec_rsp_ran[k,,3]/sqrt(t(Vec_rsp_ran[k,,2])%*%Vec_rsp_ran[k,,2]*t(Vec_rsp_ran[k,,3])%*%Vec_rsp_ran[k,,3])
        vect_cor_ran[k,6]<-t(Vec_rsp_ran[k,,2])%*%Vec_rsp_ran[k,,4]/sqrt(t(Vec_rsp_ran[k,,2])%*%Vec_rsp_ran[k,,2]*t(Vec_rsp_ran[k,,4])%*%Vec_rsp_ran[k,,4])
        vect_cor_ran[k,7]<-t(Vec_rsp_ran[k,,2])%*%Vec_rsp_ran[k,,5]/sqrt(t(Vec_rsp_ran[k,,2])%*%Vec_rsp_ran[k,,2]*t(Vec_rsp_ran[k,,5])%*%Vec_rsp_ran[k,,5])
        vect_cor_ran[k,8]<-t(Vec_rsp_ran[k,,3])%*%Vec_rsp_ran[k,,4]/sqrt(t(Vec_rsp_ran[k,,3])%*%Vec_rsp_ran[k,,3]*t(Vec_rsp_ran[k,,4])%*%Vec_rsp_ran[k,,4])
        vect_cor_ran[k,9]<-t(Vec_rsp_ran[k,,3])%*%Vec_rsp_ran[k,,5]/sqrt(t(Vec_rsp_ran[k,,3])%*%Vec_rsp_ran[k,,3]*t(Vec_rsp_ran[k,,5])%*%Vec_rsp_ran[k,,5])
        vect_cor_ran[k,10]<-t(Vec_rsp_ran[k,,4])%*%Vec_rsp_ran[k,,5]/sqrt(t(Vec_rsp_ran[k,,4])%*%Vec_rsp_ran[k,,4]*t(Vec_rsp_ran[k,,5])%*%Vec_rsp_ran[k,,5])
        
        #Store the mean vector correlation for each vector in the table
        for (a in 1:10){
          ranskew_stats[x,a]<-mean(vect_cor_ran[,a])
        }
      }  
    }
    
    write.table(ranskew_stats_ranG, file=paste("ranskew_stats_ranG",y,".txt"), sep="\t")
    
    print("end of random skewers")
    
    
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
    
    nnonzero <- min(n*(n+1)/2,m-1)
    
    MCMC.covtensor.rand <- covtensor(Garray)
    
    colMeans(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero])
    
    Cov.ran<-colMeans(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero])
    Cov.ran<-as.vector(Cov.ran)
    
    four_stats[x,1:4]<-Cov.ran
    
    write.table(four_stats_ranG, file=paste("four_stats_ranG",y,".txt"), sep="\t")
    
  } # end of looping through individual random G matrix generation based on MCMC estiamtes and summary stat generation
  
} # end of parallel

stopImplicitCluster()


