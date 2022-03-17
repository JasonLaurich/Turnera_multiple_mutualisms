## This file runs selection skewers through the G matrices of our 5 largest populations to determine how differences in G could lead to divergence among them.
## It also estimates the degree to which our largest populations experience constraint or facilitation in their responses to multivariate selection through the
## calculation of R values.
## It generates Figure A4 (R-values) and Figure 5 (Response to selection)

#Contains code for generation of figure 5 and supplemental figure 4

library(MCMCglmm)
library(cowplot)
library(grid)
library(gridExtra)
library(ggplot2)
library(forcats)

#Will require changing directory
setwd("C:/Users/Jason/Desktop/Data/TurnChap1/Dryad Files")

#this I need output arranged in [trait,trait,pop,mcmc] fasion
make.array <- function(model.list,ntraits,nsamples,whichtype){ 
  mcmc.array <- array(,dim=c(ntraits,ntraits,length(model.list),nsamples))
  if(whichtype=="DS"){
    for(i in 1:length(model.list)){
      mcmc.array[,,i,] <- model.list[[i]]$G
    }
  } else if(whichtype=="animal"){
    gcols <- grep(colnames(model.list[[1]]$VCV),pattern="animal")
    for(i in 1:length(model.list)){
      for(j in 1:nsamples){
        mcmc.array[,,i,j] <- matrix(model.list[[i]]$VCV[j,gcols],ncol=ntraits)
      }
    }
  } else {print("whichtype should equal \"DS\" or \"Animal\"")}
  return(mcmc.array)
}

# A loop to extract G matrix data for each posterior estimate of each population. Full matrices assigned to GM[i], posterior estimates with real correlation
# assigned to G[i], and posterior estimates with correlations set to zero assigned to Gnull[i].

for (i in 1:5){
  load(paste('GMpop_final',i,".Rdata",sep=""))
  assign(paste("GM",i,sep=""),GM)
  
  G<-make.array(list(GM) , 4 , nrow(GM$VCV) , whichtype="animal" )
  assign(paste("G", i, sep=""),G)
  
  Gnull<-make.array(list(GM) , 4 , nrow(GM$VCV) , whichtype="animal" )
  Gnull[1,2:4,,]<-0
  Gnull[2,c(1,3:4),,]<-0
  Gnull[3,c(1:2,4),,]<-0
  Gnull[4,1:3,,]<-0
  assign(paste("Gnull", i, sep=""),Gnull)

}

load("GM_Jam.Rdata")

MCMCarray <- array(,c(10000,(4^2)))
MCMCarray<-GMpr$VCV

GJam<-array(,c(4,4,1,10000))

for (j in 1:10000){
  GJam[,,,j] <- matrix(MCMCarray[j,2:17],ncol= 4)
}

GnullJam<-GJam
GnullJam[1,2:4,,]<-0
GnullJam[2,c(1,3:4),,]<-0
GnullJam[3,c(1:2,4),,]<-0
GnullJam[4,1:3,,]<-0

#Next we define our selection scenarios

# Weak universal positive selection scenario
sel_1<-c(0.1,0.1,0.1,0.1)

# Weak positive selection on dispersal and defence (elaiosome mass and EFN), negative selection on FN and stigma-anther separation
sel_2<-c(0.1,0.1,-0.1,-0.1)

# Weak negative selection against out-crossing (stigma-anther separation), positive selection on dispersal and defence, with neutral selection on floral nectar
sel_3<-c(0.1,0.1,0,-0.1)

# Calculate R for the first selection scenario for each population. 

Glist<-list(GJam,G1,G2,G3,G4,G5)
Gnulllist<-list(GnullJam,Gnull1,Gnull2,Gnull3,Gnull4,Gnull5)

for (i in 1:6){
  R<-vector()
  G<-Glist[[i]]
  R0<-vector()
  Gnull<-Gnulllist[[i]]
  for (j in 1:10000){
    r<-G[,,,j]%*%sel_1
    r<-t(r)
    R<-rbind(R,r) # don't save individual posterior responses to selection
  }
  
  R_sel<-array(dim=c(4,3)) 
  for (k in 1:4){
    R_sel[k,1]<-mean(R[,k])
    LQ<-vector()
    UQ<-vector()
    Q<-quantile(R[,k], seq(0,1,by=0.025))
    LQ<-c(LQ, Q[2])
    UQ<-c(UQ, Q[40])
    R_sel[k,2]<-LQ
    R_sel[k,3]<-UQ
  }
  assign(paste("R_sel1_", i, sep=""),R_sel) #save hpd intervals for response to selection for each population.
  
  for (j in 1:10000){
    r<-Gnull[,,,j]%*%sel_1
    r<-t(r)
    R0<-rbind(R0,r) # don't save the null responses. 
  }
  R_v<-R/R0
  
  R_tab<-array(dim=c(4,3)) 
  for (k in 1:4){
    R_tab[k,1]<-mean(R_v[,k])
    LQ<-vector()
    UQ<-vector()
    Q<-quantile(R_v[,k], seq(0,1,by=0.025))
    LQ<-c(LQ, Q[2])
    UQ<-c(UQ, Q[40])
    R_tab[k,2]<-LQ
    R_tab[k,3]<-UQ
  }
  assign(paste("R_tab1_", i, sep=""),R_tab) # save the hpd intervals for the R values of each population. 
}


#Make summary tables to generate figures with.
pop <- c(rep("Jamaica", 4), rep("Brown's Town",4), rep("Cave",4), rep("Mosquito Cove",4), rep("Murdock",4), rep("St. Ann's Bay",4))
pop<-as.factor(pop)
pop<-fct_relevel(pop, "Jamaica")

trait<-rep(c("Elaiosome mass", "EFN", "FN", "Stigma-anther separation"),6)
trait <- factor(trait, levels = c("Elaiosome mass", "EFN", "FN", "Stigma-anther separation"))

R_sel1<-rbind(R_sel1_1, R_sel1_2,R_sel1_3,R_sel1_4,R_sel1_5,R_sel1_6)
R_sel1<-as.data.frame(R_sel1)
names(R_sel1)<-c("mean","LQ","UQ")
R_sel1<-cbind(R_sel1,pop,trait)

R_sel1$pop<-as.factor(R_sel1$pop)
R_sel1$pop<-fct_relevel(R_sel1$pop, "Jamaica")

R_tab1<-rbind(R_tab1_1, R_tab1_2,R_tab1_3,R_tab1_4,R_tab1_5,R_tab1_6)
R_tab1<-as.data.frame(R_tab1)
names(R_tab1)<-c("mean","LQ","UQ")
R_tab1<-cbind(R_tab1, pop, trait)

# Weak positive selection on dispersal and defence (elaiosome mass and EFN), negative selection on FN and stigma-anther separation
for (i in 1:6){
  R<-vector()
  G<-Glist[[i]]
  R0<-vector()
  Gnull<-Gnulllist[[i]]
  for (j in 1:10000){
    r<-G[,,,j]%*%sel_2
    r<-t(r)
    R<-rbind(R,r)
  }
  
  R_sel<-array(dim=c(4,3)) 
  for (k in 1:4){
    R_sel[k,1]<-mean(R[,k])
    LQ<-vector()
    UQ<-vector()
    Q<-quantile(R[,k], seq(0,1,by=0.025))
    LQ<-c(LQ, Q[2])
    UQ<-c(UQ, Q[40])
    R_sel[k,2]<-LQ
    R_sel[k,3]<-UQ
  }
  assign(paste("R_sel2_", i, sep=""),R_sel)
  
  for (j in 1:10000){
    r<-Gnull[,,,j]%*%sel_2
    r<-t(r)
    R0<-rbind(R0,r)
  }
  R_v<-R/R0
  
  R_tab<-array(dim=c(4,3)) 
  for (k in 1:4){
    R_tab[k,1]<-mean(R_v[,k])
    LQ<-vector()
    UQ<-vector()
    Q<-quantile(R_v[,k], seq(0,1,by=0.025))
    LQ<-c(LQ, Q[2])
    UQ<-c(UQ, Q[40])
    R_tab[k,2]<-LQ
    R_tab[k,3]<-UQ
  }
  assign(paste("R_tab2_", i, sep=""),R_tab)
}

R_sel2<-rbind(R_sel2_1, R_sel2_2,R_sel2_3,R_sel2_4,R_sel2_5,R_sel2_6)
R_sel2<-as.data.frame(R_sel2)
names(R_sel2)<-c("mean","LQ","UQ")
R_sel2<-cbind(R_sel2,pop,trait)

R_tab2<-rbind(R_tab2_1, R_tab2_2,R_tab2_3,R_tab2_4,R_tab2_5,R_tab2_6)
R_tab2<-as.data.frame(R_tab2)
names(R_tab2)<-c("mean","LQ","UQ")
R_tab2<-cbind(R_tab2, pop, trait)

# Weak negative selection against out-crossing (stigma-anther separation), positive selection on dispersal and defence, with neutral selection on floral nectar
for (i in 1:6){
  R<-vector()
  G<-Glist[[i]]
  R0<-vector()
  Gnull<-Gnulllist[[i]]
  for (j in 1:10000){
    r<-G[,,,j]%*%sel_3
    r<-t(r)
    R<-rbind(R,r)
  }
  
  R_sel<-array(dim=c(4,3)) 
  for (k in 1:4){
    R_sel[k,1]<-mean(R[,k])
    LQ<-vector()
    UQ<-vector()
    Q<-quantile(R[,k], seq(0,1,by=0.025))
    LQ<-c(LQ, Q[2])
    UQ<-c(UQ, Q[40])
    R_sel[k,2]<-LQ
    R_sel[k,3]<-UQ
  }
  assign(paste("R_sel3_", i, sep=""),R_sel)
  
  for (j in 1:10000){
    r<-Gnull[,,,j]%*%sel_3
    r<-t(r)
    R0<-rbind(R0,r)
  }
  R_v<-R/R0
  
  R_tab<-array(dim=c(4,3)) 
  for (k in 1:4){
    R_tab[k,1]<-mean(R_v[,k])
    LQ<-vector()
    UQ<-vector()
    Q<-quantile(R_v[,k], seq(0,1,by=0.025))
    LQ<-c(LQ, Q[2])
    UQ<-c(UQ, Q[40])
    R_tab[k,2]<-LQ
    R_tab[k,3]<-UQ
  }
  assign(paste("R_tab3_", i, sep=""),R_tab)
}

R_sel3<-rbind(R_sel3_1, R_sel3_2,R_sel3_3,R_sel3_4,R_sel3_5,R_sel3_6)
R_sel3<-as.data.frame(R_sel3)
names(R_sel3)<-c("mean","LQ","UQ")
R_sel3<-cbind(R_sel3,pop,trait)

R_tab3<-rbind(R_tab3_1, R_tab3_2,R_tab3_3,R_tab3_4,R_tab3_5,R_tab3_6)
R_tab3<-as.data.frame(R_tab3)
names(R_tab3)<-c("mean","LQ","UQ")
R_tab3<-cbind(R_tab3, pop, trait)
R_tab3<-R_tab3[-c(3,7,11,15,19,23),]
R_tab3<-droplevels(R_tab3)

#Generate Figure A4 : estimates of constraint and facilitation through G (R-values)

p1<-ggplot(data=R_tab1,aes(y=mean,x=trait, group=pop, colour=pop)) 
p1<-p1+geom_point(position=position_dodge(width=0.4), size=3)
p1<-p1+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
p1<-p1+ geom_errorbar(aes(ymin=LQ, ymax=UQ), position=position_dodge(width=0.4), size=1,width=0.1)
p1<-p1+theme(axis.title.x=element_blank())
p1<-p1+labs(y= "R-value")
p1<-p1 +ggtitle('A') +theme(plot.title=element_text(hjust=-0.01))
p1<-p1 +theme(legend.position = "none")
p1<-p1 +scale_color_manual(values = c("black","red3","dodgerblue3","chartreuse4","darkorange1","darkorchid4"))
p1<-p1 +geom_hline(yintercept=1, linetype="dashed")
p1<-p1 +ylim(-1.5, 3.5)
p1<-p1 +scale_x_discrete(labels=c(expression(paste("Elaiosome ", italic(m))), "EFN", "FN", "S-A separation"))
p1

p2<-ggplot(data=R_tab2,aes(y=mean,x=trait, group=pop, colour=pop)) 
p2<-p2+geom_point(position=position_dodge(width=0.4), size=3)
p2<-p2+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
p2<-p2+ geom_errorbar(aes(ymin=LQ, ymax=UQ), position=position_dodge(width=0.4), size=1,width=0.1)
p2<-p2+theme(axis.title.x=element_blank())
p2<-p2+theme(axis.title.y=element_blank())
p2<-p2+theme(legend.position = "none")
p2<-p2 +ggtitle('B') +theme(plot.title=element_text(hjust=-0.01))
p2<-p2 +scale_color_manual(values = c("black","red3","dodgerblue3","chartreuse4","darkorange1","darkorchid4"))
p2<-p2 +geom_hline(yintercept=1, linetype="dashed")
p2<-p2 +ylim(-1.5, 3.5)
p2<-p2 +scale_x_discrete(labels=c(expression(paste("Elaiosome ", italic(m))), "EFN", "FN", "S-A separation"))
p2

p3<-ggplot(data=R_tab3,aes(y=mean,x=trait, group=pop, colour=pop)) 
p3<-p3+geom_point(position=position_dodge(width=0.4), size=3)
p3<-p3+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
p3<-p3+ geom_errorbar(aes(ymin=LQ, ymax=UQ), position=position_dodge(width=0.4), size=1,width=0.1)
p3<-p3+theme(axis.title.x=element_blank())
p3<-p3+theme(axis.title.y=element_blank())
p3<-p3+theme(legend.title = element_blank())
p3<-p3+theme(legend.position = c(0.15,0.90))
p3<-p3 +ggtitle('C') +theme(plot.title=element_text(hjust=-0.01))
p3<-p3 +scale_color_manual(values = c("black","red3","dodgerblue3","chartreuse4","darkorange1","darkorchid4"))
p3<-p3 +geom_hline(yintercept=1, linetype="dashed")
p3<-p3 +ylim(-1.5, 3.5)
p3<-p3 +scale_x_discrete(labels=c(expression(paste("Elaiosome ", italic(m))), "EFN", "S-A separation"))
p3

x.grob <- textGrob("Trait", gp=gpar(fontsize=15))
plot_R<-plot_grid(p1,p2,p3, ncol=3,nrow=1)

plot_R2<-grid.arrange(arrangeGrob(plot_R, bottom=x.grob))

#Generate Figure 5: response to selection scenarios

sel1<-ggplot(data=R_sel1,aes(y=mean,x=trait, group=pop, colour=pop)) 
sel1<-sel1+geom_point(position=position_dodge(width=0.4), size=3)
sel1<-sel1+ geom_errorbar(aes(ymin=LQ, ymax=UQ), position=position_dodge(width=0.4), size=1,width=0.1)
sel1<-sel1+theme(axis.title.x=element_blank())
sel1<-sel1+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
sel1<-sel1+labs(y= "Response to selection")
sel1<-sel1+theme(axis.title.x=element_blank())
sel1<-sel1+ggtitle('A') +theme(plot.title=element_text(hjust=-0.01))
sel1<-sel1+theme(legend.position = "none")
sel1<-sel1+scale_color_manual(values = c("black","red3","dodgerblue3","chartreuse4","darkorange1","darkorchid4"))
sel1<-sel1+geom_hline(yintercept=0, linetype="dashed")
sel1<-sel1+ylim(-0.2, 0.2)
sel1<-sel1+scale_x_discrete(labels=c(expression(paste("Ela ", italic(m))), "EFN", "FN", "S-A sep"))
sel1

sel2<-ggplot(data=R_sel2,aes(y=mean,x=trait, group=pop, colour=pop)) 
sel2<-sel2+geom_point(position=position_dodge(width=0.4), size=3)
sel2<-sel2+ geom_errorbar(aes(ymin=LQ, ymax=UQ), position=position_dodge(width=0.4), size=1,width=0.1)
sel2<-sel2+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
sel2<-sel2+theme(axis.title.x=element_blank())
sel2<-sel2+theme(axis.title.y=element_blank())
sel2<-sel2+ggtitle('B') +theme(plot.title=element_text(hjust=-0.01))
sel2<-sel2+theme(legend.position = "none")
sel2<-sel2+scale_color_manual(values = c("black","red3","dodgerblue3","chartreuse4","darkorange1","darkorchid4"))
sel2<-sel2+geom_hline(yintercept=0, linetype="dashed")
sel2<-sel2+ylim(-0.2, 0.2)
sel2<-sel2+scale_x_discrete(labels=c(expression(paste("Ela ", italic(m))), "EFN", "FN", "S-A sep"))
sel2

sel3<-ggplot(data=R_sel3,aes(y=mean,x=trait, group=pop, colour=pop)) 
sel3<-sel3+geom_point(position=position_dodge(width=0.4), size=3)
sel3<-sel3+ geom_errorbar(aes(ymin=LQ, ymax=UQ), position=position_dodge(width=0.4), size=1,width=0.1)
sel3<-sel3+theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), axis.line=element_line(color='black'))
sel3<-sel3+theme(axis.title.x=element_blank())
sel3<-sel3+theme(axis.title.y=element_blank())
sel3<-sel3+ggtitle('C') +theme(plot.title=element_text(hjust=-0.01))
sel3<-sel3+theme(legend.position = c(0.7,0.88))
sel3<-sel3+theme(legend.title = element_blank())
sel3<-sel3+scale_color_manual(values = c("black","red3","dodgerblue3","chartreuse4","darkorange1","darkorchid4"))
sel3<-sel3+geom_hline(yintercept=0, linetype="dashed")
sel3<-sel3+ylim(-0.2, 0.2)
sel3<-sel3+scale_x_discrete(labels=c(expression(paste("Ela ", italic(m))), "EFN", "FN", "S-A sep"))
sel3

plot_sel<-plot_grid(sel1,sel2,sel3, ncol=3,nrow=1)

plot_sel_t<-grid.arrange(arrangeGrob(plot_sel, bottom=x.grob))

plot_sel_t
