# This file generates Figure 3 : a depiction of variation in phenotypic trait means among our 17 population of Turnera ulmifolia. 

library(gridExtra)
library(ggplot2)

#Will require changing directory
setwd("C:/Users/Jason/Desktop/Dryad 1")

# We're going to start making the MCMC objects for the Jamaican metapopulation

# Upload data files and reorganize them 

pheno<- read.table("pheno.txt",stringsAsFactors=F,header=F)#phenotypes
names(pheno)<-c('id','site','mat','elave','efn','fn','herk','growth','efn_nec')

data<-pheno[!(pheno$site=='Extra'),]
data<-data[(data$efn>0),]
data<-droplevels(data)

dataG<-data[(!is.na(data$efn)&!is.na(data$elave)&!is.na(data$fn)&!is.na(data$herk)),]

dataG$lnfn<-log(dataG$fn+1)

#Function to get SEs for trait data. 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Create figure A : variation in elaisome mass among populations. 
ELA<-summarySE(dataG, measurevar="elave", groupvars=c("site"), na.rm=TRUE)

ela<-ggplot(ELA,aes(y=elave,x=site))+geom_errorbar(aes(ymin=elave-se,ymax=elave+se),width=1)+geom_point(size=3)
ela<-ela+labs(y= "Elaiosome \nmass (mg)")
ela<-ela + ggtitle('a') +theme(plot.title=element_text(hjust=-0.01))
ela<-ela+labs(x= "Population")
ela<-ela +theme(axis.title.x=element_blank())
ela<-ela +theme(axis.text.x=element_blank())
ela

#Create figure B : variation in EFN production among populations. 
EFN<-summarySE(dataG, measurevar="efn_nec", groupvars=c("site"), na.rm=TRUE)

efn<-ggplot(EFN,aes(y=efn_nec,x=site))+geom_errorbar(aes(ymin=efn_nec-se,ymax=efn_nec+se),width=1)+geom_point(size=3)
efn<-efn+labs(y= "Extrafloral nectar production \n(mg sucrose/nectary)")
efn<-efn + ggtitle('b') +theme(plot.title=element_text(hjust=-0.01))
efn<-efn+labs(x= "Population")
efn<-efn +theme(axis.title.x=element_blank())
efn<-efn +theme(axis.text.x=element_blank())
efn

#Create figure C : variation in FN production among populations. 
FN<-summarySE(dataG, measurevar="fn", groupvars=c("site"), na.rm=TRUE)

fn<-ggplot(FN,aes(y=fn,x=site))+geom_errorbar(aes(ymin=fn-se,ymax=fn+se),width=1)+geom_point(size=3)
fn<-fn+labs(y= "Floral nectar \nproduction (mg sucrose)")
fn<-fn+scale_x_discrete(labels=c("Brown's Town", "Cave","Galina","Haining","Hector's River","Hopewell","Mosquito Cove","Murdock","Nonsuch","St. Ann's Bay","St. Elizabeth","Salt Gut","Spicy Hill","Salt Marsh", "Spur Tree", "Titchfield","Whithorn"))
fn<-fn + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))
fn<-fn + ggtitle('c') +theme(plot.title=element_text(hjust=-0.01))
fn<-fn+theme(axis.title.x=element_blank())
fn

#Create figure D : variation in stigma-anther separation among populations. 
HRK<-summarySE(dataG, measurevar="herk", groupvars=c("site"), na.rm=TRUE)

hrk<-ggplot(HRK,aes(y=herk,x=site))+geom_errorbar(aes(ymin=herk-se,ymax=herk+se),width=1)+geom_point(size=3)
hrk<-hrk+labs(y= "Stigma-anther \nseparation (mm)")
hrk<-hrk+scale_x_discrete(labels=c("Brown's Town", "Cave","Galina","Haining","Hector's River","Hopewell","Mosquito Cove","Murdock","Nonsuch","St. Ann's Bay","St. Elizabeth","Salt Gut","Spicy Hill","Salt Marsh", "Spur Tree", "Titchfield","Whithorn"))
hrk<-hrk + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))
hrk<-hrk + ggtitle('d') +theme(plot.title=element_text(hjust=-0.01))
hrk<-hrk+theme(axis.title.x=element_blank())
hrk

plot_pop<-grid.arrange(ela,efn,fn,hrk,ncol=2,nrow=2)

plot_grid(ela,efn,fn,hrk, ncol=2,nrow=2, axis = "tblr", align="h", rel_heights = c(1, 1))
plot_grid(ela,efn,fn,hrk, ncol=2,nrow=2, axis = "tblr", align="hv", rel_heights = c(1, 1))
plot_grid(ela,efn,fn,hrk, ncol=2,nrow=2, axis = "tblr", align="v", rel_heights = c(1, 1.3))

ela + efn + fn + hrk + plot_layout(ncol = 2, nrow = 2)
