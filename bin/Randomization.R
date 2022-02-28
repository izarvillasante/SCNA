
library(doParallel)
library(parallel)
library(dplyr)
library(data.table)
source("Data/functions.R")

all_samples<-ss<-fread("Data/randomization_sample_sheet.txt")
CancerGenes<-readxl::read_xlsx("Supplementary_table_S3.xlsx")%>%as.data.frame()
dd <- toGRanges(CancerGenes)

for (cnstate in names(ss)[5:10]){#c("ASCAT_Gain","ASCAT_AMP","ASCAT_AMP10")){
  cl <- makeCluster(detectCores()-1,outfile="")
  registerDoParallel(cl)
  
  all_scnas <-
    foreach(sample=all_samples$Sample_Name,
            .combine='rbind', 
            .inorder=FALSE,
            .errorhandling = "pass",
            .packages =c("data.table","dplyr","regioneR")
    ) %dopar% {
      scna(cn = cnstate, ID=sample,samp=all_samples[Sample_Name==sample,][1,],dd=dd)
      
      #data.frame(x=obj$x, a=a, b=b, err=obj$err)
    }
  stopCluster(cl)
  saveRDS(all_scnas,paste0("all_scnas_",cnstate,".rds"))

  sample_pools<-lapply(1:10000,function(x) sample(unique(all_scnas$ID),size=422))

  # cl <- makeCluster(detectCores()-1,outfile="")

  Kcs<-
    foreach(samples=sample_pools,
            .combine='rbind', 
            .inorder=TRUE,
            .errorhandling = "pass"
    )%dopar%{
	    
              x<-all_scnas[all_scnas$ID %in% samples,"X"]-all_scnas[all_scnas$ID %in% samples,"Int"]
              y<-all_scnas[all_scnas$ID %in% samples,"Var"]
              fit <-lm(y~0+x)
              coeff <- fit$coefficients
              K=1/coeff
              
    }
   stopCluster(cl)
  
  saveRDS(Kcs,paste0("Kcs_",n,".rds"))
  message("saved")
}

K=list(Amp10=3.666485,Amp=1.495666,Gains=1.096756,HetLoss=-1.887569,HomDel=-5.124848)
# 3.6737, 1.7056, 1.118, -1.8844 and -4.908 for c = {Amp10, Amp, Gain,
#   HetLoss or HomDel}, r
mean_list<-list()
Kcs_ASCAT_AMP10 <- readRDS("~/Documents/Projects/20211210CUP/Kcs_ASCAT_AMP10.rds")

library(ggplot2)
AMP10<-as.data.frame(Kcs_ASCAT_AMP10)
names(AMP10)<-"Kc_AMP10"
q<-quantile(AMP10$Kc_AMP10,probs=c(0.25,0.5,0.75))
mean_list[["Kc_AMP10"]]<-q[2]

o<-quantile(AMP10$Kc_AMP10,probs=c(0.05,0.95))
ggplot(AMP10, aes(x = Kc_AMP10)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), bins = 100) +
  scale_y_continuous(labels = scales::percent)+
  geom_vline(xintercept = K[["Amp10"]],colour="blue")+
  geom_vline(xintercept = q,colour="green",linetype="dotted")+
  geom_vline(xintercept = o,colour="red",linetype="dotted")+
  ggtitle(paste0("Kc's randomization for AMP10"))

Kcs_ASCAT_AMP <- readRDS("~/Documents/Projects/20211210CUP/Kcs_ASCAT_AMP.rds")

AMP<-as.data.frame(Kcs_ASCAT_AMP)
names(AMP)<-"Kc_AMP"
q<-quantile(AMP$Kc_AMP,probs=c(0.25,0.5,0.75))
mean_list[["Kc_AMP"]]<-q[2]

o<-quantile(AMP$Kc_AMP,probs=c(0.05,0.95))
ggplot(AMP, aes(x = Kc_AMP)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), bins = 100) +
  scale_y_continuous(labels = scales::percent)+
  geom_vline(xintercept = K$Amp,colour="blue")+
  geom_vline(xintercept = q,colour="green",linetype="dotted")+
  geom_vline(xintercept = o,colour="red",linetype="dotted")+
  ggtitle(paste0("Kc's randomization for AMP"))
Kcs_ASCAT_Gain <- readRDS("~/Documents/Projects/20211210CUP/Kcs_ASCAT_Gain.rds")

Gain<-as.data.frame(Kcs_ASCAT_Gain)
names(Gain)<-"Kc_Gain"
q<-quantile(Gain$Kc_Gain,probs=c(0.25,0.5,0.75))
mean_list[["Kc_Gain"]]<-q[2]

o<-quantile(Gain$Kc_Gain,probs=c(0.05,0.95))
ggplot(Gain, aes(x = Kc_Gain)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), bins = 100) +
  scale_y_continuous(labels = scales::percent)+
  geom_vline(xintercept = K$Gains,colour="blue")+
  geom_vline(xintercept = q,colour="green",linetype="dotted")+
  geom_vline(xintercept = o,colour="red",linetype="dotted")+
  ggtitle(paste0("Kc's randomization for Gain"))

Kcs_ASCAT_diploid <- readRDS("~/Documents/Projects/20211210CUP/Kcs_ASCAT_diploid.rds")

diploid<-as.data.frame(Kcs_ASCAT_diploid)
names(diploid)<-"Kc_diploid"
q<-quantile(diploid$Kc_diploid,probs=c(0.25,0.5,0.75))
mean_list[["Kc_diploid"]]<-q[2]

o<-quantile(diploid$Kc_diploid,probs=c(0.05,0.95))
ggplot(diploid, aes(x = Kc_diploid)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), bins = 100) +
  scale_y_continuous(labels = scales::percent)+
  geom_vline(xintercept = 0,colour="blue")+
  geom_vline(xintercept = q,colour="green",linetype="dotted")+
  geom_vline(xintercept = o,colour="red",linetype="dotted")+
  ggtitle(paste0("Kc's randomization for diploid"))

Kcs_ASCAT_HetLoss <- readRDS("~/Documents/Projects/20211210CUP/Kcs_ASCAT_HetLoss.rds")

HetLoss<-as.data.frame(Kcs_ASCAT_HetLoss)
names(HetLoss)<-"Kc_HetLoss"
q<-quantile(HetLoss$Kc_HetLoss,probs=c(0.25,0.5,0.75))
mean_list[["Kc_HetLoss"]]<-q[2]

o<-quantile(HetLoss$Kc_HetLoss,probs=c(0.05,0.95))
ggplot(HetLoss, aes(x = Kc_HetLoss)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), bins = 100) +
  scale_y_continuous(labels = scales::percent)+
  geom_vline(xintercept = K$HetLoss,colour="blue")+
  geom_vline(xintercept = q,colour="green",linetype="dotted")+
  geom_vline(xintercept = o,colour="red",linetype="dotted")+
  ggtitle(paste0("Kc's randomization for HetLoss"))

Kcs_ASCAT_Homdel <- readRDS("~/Documents/Projects/20211210CUP/Kcs_ASCAT_Homdel.rds")

Homdel<-as.data.frame(Kcs_ASCAT_Homdel)
names(Homdel)<-"Kc_Homdel"
q<-quantile(Homdel$Kc_Homdel,probs=c(0.25,0.5,0.75))
mean_list[["Kc_Homdel"]]<-q[2]

o<-quantile(Homdel$Kc_Homdel,probs=c(0.05,0.95))
ggplot(Homdel, aes(x = Kc_Homdel)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), bins = 100) +
  scale_y_continuous(labels = scales::percent)+
  geom_vline(xintercept = K$HomDel,colour="blue")+
  geom_vline(xintercept = q,colour="green",linetype="dotted")+
  geom_vline(xintercept = o,colour="red",linetype="dotted")+
  ggtitle(paste0("Kc's randomization for Homdel"))


K_random<-mean_list
K_random$Kc_diploid<-NULL
names(K_random)<-names(K)
saveRDS(K_random,"K_random.rds")
