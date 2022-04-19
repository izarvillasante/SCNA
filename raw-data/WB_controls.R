## Download Whole Blood samples.

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GEOquery")
#install.packages("data.table")
library(data.table)
setDTthreads(0L)
library(GEOquery)
library(filesstrings)
## Whole Blood:
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73103


cdir="./data-raw/WholeBlood_Controls"
if(!dir.exists(cdir))create_dir(cdir)
cdirf=paste(cdir,"Females",sep = "/")
if(!dir.exists(cdirf))create_dir(cdirf)
cdirm=paste(cdir,"Males",sep = "/")
if(!dir.exists(cdirm))create_dir(cdirm)

#Fetch all 48 Females Age >= 18


count=0
for(k in 1886364:1886803){
  gsm <- getGEO(sprintf("GSM%i",k));
  Sex <- gsub(" ","",strsplit(Meta(gsm)[["characteristics_ch1"]][c(1,3)][1],split=":")[[1]][2]);
  Age <- gsub(" ","",strsplit(Meta(gsm)[["characteristics_ch1"]][c(1,3)][2],split=":")[[1]][2]);
  if(Sex=="Female" & Age >=18){GEOquery::getGEOSuppFiles(sprintf("GSM%i",k),fetch_files = T,baseDir = cdirf);count=count+1;
  print(count)}
}

#Fetch  48 Males Age >= 18
dir_females <- gsub("./","",list.dirs(cdirf)[-1])
GSM_all <- paste0("GSM",c(1886364:1886588))
GSM_males <- GSM_all[!GSM_all%in%dir_females]


count=0
for(GSM in GSM_males){
  count=count+1;
  print(count);
  if(count==48)break;
  gsm <- getGEO(GSM);
  Sex <- gsub(" ","",strsplit(Meta(gsm)[["characteristics_ch1"]][c(1,3)][1],split=":")[[1]][2]);
  Age <- gsub(" ","",strsplit(Meta(gsm)[["characteristics_ch1"]][c(1,3)][2],split=":")[[1]][2]);
  if(Sex=="Male" & Age >=18)GEOquery::getGEOSuppFiles(GSM,fetch_files = T,baseDir = cdirm)
}

# Dump all  Male/ and Female/ dirs .idat files in WholeBlood_Controls.
setwd("../../")

#usethis::use_data(TrainingSet_Sample_sheet, overwrite = TRUE)
