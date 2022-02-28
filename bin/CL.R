################# Declare functions:
library(filesstrings)
library(preprocessCore)
library(ChAMP)
library(ChAMPdata)
library(DNAcopy)
library(data.table)

library(minfi)
library("minfiData")
## Exclude cross reactive probes
library(devtools)
#install_github("markgene/maxprobes", force=T)
library(conumee)
library(minfi)
library("minfiData")
library(regioneR)
library(filesstrings)




# genes ranges:

CancerGenes<-readxl::read_xlsx("Supplementary_table_S3.xlsx")%>%as.data.frame()
dd <- toGRanges(CancerGenes)



#################### Download GEO:
#ss_validation_CL<- read.csv("~/Documents/Projects/20211210CUP/analysis/CL/Sample_Sheet_ValidationCohort_CNV.Annotated_5copies.csv")
ss_validation_CL<-fread("analysis/ChAMP/CL_validation/Sample_Sheet.txt")
sampinfo<- readxl::read_excel("metadata/Samples_on_Array_20191213.xlsx",skip = 7)
info_ss_validation<-sampinfo[  sampinfo$barcode %in% ss_validation_CL$Basename,]
library(GEOquery)
gse <- getGEO("GSE68379",GSEMatrix=TRUE)
gse_sample_names<-str_remove_all( pData(phenoData(gse[[1]]))[,8],"-")
idx<-gse_sample_names %in% ss_validation_CL$Sample_Name
gsms<-colnames(gse[[1]]@assayData$exprs)[idx]
?getGEO()
for (GSM in gsms){
  if (!(file.exists(paste0("./analysis/CL",GSM)))){
  # gsm<-getGEO(GSM)
    e<-simpleError(paste0(GSM,"not found"))
  tryCatch( GEOquery::getGEOSuppFiles(GSM,fetch_files = T,baseDir = "./analysis/CL"),
             print(GSM),
             error=function(e)e)

  }
}
getwd()


ss_train_CL<- read.csv("~/Documents/Projects/20211210CUP/analysis/CL/Sample_Sheet_TrainingCohort_CNV_Annotated_5copies.csv")
# sampinfo<- read_excel("metadata/Samples_on_Array_20191213.xlsx",skip = 7)
info_ss_train<-sampinfo[  sampinfo$barcode %in% ss_train_CL$Basename,]
library(GEOquery)
gse <- getGEO("GSE68379",GSEMatrix=TRUE)
gse_sample_names<-str_remove_all( pData(phenoData(gse[[1]]))[,8],"-")
idx<-gse_sample_names %in% ss_train_CL$Sample_Name
gsms<-colnames(gse[[1]]@assayData$exprs)[idx]
?getGEO()
for (GSM in gsms){
  if (!(file.exists(paste0("./analysis/CL",GSM)))){
    # gsm<-getGEO(GSM)
    e<-simpleError(paste0(GSM,"not found"))
    tryCatch( GEOquery::getGEOSuppFiles(GSM,fetch_files = T,baseDir = "./analysis/CL"),
              print(GSM),
              error=function(e)e)

  }
}
getwd()

ss_validation_CL$filenames

ss_train_CL$filenames
gsms
pheno_all<-pData(phenoData(gse[[1]]))
setDT(pheno_all)
pheno_all$Sample_Name<-str_remove_all( pheno_all$source_name_ch1,"-")
#"TT and T-T problem, remove row 988. check Basename in pheno_all$supplementary_file "
pheno_all<-pheno_all[-988,]
setkey(pheno_all,"Sample_Name")
pheno_train<-pheno_all[ss_train_CL$Sample_Name,]

ss_train<-merge(pheno_train,ss_train_CL)
identical(ss_train$Sample_Name,ss_train_CL$Sample_Name)
ss_train$filenames<-paste0(getwd(),"/analysis/CL/",ss_train$geo_accession,"/",
                           substr(basename(ss_train$supplementary_file),1,28))
ss_train[,c("title","Sample_Name","geo_accession","Basename","filenames")]

setDT(ss_train_CL)
setkey(ss_train_CL,"Sample_Name")
#ss_train_CL[ss_train$Sample_Name,filenames:=ss_train$filenames]
#ss_train_CL[ss_train$Sample_Name,Basename:=ss_train$Basename]

samps_train<-data.table(
  Sample_Name=ss_train$Sample_Name,
  filenames=ss_train$filenames,
  Cancer=ss_train$`primary site:ch1`,
  Purity_Impute_RFPurity.Absolute.=ss_train$Purity_Impute_Absolute,
  ASCAT_Homdel= ss_train$HomDel,
  ASCAT_HetLoss= ss_train$HetLoss,
  ASCAT_diploid= "",
  ASCAT_Gain= ss_train$Gain,
  ASCAT_AMP= ss_train$Amp1,
  ASCAT_AMP10 = ss_train$Amp2,
  Sample_Plate =NA,
  Sample_Group ="Cell Line",
  Pool_ID=NA,
  Project=NA,
  Sample_Well =NA

)

if(!file.exists("CLtrain_intensities.fst")){
  intensities <- pre_process(samps_train,subf = "CL_train")
  intensities$probeid<-rownames(intensities)
  fst::write_fst(intensities,"CLtrain_intensities.fst")
}
pheno_validation<-pheno_all[ss_validation_CL$Sample_Name,]

ss_validation<-merge(pheno_validation,ss_validation_CL)
identical(ss_validation$Sample_Name,ss_validation_CL$Sample_Name)
# ss_validation$filenames<-paste0(getwd(),"/analysis/CL/",ss_validation$geo_accession,"/",
#                            substr(basename(ss_validation$supplementary_file),1,28))
samps_validation<-data.table(
  Sample_Name=ss_validation$Sample_Name,
  filenames=ss_validation$filenames,
  Cancer=ss_validation$`primary site:ch1`,
  Purity_Impute_RFPurity.Absolute.=ss_validation$Purity_Impute_Absolute,
  ASCAT_Homdel= ss_validation$HomDel,
  ASCAT_HetLoss= ss_validation$HetLoss,
  ASCAT_diploid= "",
  ASCAT_Gain= ss_validation$Gain,
  ASCAT_AMP= ss_validation$Amp1,
  ASCAT_AMP10 = ss_validation$Amp2,
  Sample_Plate =NA,
  Sample_Group ="Cell Line",
  Pool_ID=NA,
  Project=NA,
  Sample_Well =NA
  
)

if(!file.exists("CLvalidation_intensities.fst")){
  intensities <- pre_process(samps_validation,subf = "CL_validation",query = T)
  intensities$probeid<-rownames(intensities)
  fst::write_fst(intensities,"CLvalidation_intensities.fst")
}


# Calculate TPR/FPR with new Kc values:

# We are using validation set:

 ss_purity_val<-"analysis/ChAMP/CL_validation/Sample_Sheet.txt"
 ss_validation_purity<-fread(ss_purity_val)
 ss_validation_purity$Purity_Impute_RFPurity.Absolute.
 ss_validation$Purity_Impute_RFPurity.Absolute.<-ss_validation_purity$Purity_Impute_RFPurity.Absolute.

all_samples<-ss<-ss_validation
ss$GAINS<-apply(ss,1,function(x){
  a<-paste(x["Gain"],x["Amp1"],x["Amp2"],collapse=";",sep=";")
  b<-unlist(strsplit(a,";"))
  c<-intersect(b,dd$name)
  d<-paste(c,collapse = ";")
  return(d)
})
ss$LOSS<-apply(ss,1,function(x){
  a<-paste(x["HetLoss"],x["HomDel"],collapse=";",sep=";")
  b<-unlist(strsplit(a,";"))
  c<-intersect(b,dd$name)
  d<-paste(c,collapse = ";")
  return(d)
})


setkey(ss,"Sample_Name")
#Calculate predictions:
intlist<-list()

# CONUMEE KC:

###################### validation:

# cl <- makeCluster(4,outfile="")
# registerDoParallel(cl)
# start_time <- Sys.time()
# query_list<-foreach(k=1:NROW(all_samples),
#                     .combine="rbind",
#                     .errorhandling = "pass",
#                     .packages =c("conumee","fst","dplyr")
# ) %dopar%
#   {
# 
#     run_conumee(anno=anno, controls=controls, ss = all_samples,k=k,
#                 fst.file = "CLvalidation_intensities.fst"  )
#   }
# stopCluster(cl)
ss
int_list<-list()
ss<-ss_train
K_list<-list()
K_list$HomDel<-(-5.124848)
K_list$HetLoss<-(-1.887569)
K_list$Gains<-  1.096756
K_list$Amp2<- 3.666485
K_list$Amp1<- 1.495666
for (i in ss$Sample_Name){
  #i<-ss$Sample_Name[21]
  ID<-i
  
  CONUME_KC_GAINS <- conumee_kc(ID=ID,segfolder="analysis/CONUMEE/",samp= ss[i,],dd, K=K_list)
  message("done CONUME_KC")
  intlist[[i]]<-CONUME_KC_GAINS
  for (cnstate in names(CONUME_KC_GAINS)){
    col<-paste0("TP_CONUMEE_KC_",cnstate)
    ss[i,(col):=CONUME_KC_GAINS[[cnstate]][[1]]$TPR]
    col_fp<-paste0("FP_CONUMEE_KC_",cnstate)
    ss[i,(col_fp):=CONUME_KC_GAINS[[cnstate]][[1]]$FPR]
    CONUME_Std <-conumee_st(ID=ID,segfolder="analysis/CONUMEE/",samp= ss[i,],dd=dd) 
    ss[i,TP_CONUMEE_Std_all:=CONUME_Std$ASCAT_GAINS[[1]]$TPR ]
    ss[i,FP_CONUMEE_Std_all:=CONUME_Std$ASCAT_GAINS[[1]]$FPR ]
  }
}


saveRDS(ss,"ss_CELL_LINE_TPR.rds")

ss_plot<-select(ss,all_of(names(ss)[startsWith(names(ss), "TP")| startsWith(names(ss), "FP")]))
#ss_plot<-ss_plot%>%select(names(ss_plot)[!names(ss_plot)%like% "all"])

pdata<-data.frame(
  val=as.numeric(ss_plot[,lapply(.SD,mean,na.rm=T)]),
  test=ifelse(startsWith(names(ss_plot), "TP"),"TP","FP"),
  prog=substr(names(ss_plot),4,nchar(names(ss_plot)))
)
writexl::write_xlsx(pdata,"fig1E.xlsx")
library(ggplot2)
ggplot2::ggplot(data=pdata, aes(x=test,y=val,fill=test))+
  geom_bar(stat="identity",position = "dodge")+
  ggtitle("Fig.1E: CONUMEE vs ASCAT algorithm") +
  facet_wrap(~prog)

ggsave("Fig_1E_facet.png")

