## Sample sheet for 442 TCGA TrainingSet arrays.
library(data.table)
setDTthreads(0L)
cdir<-"/data/SCNA/data-raw/TrainingSet_Arrays"
ss_train<-TrainingSet_Sample_sheet <- fread("Data/Additional.File.2_TableS1.csv")
setDT(TrainingSet_Sample_sheet)
setkey(TrainingSet_Sample_sheet,"Sample_Name")

files<-list.files(cdir,recursive = T)
idx<-1:(length(files)/2)*2 #Pairs, there are other better ways though

Sample_Name <- files%>%
  strsplit(.,"/")%>%
  sapply(.,"[",2)%>%
  unique

Bname<-files%>%
  strsplit(.,"/")%>%
  sapply(.,"[",3)%>%
  strtrim(.,17)%>%
  unique
fnames<-paste(cdir,dirname(files)[idx],Bname,sep="/")
#
# ta_files <- list.files("/data/20211210CUP/TrainingSet_Arrays/",pattern = "\\.idat$",recursive = T)
# ta_ids <-sapply(ta_files,function(x)strsplit(x = x,"/")[[1]][[2]])
#
# for (ID in ss_train$Sample_Name){
#   if (ID %in% ta_ids){
#     file <- names(ta_ids[ta_ids == ID])[1]
#     name <- substr(file,1,nchar(file)-9)
#     filename <- paste0(cdir,name)
#   }else{
#     fname <- ss_train$filenames[ss_train$Sample_Name == ID]
#     barcode <-basename(fname)
#     if( barcode %in% dori_ids){
#       file <- dori_files[dori_ids == barcode][1]
#       name <- substr(file,1,nchar(file)-9)
#       filename <- paste0("/home/idevillasante/shares/BACKUP/#recycle/raw/Arrays/TheCancerGenomeAtlas/",name)
#     }else{
#       message(paste0("missing file: ", ID))
#       filename<-"NA"
#       download <- c(download,ID)
#     }
#   }
#   print(filename)
#   setkey(ss_train,"Sample_Name")
#   ss_train[ID,filenames:=filename]
# }
setnames(TrainingSet_Sample_sheet,
         c("ASCAT_HomDel","ASCAT_HetLoss",
           "ASCAT_Gain","ASCAT_AMP","ASCAT_AMP10"),
         c("HomDel","HetLoss","Gains","Amp","Amp10")
         )

TrainingSet_Sample_sheet<-TrainingSet_Sample_sheet[Sample_Name,filenames:=(fnames)]
TrainingSet_Sample_sheet[Sample_Name,Basename:=Bname]
TrainingSet_Sample_sheet$Purity_RFPurify=NA

TrainingSet_Sample_sheet$Sample_Plate =NA
TrainingSet_Sample_sheet$Sample_Group = "Cancer"
TrainingSet_Sample_sheet$Pool_ID=NA
TrainingSet_Sample_sheet$Project=NA
TrainingSet_Sample_sheet$Sample_Well =NA
usethis::use_data(TrainingSet_Sample_sheet, overwrite = TRUE)
