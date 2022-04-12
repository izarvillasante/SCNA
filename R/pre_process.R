# pre_process<-function(ss,subf,purity=T,query=T,out="analysis/ChAMP/"){
#   ##Create dir ChAMP/n/. put all the sample sheet idat files here
#   ChAMP_folder<-paste0(out,subf)
#   #ChAMP_folder<-paste(getwd(),ChAMP_folder,sep="/")
#
#   dir.create(ChAMP_folder)
#   ## this directory, you should copy the file "Sample_Sheet.csv", as ChAMP will look for it.
#   message("working directory: ",ChAMP_folder)
#   file.remove(list.files(ChAMP_folder, full.names = TRUE))
#   file.copy(paste0(ss$filenames,"_Grn.idat"),ChAMP_folder)
#   file.copy(paste0(ss$filenames,"_Red.idat"),ChAMP_folder)
#   message(ss$Basename)
#
#   write.csv(ss,paste0(ChAMP_folder,"/Sample_sheet.csv"))
#   data(hm450.manifest.hg19 )
#   #Run ChAMP with method for QC minfi
#   myLoad <- ChAMP::champ.load(ChAMP_folder, method="minfi")
#   #read.metharray.exp(base = ChAMP_folder)->a
#   ss<-myLoad$pd
#   #save(myLoad,file = paste0(ChAMP_folder,"TCGA_New.RData"))#Give it the name you want.
#   message("Champ does load data")
#   if(purity==T){
#     purity<-purify(myLoad=myLoad)
#     message("Rfpurifies")
#     ss$Purity_Impute_RFPurity.Absolute. <- purity
#     write.table(ss,paste(ChAMP_folder,"Sample_Sheet.txt",sep="/"), col.names = T, row.names = F, quote = F, sep="\t")
#     message(paste("ss is saved ",subf))
#   }
#   if(query==T){
#     query <- queryfy(myLoad$rgSet,ss=ss,ChAMP_folder=ChAMP_folder)
#     saveRDS(query,paste0(ChAMP_folder,"/intensities.rds"),compress = FALSE)
#     message("Pre-processing Completed successfully!")
#     return(query)
#   }
# }
