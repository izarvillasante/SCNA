## code to download `CN_calls` folder and prepare `all_CN_ASCAT` dataset:
library(data.table)
setDTthreads(0L)
#cancers used in paper:
cancer_types<- c("BRCA", "BLCA", "CESC", "HNSC", "LUSC", "PRAD", "STAD", "UCEC", "KIRC", "LGG",  "LIHC", "LUAD", "PCPG", "GBM",  "PAAD", "SARC", "KICH", "SKCM")

page <- "https://github.com/riazn/biallelic_hr/raw/master/Supplementary_Files/CopyNumberData.tgz"
download.file(url=page,destfile = "./data-raw/CopyNumberData.tgz")
tmp="./data-raw/CopyNumberData"
dest="./data-raw/CN_calls"
untar("./data-raw/CopyNumberData.tgz",exdir=tmp)
unlink(dest,recursive = T)
file.rename("./data-raw/CopyNumberData/CN_calls/","./data-raw/CN_calls")
unlink(tmp,recursive = T)

#From the 25 cancer types present on this study only the previously described 18 are chosen.
DTlist <- lapply(list.files(dest,full.names = T),
                 function(x) {
                   # Only want cancers present in paper defined above:
                   can<-strsplit(x,"_")[[1]][3]
                   if ( can %in% cancer_types){
                     CNfile=fread(x)
                     names(CNfile)[1]<-"barcodes"
                     CNfile$Project =  paste("TCGA",can,sep = "-")
                     CNfile$ASCAT.path = x
                     start <- c("Project","barcodes","ASCAT.path")
                     setcolorder(
                       x = CNfile,
                       neworder = c(start,names(CNfile)[!(names(CNfile) %in% start)])
                     )
                     return(CNfile)
                   }
                 }
)
# Transform list of data.tables into a single DT:
all_CN_ASCAT_raw <- data.table::rbindlist(DTlist)
#all_CN_ASCAT <- all_CN_ASCAT_raw[apply(all_CN_ASCAT_raw[,-c(1:3)],1,function(x)!(all(is.na(x)))),]

genenames<- setdiff(names(all_CN_ASCAT_raw),c("Project","barcodes","ASCAT.path"))
empty_rows <- apply(all_CN_ASCAT_raw[,.SD,.SDcols=genenames],1,function(x)all(is.na(x)))
all_CN_ASCAT<-all_CN_ASCAT[!(empty_rows),]
usethis::use_data(all_CN_ASCAT, overwrite = TRUE)
