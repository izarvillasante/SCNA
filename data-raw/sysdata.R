source("./data-raw/anno.R")
control<-readRDS("/data/SCNA/control.rds")
control<-control[1:25,]
<<<<<<< HEAD
library("RFpurify")
RFpurify_ABSOLUTE<-RFpurify::RFpurify_ABSOLUTE
library(ChAMP)
data(hm450.manifest.hg19)

usethis::use_data(anno,control,RFpurify_ABSOLUTE,hm450.manifest.hg19, overwrite = TRUE,internal = TRUE)
=======
usethis::use_data(anno,control, overwrite = TRUE,internal = TRUE)
>>>>>>> package
