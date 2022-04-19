source("./data-raw/anno.R")
control<-readRDS("/data/SCNA/control.rds")
control<-control[1:25,]
usethis::use_data(anno,control, overwrite = TRUE,internal = TRUE)
