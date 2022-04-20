#' Generate segment & log2r intensity files calculated by CONUMEE.
#' @param intensities dataset with intensities. either dataframe or path to file.
#' For big datasets '.fst' file is recomended.
#' @param anno_file anno file CNV.anno object saved as .rds. default
#' CNV_Germline_GSITIC2_BROAD_SNP6.merged.151117.hg19.CNV.txt
#' @param ctrl_file CNV data object with intensities from the controls group.
#' must be compressed as '.rds' file format, default uses 96 WB samples.
#' Default='WB'
#' @param out Parent working directory where you want to have your results.
#' @param seg.folder  Subfolder inside results where segment files are saved.
#' The segments are generated with CONUMEE CNV.segment(CNV.detail(CNV.bin(fit)))
#' default = "Segments"
#' @param log2r.folder Subfolder inside results where log2r values are saved.
#' log2r values are the fitted log intensity value against the reference set
#' returned by CNV.fit function from CONUMEE.
#' default = "log2"
#' @param Sample_Name Samples to be analysed. If is NULL all columns in input
#' file will be used. Accepts numbered index and names. Default=NULL
#' @param probeid name of column with probe ids.
#' @return Log2r intensities and segment files.
#' @export
#'
#' @examples
#'
#' #data("anno")
#' #data("controls")
#' #library(conumee)
#' library(SummarizedExperiment)
#' intensity<-readRDS("./analysis/ChAMP/intensities.rds")
#' intensity<-conumee::CNV.load(intensity)@intensity
#' run_conumee(intensities=intensity)



run_conumee<-function(intensities, anno_file=NULL, ctrl_file='WB', Sample_Name=NULL,
                      seg.folder = "Segments", log2r.folder = "log2",
                      out="analysis/CONUMEE/", probeid="probeid"){
  requireNamespace("conumee")
  if(!is.null(anno_file ))anno <- readRDS(anno_file)
  message("anno")

  if(ctrl_file != "WB" )control <- readRDS(ctrl_file)
  message("controls")
  if(is.null(intensities) ){
    stop("Provide conumee with either a valid path to file or a data.frame with intensities",call. = F)
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  intensities<- read_intensity(infile = intensities,Sample_Name=Sample_Name,probeid=probeid)
  data.table::setkey(intensities,probeid)
  message("intensities")

  # Intersect common probes
  cg1<-intersect(intensities[[probeid]],rownames(control@intensity))
  cgcommon<-intersect(cg1,names(anno@probes@ranges))
  message(head(cgcommon))
  # Subset
  intensities<-intensities[cgcommon,]
  control@intensity <-control@intensity[cgcommon,]
  anno@probes<-anno@probes[names(anno@probes)%in% cgcommon,]
  # conumee:

  idx<-intensities[[probeid]]
  intensities[,(probeid):=NULL]
  # Cluster:
  ncores<-parallel::detectCores()-2
  if (!is.na(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))))ncores=as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  cl<- parallel::makePSOCKcluster(ncores, outfile="")
  doParallel::registerDoParallel(cl)

  res<-foreach::foreach(i=1:NCOL(intensities),#isplitIndices(1400,chunks=ncores),
                        .combine='c',
                        .multicombine = F,
                        .inorder=F,
                        .errorhandling = "pass"
  )%dopar%{
    Sample_Name <- names(intensities)[i]
    my.datacnv<-conumee::CNV.load(intensities[,..i])
    rownames(my.datacnv@intensity)<-idx
    fit  <- conumee::CNV.fit(my.datacnv, control, anno)
    #Log2 <- list(log2ratio=fit@fit$ratio, Purity=ss$Purity_Impute_RFPurify.Absolute.[ss$Sample_Name%in%Sample_Name])
    fit2 <- conumee::CNV.segment(conumee::CNV.detail(conumee::CNV.bin(fit)))
    write.table(conumee::CNV.write(fit2, what = "segments"),sprintf("%sSegments_%s.txt",out,Sample_Name),col.names = T, row.names = F, quote = F, sep = "\t")
    message("Sgements file" ," saved")
    log2ratio<-as.data.frame(fit@fit$ratio)
    names(log2ratio) <- Sample_Name
    #log2ratio$purity <- ss[ss$Sample_Name== Sample_Name,Purity_Impute_RFPurity.Absolute.]

    log2_file<-paste(out,log2r.folder,"/",Sample_Name,"_log2.txt", sep="")
    dir.create(dirname(log2_file),recursive = TRUE,showWarnings=FALSE)
    message(log2_file)
    write.table(log2ratio,log2_file)
    message(log2_file ," saved")
    return(log2ratio)
  }
  parallel::stopCluster(cl)
return(res)
  }
#
#   # rownames(intensities)<-intensities$probeid
#   # names(intensities)[2]<-"intensity"
#
#   my.datacnv<-conumee::CNV.load(input= intensities[,2])
#   rownames(my.datacnv@intensity)<-intensities$probeid
#   fit  <- conumee::CNV.fit(my.datacnv, control, anno)
#   #Log2 <- list(log2ratio=fit@fit$ratio, Purity=ss$Purity_Impute_RFPurify.Absolute.[ss$Sample_Name%in%Sample_Name])
#   fit2 <- conumee::CNV.segment(conumee::CNV.detail(conumee::CNV.bin(fit)))
#   write.table(conumee::CNV.write(fit2, what = "segments"),sprintf("%sSegments_%s.txt",out,Sample_Name),col.names = T, row.names = F, quote = F, sep = "\t")
#   message("Sgements file" ," saved")
#   log2ratio<-as.data.frame(fit@fit$ratio)
#   names(log2ratio) <- Sample_Name
#   #log2ratio$purity <- ss[ss$Sample_Name== Sample_Name,Purity_Impute_RFPurity.Absolute.]
#   write.table(log2ratio,log2_file)
#   message(log2_file ," saved")
#   return(log2ratio)
#   #}else return(NULL)
# }

#' @export
read_intensity<-function(infile,Sample_Name=NULL,probeid="probeid"){
  requireNamespace("data.table")
if(!is.character(probeid)&length(probeid)==1){
  stop("probeid must be the name of a single column. ")
}
tryCatch(
  if(is.character(infile)){
    if(endsWith(infile,".rds")){
      my.data <- readRDS(infile)
      data.table::setDT(my.data)
      n<-colnames(my.data)
      if(is.null(Sample_Name))Sample_Name<-setdiff(n,probeid)
      n_ok<-intersect(Sample_Name,n)
      n_missing<- setdiff(Sample_Name,n_ok)
      if(length(n_missing) > 0) warning(n_missing, " not found in file")
      my.data <- my.data[,.SD,.SDcols=c(probeid,n_ok)]

    }else if (endsWith(infile,".fst")){
      meta<-fst::metadata_fst(infile)
      n<-meta$columnNames
      if(is.null(Sample_Name))Sample_Name<-setdiff(n,probeid)
      n_ok<-intersect(Sample_Name,n)
      n_missing<- setdiff(Sample_Name,n_ok)
      if(length(n_missing) > 0) warning(n_missing, " not found in file")
      my.data <- fst::read_fst(infile,columns=c(probeid,n_ok),as.data.table = T)
    }else{
      my.data <- data.table::fread(infile)
      if(!probeid %in% names(my.data)) data.table::setnames(my.data,"V1",probeid,skip_absent = TRUE)
      n<-colnames(my.data)
      if(is.null(Sample_Name))Sample_Name<-setdiff(n,probeid)
      n_ok<-intersect(Sample_Name,n)
      n_missing<- setdiff(Sample_Name,n_ok)
      if(length(n_missing) > 0) warning(n_missing, " not found in file")
      my.data <- my.data[,.SD,.SDcols=c(probeid,n_ok)]
    }
  }else{
    message("dataset")
      rn<-rownames(infile)
      my.data<-data.table::as.data.table(infile)
      if(!probeid %in% names(my.data)) my.data[,(probeid):=rn]
      n<-colnames(my.data)
      if(is.null(Sample_Name))Sample_Name<-setdiff(n,probeid)
      n_ok<-intersect(Sample_Name,n)
      n_missing<- setdiff(Sample_Name,n_ok)
      if(length(n_missing) > 0) warning(n_missing, " not found in file")
      my.data <- my.data[,.SD,.SDcols=c(probeid,n_ok)]
      }, error = function(e) "invalid file"
  )
  return(my.data)

}

