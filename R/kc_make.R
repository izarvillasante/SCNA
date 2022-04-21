#' Set of functions to obtain Kc from reference copy number state.
#' @param ss Sample sheet containing purities, cnstate and Sample_Name ids.
#' @param fname identifier for your Kc reference data used in the output names.
#' @param cncols column names containing reference set of genes for each cn state
#' @inheritParams run_conumee
#' @examples
#' @export
#' ss<-data.table::fread("analysis/ChAMP/Sample_Sheet.txt")
#' make_Kc(ss=ss,fname="test",cncols=c("Amp","Amp10","Gains"))

make_Kc<-function(ss,fname,cncols=c("Amp10","Amp","Gains","Hetloss","HomDel"),
                  conumee.folder="analysis/CONUMEE",seg.folder = "Segments",
                  log2r.folder = "log2r"){
  K_list<-list()
  segfile_folder <- paste0(conumee.folder,"/",seg.folder)
  std_segfile <- rlang::expr(paste0(folder,ID,'_Segments.txt'))

  ss$segfile <- check_input_files(Sample_Name = ss$Sample_Name,folder = segfile_folder, std_file = std_segfile)

  # segfile_folder <- paste0(conumee.folder,"/",seg.folder)
  # segfile_list<-list.files(segfile_folder)
  # data.table::setDT(ss)
  # setkey(ss,"Sample_Name")
  # sapply(ss$Sample_Name,function(ID){
  #   match=FALSE
  #   match<-sapply(segfile_list, function(file){
  #     #print(match)
  #     if(like(file,ID)){
  #     match=TRUE
  #     ss[ID,segfile:=file]
  #     return(match)}
  #     return(match)})
  #   #print(sum(match))
  #   if(sum(match) == 0) warning("No segment file found for sample: ",ID)
  #   if(sum(match) > 1){
  #     std_segfile <- paste0(conumee.folder,"/",seg.folder,"/",ID,"_Segments.txt")
  #     match_segfiles <- segfile_list[match]
  #     if(std_segfile %in% match_segfiles){ok<-std_segfile}else{
  #       ok<-segfile_list[which.max(match)]}
  #     overwrite<-setdiff(match_segfiles,ok)
  #     warning("multiple files match ",ID, ". \n ", paste(overwrite,collapse = ", ")," ignored.")
  #     return(ok)
  #   }
  #   match_files <- segfile_list[match]
  #   return(match_files)
  # })

  log2rfile_folder <- paste0(conumee.folder,"/",log2r.folder)
  std_log2rfile <- rlang::expr(paste0(folder,ID,'_log2r.txt'))
  ss$log2file <- check_input_files(Sample_Name = ss$Sample_Name,folder = log2rfile_folder, std_file = std_log2rfile)
  # log2rfile_list<-list.files(log2rfile_folder)
  # data.table::setDT(ss)
  # setkey(ss,"Sample_Name")
  # a<-sapply(ss$Sample_Name,function(ID){
  #   match=FALSE
  #   match<-sapply(log2rfile_list, function(file){
  #     #print(match)
  #     if(like(file,ID)){
  #       match=TRUE
  #       ss[ID,log2rfile:=file]
  #       return(match)}
  #     return(match)})
  #   #print(sum(match))
  #   if(sum(match) == 0) warning("No segment file found for sample: ",ID)
  #   if(sum(match) > 1){
  #     std_log2rfile <- paste0(conumee.folder,"/",log2r.folder,"/",ID,"_log2r.txt")
  #     match_files <- log2rfile_list[match]
  #     if(std_log2rfile %in% match_files){ok<-std_log2rfile}else{
  #       ok<-log2rfile_list[which.max(match)]}
  #     overwrite<-setdiff(match_files,ok)
  #     warning("multiple files match ",ID, ". \n ", paste(overwrite,collapse = ", ")," ignored.")
  #     return(ok)
  #   }
  #   match_files <- log2rfile_list[match]
  #   return(match_files)
  # })

  ncores<-parallel::detectCores()-2
  if (!is.na(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))){
    ncores=as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))}
  cl<- parallel::makePSOCKcluster(ncores, outfile="")
  doParallel::registerDoParallel(cl)
  for (cnstate in cncols){
    Kc_file<-paste0(fname,"_",cnstate,".rds")
    if(file.exists(Kc_file)){K<-readRDS(Kc_file)}else{
      scnas_file<-paste0("scnas_",fname,"_",cnstate,".rds")
      #if(file.exists(scnas_file)){readRDS(scnas_file)}
      all_scnas <-
        foreach::foreach(sname=ss$Sample_Name,
                .combine='rbind',
                .inorder=FALSE,
                .errorhandling = "pass",
                .packages =c("data.table","regioneR","SCNA")
        ) %dopar% {
          segfile <- ss[sname,segfile]#paste0(conumee.folder,"/",seg.folder,"/",sname,"_Segments.txt")
          log2file <- ss[sname,log2file]#paste0(conumee.folder,"/",log2r.folder,"/",sname,"_log2r.txt")

          SCNA::get_scna(cn = cnstate, ID=sname,ss=ss[sname,],ref_genes="paper",log2file = log2file,segfile = segfile)
        }
      saveRDS(all_scnas,scnas_file)
      all_scnas<-all_scnas[complete.cases(all_scnas),]
      fit <-lm(Var~0+I(X-Int),data=all_scnas)
      coeff <- fit$coefficients
      K=1/coeff
      saveRDS(K,Kc_file)
    }
    K_list[[cnstate]]<-K
    saveRDS(K_list,paste0(fname,"_Kc_list.rds"))
  }
  stopCluster(cl)
  Ks<-as.data.frame(t(data.frame(cn=fnames(K_list),K=unlist(unname(K_list)))))
  writexl::write_xlsx(Ks, paste0(fname,"_Ks.xlsx"))
}




#' @rdname make_Kc
#' @param cn Sample sheet column name containing reference gene set for that cn.
#' @param ref_genes Granges object with subset of genes to use for cnv analysis.
#' default = "all"
#' @param segfile  file containing segment data from conumee.
#' The segments are generated with CONUMEE CNV.segment(CNV.detail(CNV.bin(fit)))
#' @param log2file file with log2r values, which are the fitted log intensity
#' value against the reference set returned by CNV.fit function from CONUMEE.
#' @param ID name of sample. Must be in Sample_Name column.
#' @return dataframe with input data for linear model. Int, mean and Var .
#' @export
#' @examples
#' ss<-data.table::fread("analysis/ChAMP/Sample_Sheet.txt")
#' res<-scna(ID="TCGA-05-4405-01A-21D-1856-05",ss=ss,ref_genes="paper",cn="Amp",
#' log2file="analysis/CONUMEE/log2/TCGA-05-4405-01A-21D-1856-05_log2.txt",
#' segfile="analysis/CONUMEE/Segments_TCGA-05-4405-01A-21D-1856-05.txt")
#' res
get_scna<-function(ID,log2file,segfile,ss,ref_genes="all",cn){
  bp_stopifnot = getFromNamespace("stopifnot", "backports")

  # Check params:
  bp_stopifnot("log2file must contain a valid path "             = is.character(log2file) & file.exists(log2file))
  bp_stopifnot("log2file must contain path to a single file"     = length(log2file)==1)
  bp_stopifnot("segfile must contain a valid path "              = is.character(segfile) & file.exists(segfile))
  bp_stopifnot("segfile must contain path to a single file"      = length(segfile)==1)
  bp_stopifnot("cn must be a valid column in sample_sheet "      = cn %in% names(ss))

  bp_stopifnot("sample_sheet must contain Sample_Name column"    = "Sample_Name" %in% names(ss))
  #bp_stopifnot("purities must be given in "    = any( sapply(names(ss),function(x)startsWith(x,"Purity_Impute_RFPurify") )))
  bp_stopifnot("ID must be of length 1"                  = length(ID)==1)
  bp_stopifnot("ID not found in Sample_Name column"      = ID %in% ss$Sample_Name)


  log2_message <- paste0("please provide a valid log2file for ",ID," you can use run_conumee.")
  seg_message <- paste0("please provide a valid segments file for ",ID," you can use run_conumee.")

  # prepare data:
  data.table::setDT(ss)
  ss<-ss[Sample_Name %in% ID,]
  cn_genes<-strsplit(as.character(with(ss,get(cn))), split = ";")[[1]]
  cn1<-sprintf("Invalid genes input for sample: %s",ID)
  bp_stopifnot( cn1 = is.character(cn_genes))
  if(length(cn_genes)<1){
    df<-data.frame(ID=NULL,Int=NULL,X=NULL,Var=NULL)
    warning("There are no ",cn," genes for sample: ",ID)
    return(df)
    }

  # rare_genes <- setdiff(cn_genes , AllGenes$name)
  # if(length(rare_genes)>=1)warning(paste(rare_genes,sep = ", ",collapse = ", "), " genes may fall in sex chromosomes or not present in ref genome: 'Hsapiens.UCSC.hg19.knownGene' ")
  # Log2 & purity:
  suppressWarnings(Log2<-data.table::fread(log2file,col.names = c("probeid","log2r")))
  baseline <- mean(Log2$log2r, na.rm=T)
  purity<-names(ss)[names(ss) %ilike% "impute" & names(ss) %ilike% "purity" ]
  p<-with(ss,get(purity))
  p<-ss[,..purity]
  Var <- p*sd(Log2$log2r,na.rm=T)

  # segments:
  seg<-read.delim(segfile,header=T,sep="\t")
  int<-get_int(seg,ref_genes=ref_genes,cn_genes=cn_genes)
  X=mean(int$log2r,na.rm=T)
  df<-data.frame(
    ID=ID,
    Int=baseline,
    X=X,
    Var=Var)
  return(df)
}

#' @rdname make_Kc
get_int<-function(seg,ref_genes="all",cn_genes){
  if(ref_genes=="all"){interest_geneset <- AllGenes}else if(ref_genes == "paper"){
    interest_geneset <- CancerGenes}else{interest_geneset<-read.table(ref_genes)}
  interest_genes<-interest_geneset$name
  genes <- intersect(cn_genes,interest_genes)
  if(length(genes)<1){
    df<-data.frame(ID=NULL,Int=NULL,X=NULL,Var=NULL)
    warning("No genes of interest present in cnstate: ",cn)
    return(df)
  }
  segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean)
  seggr <- regioneR::toGRanges(segb)
  grset <- regioneR::toGRanges(interest_geneset)
  int <- suppressWarnings(SummarizedExperiment::findOverlaps(seggr,grset[grset$name %in% genes,]))
  seggr.matched <- seggr[S4Vectors::queryHits(int)];
  S4Vectors::mcols(seggr.matched) <- cbind.data.frame(
    S4Vectors::mcols(seggr.matched),
    S4Vectors::mcols(grset[grset$name%in%genes][S4Vectors::subjectHits(int)]));
  int <- seggr.matched
  message("int finish")
  return(int)
}
#rlang::expr(paste0(conumee.folder,'/',seg.folder,'/',sname,'_Segments.txt'))
check_input_files<-function(Sample_Name,folder,std_file="standard_format"){

  file_list<-list.files(folder)
  path<-sapply(Sample_Name,function(ID){
    match=FALSE
    match<-sapply(file_list, function(file) ifelse(data.table::like(file,ID),TRUE,FALSE))
    if(sum(match) == 0) warning("No input file found for sample: ",ID)
    if(sum(match) > 1){
      #std_file <- paste0(conumee.folder,"/",.folder,"/",ID,"_Segments.txt")
      match_files <- file_list[match]
      std_file<-eval(std_file)
      if(std_file %in% match_files){ok<-std_file}else{
        ok<-file_list[which.max(match)]}
      overwrite<-setdiff(match_files,ok)
      warning("multiple files match ",ID, ". \n ", paste(overwrite,collapse = ", ")," ignored.")
      return(ok)
    }# otherwise match == 1
    match_files <- file_list[match]
    return(match_files)
  })
  return(path)
}
a<-"paste0(conumee.folder,'/',seg.folder,'/',sname,'_Segments.txt')"
test<-function(folder="myfolder/",std_file){
  std_file<-unquote(std_file);print(std_file)
  }

