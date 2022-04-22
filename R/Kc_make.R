#' Set of functions to obtain Kc from reference copy number state.
#' @param ss Sample sheet containing purities, cnstate and Sample_Name ids.
#' @param fname identifier for your Kc reference data used in the output names.
#' @param cncols column names containing reference set of genes for each cn state
#' @param ref_genes Granges object with subset of genes to use for cnv analysis.
#' default = "paper"
#' @inheritParams run_conumee
#' @export
#' @examples
#' data("ss")
#' Kc_make(ss=ss,fname="test",cncols=c("Amp","Amp10","Gains"))

Kc_make<-function(ss,fname,cncols=c("Amp10","Amp","Gains","Hetloss","HomDel"),
                  conumee.folder="analysis/CONUMEE",seg.folder = "Segments",
                  log2r.folder = "log2r",ref_genes="paper"){
  K_list<-list()
  ss<-data.table::setDT(ss)
  data.table::setkey(ss,"Sample_Name")
  segfile_folder <- paste0(conumee.folder,"/",seg.folder)
  std_segfile <- rlang::expr(paste0(folder,ID,'_Segments.txt'))#rlang::expr(paste0(folder,"/","Segments_",ID,'.txt'))
  ss$segfile <- check_input_files(Sample_Name = ss$Sample_Name,folder = segfile_folder, std_file = std_segfile)

  log2rfile_folder <- paste0(conumee.folder,"/",log2r.folder)
  std_log2rfile <- rlang::expr(paste0(folder,ID,'_log2r.txt'))
  ss$log2file <- check_input_files(Sample_Name = ss$Sample_Name,folder = log2rfile_folder, std_file = std_log2rfile)

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
        foreach::foreach(ID=ss$Sample_Name,
                .combine='c',
                .inorder=FALSE,
                .errorhandling = "pass",
                .export=c("get_int","scna"),
                .packages =c("data.table","regioneR")
        ) %dopar% {
          ID<-"TCGA-19-A6J4-01A-11D-A33U-05"
          segfile <- ss[ID,segfile]#paste0(conumee.folder,"/",seg.folder,"/",ID,"_Segments.txt")
          log2file <- ss[ID,log2file]#paste0(conumee.folder,"/",log2r.folder,"/",ID,"_log2r.txt")
          ss<-ss[ID,]
          scna(cn = cnstate, ID=ID,ss=ss,ref_genes=ref_genes,log2file = log2file,segfile = segfile)
        }
      saveRDS(all_scnas,scnas_file)
      all_scnas<-all_scnas[stats::complete.cases(all_scnas),]
      fit <-stats::lm(Var~0+I(X-Int),data=all_scnas)
      coeff <- fit$coefficients
      K=1/coeff
      saveRDS(K,Kc_file)
    }
    K_list[[cnstate]]<-K
    saveRDS(K_list,paste0(fname,"_Kc_list.rds"))
  }
  parallel::stopCluster(cl)
  Ks<-as.data.frame(t(data.frame(cn=names(K_list),K=unlist(unname(K_list)))))
  writexl::write_xlsx(Ks, paste0(fname,"_Ks.xlsx"))

}




#' Function to obtain mean intercept and variance from a
#' reference gene set and a particular copy number state
#' @param cn Sample sheet column name containing reference gene set for that cn.
#' @param segfile  file containing segment data from conumee.
#' The segments are generated with CONUMEE CNV.segment(CNV.detail(CNV.bin(fit)))
#' @param log2file file with log2r values, which are the fitted log intensity
#' value against the reference set returned by CNV.fit function from CONUMEE.
#' @param ID name of sample. Must be in Sample_Name column.
#' @inheritParams Kc_make
#' @return dataframe with input data for linear model. Int, mean and Var .
#' @export
#' @examples
#' ss<-data.table::fread("analysis/ChAMP/Sample_Sheet.txt")
#' res<-scna(ID="TCGA-05-4405-01A-21D-1856-05",ss=ss,ref_genes="paper",cn="Amp",
#' log2file="analysis/CONUMEE/log2/TCGA-05-4405-01A-21D-1856-05_log2.txt",
#' segfile="analysis/CONUMEE/Segments_TCGA-05-4405-01A-21D-1856-05.txt")
#' res
scna<-function(ID,log2file,segfile,ss,cn,ref_genes="all"){
  bp_stopifnot = getFromNamespace("stopifnot", "backports")
  df<-data.frame(ID=NULL,Int=NULL,X=NULL,Var=NULL)
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
  if(length(int)<1){
    df<-data.frame(ID=NULL,Int=NULL,X=NULL,Var=NULL)
    warning("There are no ",cn," genes for sample: ",ID)
    return(df)
    }
  X=mean(int$log2r,na.rm=T)
  df<-data.frame(
    ID=ID,
    Int=baseline,
    X=X,
    Var=Var)
  return(df)
}



#rlang::expr(paste0(conumee.folder,'/',seg.folder,'/',sname,'_Segments.txt'))

# a<-"paste0(conumee.folder,'/',seg.folder,'/',sname,'_Segments.txt')"
# test<-function(folder="myfolder/",std_file){
#   std_file<-unquote(std_file);print(std_file)
#   }

