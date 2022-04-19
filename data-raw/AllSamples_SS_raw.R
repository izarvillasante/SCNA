## Code to generate Sample sheet for all samples with raw data `all_CN_ASCAT` dataset:
library(parallel)
library(doParallel)
library(TCGAbiolinks)
library(data.table)

# Query for Primary Tumors:

data("all_CN_ASCAT")

ncores_avail<-detectCores()
if(ncores_avail > 18) ncores <- 18 else(ncores <- ncores_avail -2)
cl <- makeCluster(ncores,outfile="")
registerDoParallel(cl)
qfiles_all <- list()
query_all <-list()
TCGA_Project<-unique(all_CN_ASCAT$Project)
a<-foreach(can=TCGA_Project,
           .combine="rbind",
           .errorhandling = "pass",
           .packages =c("TCGAbiolinks","data.table")

) %dopar% {
  query <- tryCatch(GDCquery(project = can,
                             data.category = "Raw microarray data",
                             data.type = "Raw intensities",
                             experimental.strategy = "Methylation array",
                             legacy = TRUE,
                             file.type = ".idat",
                             platform = "Illumina Human Methylation 450",
                             barcode = all_CN_ASCAT[Project==can, barcodes],
                             #sample.type = "Primary Tumor"
                             ),
                    error = function(e) query=NULL
  )
  if(!is.null(query)){
    # Get results:

    qfiles <- getResults(query,cols=c("cases","file_name","sample_type"))
    query_all[[can]]<-qfiles
    # Add cancer type:
    qfiles$Project <- can
    # DOWNLOAD FOLDER:
    cdir <- normalizePath(paste("data-raw","GDC",can,sep="/"))
    if ( ! file.exists(cdir) ) {
      dir.create(cdir,recursive = TRUE);
    }
    tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20,directory = cdir),
             error = function(e) print("no files from api; method=client fails"))
    flist<-list.files(cdir,recursive = T)
    filepaths<-flist[basename(flist[endsWith(x = flist,suffix = ".idat")]) %in% qfiles$file_name]
    setDT(qfiles)
    setkey(qfiles,"file_name")
    # Add the variable old path with the current paths of files:
    qfiles[basename(flist),oldpath:=paste(cdir,flist,sep=.Platform$file.sep)]
    # New path to make it more comprehensible:
    qfiles[,newpath:=paste(cdir,cases,file_name,sep=.Platform$file.sep)]

    qfiles$barcodes<-substr(qfiles$cases,1,12)
    # Check download of those samples were cn ASCAT calls are available.
    # Generate link:
    apply(qfiles,1,function(x) R.utils::createLink(link=x[6], target=x[5]))
    qfiles_all[[can]]<-qfiles
  }
}
stopCluster(cl)
# a$barcodes %in% all_CN_ASCAT$barcodes
# table(all_CN_ASCAT$barcodes %in% a$barcodes)
# all_CN_ASCAT[which(!(all_CN_ASCAT$barcodes %in% a$barcodes)),barcodes]->missing

# Amend file names:
f<-unlist(str_split(
  ifelse(endsWith(a$newpath,"_Grn.idat"),a$newpath,"")
  ,"_Grn.idat"))
filenames<-f[nzchar(f)]

ss_query<-unique(a[,-c("newpath","oldpath","file_name")])
ss_query$filenames=filenames

data("CancerGenes")
common <- intersect(CancerGenes$name,names(all_CN_ASCAT))
genes <- all_CN_ASCAT[,.SD,.SDcols=c("Project","barcodes",common)]
setkey(genes,"barcodes")
ss <- merge(ss_query,genes,by=c("barcodes","Project"))
ascat_genes <- ss[,.SD,.SDcols=common]

AllSamples_SS_raw<-data.table(
  Sample_Name=ss$cases,
  filenames=ss$filenames,
  Cancer=substr(ss$Project,6,nchar(ss$Project)),
  Purity_RFPurify=NA,
  HomDel= apply(ascat_genes,1,function(x) paste(names(ascat_genes)[x == 0],collapse = ";")),
  HetLoss= apply(ascat_genes,1,function(x) paste(names(ascat_genes)[x == 1],collapse = ";")),
  Diploid= apply(ascat_genes,1,function(x) paste(names(ascat_genes)[x == 2],collapse = ";")),
  Gains= apply(ascat_genes,1,function(x) paste(names(ascat_genes)[x %in% 3:4],collapse = ";")),
  Amp= apply(ascat_genes,1,function(x) paste(names(ascat_genes)[x %in% 5:10],collapse = ";")),
  Amp10 = apply(ascat_genes,1,function(x) paste(names(ascat_genes)[x > 10],collapse = ";")),
  Sample_Plate =NA,
  Sample_Group = ss$sample_type,
  Pool_ID=NA,
  Project=NA,
  Sample_Well =NA

)


usethis::use_data(AllSamples_SS_raw, overwrite = TRUE)
