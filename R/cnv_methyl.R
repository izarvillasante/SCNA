#' Generate segment & log2r intensity files calculated by CONUMEE.
#' @param ss Sample sheet containing purities, cnstate and Sample_Name ids.
#' @return data.frame with metadata from segmentation, genomic ranges, genes
#' and scna for each of the genes
#' @export
#'
#' @examples
#' data("TrainingSet_Sample_sheet")
#' ss<-TrainingSet_Sample_sheet[1:14,]
#' cnv_methyl(ss)
#'
cnv_methyl<-function(ss){

  intensity<-pre_process(ss)
  ss<-data.table::fread("analysis/ChAMP/Sample_Sheet.txt")
  run_conumee(intensities = intensity)
  # Cluster:
  ncores<-parallel::detectCores()-2
  if (!is.na(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))))ncores=as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  cl<- parallel::makePSOCKcluster(ncores, outfile="")
  doParallel::registerDoParallel(cl)

  res<-foreach::foreach(i=1:length(ss$Sample_Name),#isplitIndices(1400,chunks=ncores),
                        .combine='c',
                        .multicombine = F,
                        .inorder=F,
                        .packages = "data.table",
                        .export = c("AllGenes","CancerGenes"),
                        .errorhandling = "pass"
  )%dopar%{
    source("R/utils.R")
    source("R/Kc_get.R")
    #source()
     ss<-ss[i,]
     ID<-ss$Sample_Name
     cna <- Kc_get(ss=ss,ID=ID)
    # cna

    }
}
