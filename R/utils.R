check_input_files<-function(Sample_Name,folder,std_file="standard_format"){

  file_list<-list.files(folder,full.names = T)
  path<-sapply(Sample_Name,function(ID){
    match=FALSE
    match<-sapply(file_list, function(file) ifelse(data.table::like(basename(file),ID),TRUE,FALSE))
    if(sum(match) == 0) warning("No input file found for sample: ",ID)
    if(sum(match) > 1){
      #std_file <- paste0(conumee.folder,"/",.folder,"/",ID,"_Segments.txt")
      match_files <- file_list[match]
      std_file<-eval(std_file)
      if(std_file %in% match_files){
        ok <- match_files[which(std_file %in% match_files)]}else{
          ok <- file_list[which.max(match)]}
      overwrite<-setdiff(match_files,ok)
      warning("multiple files match ",ID, ". \n ", paste(overwrite,collapse = ", ")," ignored.")
      return(ok)
    }# otherwise match == 1
    match_files <- file_list[match]
    return(match_files)
  })
  return(path)
}

get_int<-function(seg,ref_genes="all",cn_genes){
  if(ref_genes=="all"){data("AllGenes");interest_geneset <- AllGenes}else if(
    ref_genes == "paper"){data("CancerGenes");interest_geneset <- CancerGenes}else{
      interest_geneset<-utils::read.table(ref_genes)}
  interest_genes<-interest_geneset$name
  genes <- intersect(cn_genes,interest_genes)
  if(length(genes)<1){
    df<-data.frame(ID=NULL,Int=NULL,X=NULL,Var=NULL)
    warning("No genes of interest present in cnstate: ")
    return(df)
  }
  #segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean)
  #seggr <- regioneR::toGRanges(segb)


  old_names<-c("chrom","loc.start","loc.end","seg.mean")
  new_names<-c("chr","start","end","log2r")
  data.table::setnames(seg, old_names,new_names)
  data.table::setcolorder(seg,new_names)
  #cols<-colSums(is.na(seg))<1
  #seg<-seg[,.SD,.SDcols=cols]
  seggr <- regioneR::toGRanges(seg)
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
