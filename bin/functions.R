################# Declare functions:
library(filesstrings)
library(preprocessCore)
library(ChAMP)
library(ChAMPdata)
library(DNAcopy)
library(data.table)

library(minfi)
library("minfiData")
## Exclude cross reactive probes
library(devtools)
#install_github("markgene/maxprobes", force=T)
library(conumee)
library(minfi)
library("minfiData")
library(regioneR)
library(filesstrings)


pre_process<-function(ss,subf,purity=T,query=T,out="analysis/ChAMP/"){
  message("enanos trabajando") 
  
  ##Create dir ChAMP/n/. put all the sample sheet idat files here
  ChAMP_folder<-paste0(out,subf)
  #ChAMP_folder<-paste(getwd(),ChAMP_folder,sep="/")

  dir.create(ChAMP_folder)
  #empty folder
  file.remove(list.files(ChAMP_folder, full.names = TRUE))
  #copy files
  file.copy(paste0(ss$filenames,"_Grn.idat"),ChAMP_folder)
  file.copy(paste0(ss$filenames,"_Red.idat"),ChAMP_folder)
  message(ss$Basename)
  ## this directory, you should copy the file "Sample_Sheet.csv", as ChAMP will look for it.
  write.csv(ss,paste0(ChAMP_folder,"/Sample_sheet.csv"))
  data(hm450.manifest.hg19 )
  #Run ChAMP with method for QC minfi
  myLoad <- ChAMP::champ.load(ChAMP_folder, method="minfi")
  #read.metharray.exp(base = ChAMP_folder)->a
  ss<-myLoad$pd
  #save(myLoad,file = paste0(ChAMP_folder,"TCGA_New.RData"))#Give it the name you want.
  message("Champ does load data")
  if(purity==T){
    purity<-purify(myLoad=myLoad)
    message("Rfpurifies")
    ss$Purity_Impute_RFPurity.Absolute. <- purity
    write.table(ss,paste(ChAMP_folder,"Sample_Sheet.txt",sep="/"), col.names = T, row.names = F, quote = F, sep="\t")
    message(paste("ss is saved ",subf))
  }
  if(query==T){
    message("start query")
    
    query <- queryfy(myLoad$rgSet,ss=ss,ChAMP_folder=ChAMP_folder)
    message("end querry")
    saveRDS(query,paste0(ChAMP_folder,"/intensities.rds"),compress = FALSE)
    message("eureka!")
    return(query)
  }
}

purify<- function(myLoad){
  library(RFpurify)
  #impute missing data
  impute <- ChAMP::champ.impute(beta = getBeta(myLoad$rgSet),pd = myLoad$pd,method="KNN")
  #?champ.impute
  dim(impute$beta)
  
  ##Apply RFPurify
  #estimate <- predict_purity_betas(impute$beta,method="ESTIMATE")
  
  absolute <- RFpurify::predict_purity_betas(impute$beta,method="ABSOLUTE")
  
  return(absolute)
}

queryfy<-function(rgSetn,ss,ChAMP_folder){
  
  #myLoad$rgSet #read.metharray(samps$Basename,force=T)
  mSetSqn <- tryCatch( preprocessQuantile(rgSetn),error=function(e) {message("failed")
    return(NULL)})
  
  # calculate p-values
  detP <- detectionP(rgSetn)
  bad <- colnames(detP)[colSums(detP >=0.01)/nrow(detP) > 0.1]
  
  ## Ensure probes are in the same order in the mSetSqn and detP objects
  detP <- detP[match(featureNames(mSetSqn),rownames(detP)),]
  
  ## Remove rows with at elast one 'NA' entry
  keep <- rowSums(detP < 0.01) == ncol(mSetSqn)
  #table(keep)
  mSetSqn <- mSetSqn[keep,]
  
  ## If your data includes males and females, remove probes on the sex chromosomes
  FDATA450 = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  keep <- !(featureNames(mSetSqn) %in% FDATA450$Name[FDATA450$chr %in% c("chrX","chrY")])
  table(keep)
  mSetSqn <- mSetSqn[keep,]
  
  ## Remove probes with SNPs at CpG site
  mSetSqn <- dropLociWithSnps(mSetSqn)
  mSetSqn <-  maxprobes::dropXreactiveLoci(mSetSqn)
  saveRDS(mSetSqn,paste0(ChAMP_folder,"/mSetSqn.rds"),compress = FALSE)
  
  mSetSqn <- CNV.load(mSetSqn)
  return(mSetSqn@intensity)
}

run_conumee<-function(anno,controls, ss, k, fst.file="intensities.fst",out="analysis/CONUMEE/"  ){
  Sample_Name<-ss$Sample_Name[k]
  
  message(k)
  #if(!(file.exists(paste(out,"log2/",Sample_Name,"_log2.txt", sep="")))){
  log2_file<-paste(out,"log2/",Sample_Name,"_log2.txt", sep="")
  message(log2_file)
  
  #############################3
  my.data<- read.fst(fst.file,c("probeid",Sample_Name))
  ##############################
  
  message(my.data[1,])
  #rownames(my.data)<-my.data$probeid
  names(my.data)[2]<-"intensity"
  
  cgcommon<-intersect(my.data$probeid,rownames(controls@intensity))%>% intersect(.,names(anno@probes@ranges))
  my.data <- my.data[my.data$probeid %in% cgcommon,]
  controls@intensity <-controls@intensity[cgcommon,]
  anno@probes<-anno@probes[names(anno@probes)%in% cgcommon,]
  
  my.datacnv<-CNV.load(my.data)
  rownames(my.datacnv@intensity)<-my.data$probeid
  fit  <- CNV.fit(my.datacnv, controls, anno)
  #Log2 <- list(log2ratio=fit@fit$ratio, Purity=ss$Purity_Impute_RFPurify.Absolute.[ss$Sample_Name%in%Sample_Name])
  fit2 <- CNV.segment(CNV.detail(CNV.bin(fit)))
  write.table(CNV.write(fit2, what = "segments"),sprintf("%sSegments_%s.txt",out,Sample_Name),col.names = T, row.names = F, quote = F, sep = "\t") 
  message("Sgements file" ," saved")
  log2ratio<-as.data.frame(fit@fit$ratio)
  names(log2ratio) <- Sample_Name
  #log2ratio$purity <- ss[ss$Sample_Name== Sample_Name,Purity_Impute_RFPurity.Absolute.]
  write.table(log2ratio,log2_file)
  message(log2_file ," saved")
  return(log2ratio)
  #}else return(NULL) 
}





scna<-function(ID,segfolder="analysis/CONUMEE/",samp,dd,cn){
  
  message("if")
  # Load complete log2r of whole array:
  log2file_list <- list.files(paste0(segfolder,"log2/"))
  log2file <- log2file_list[startsWith(log2file_list,ID)][1]
  Data<-list()
  if(isEmpty(log2file)){
    message("please provide a valid log2file you can use run_conumee.")
    message(paste0("no log2file for " , ID))
  }else{
    
    Log2file<-paste0(segfolder,"log2/",log2file)
    message("Log2file: ", Log2file)
    Log2<-fread(Log2file,col.names = c("probeid","log2r"))
    
    # Calculate CNCall:
    baseline <- mean(Log2$log2r, na.rm=T)
    purity<-names(samp)[names(samp) %ilike% "impute" & names(samp) %ilike% "purity"]
    #purity<-names(samp)[names(samp) %ilike% "Purit" & names(samp) %ilike% "CCLE"]
    p<-with(samp,get(purity))
    
    Var <- p*sd(Log2$log2r,na.rm=T)
    
    
    # Load segmented data and make Granges:
    segfile_list <- list.files(segfolder)
    segfile <- segfile_list[startsWith(segfile_list,paste0("Segments_",ID))][1]
    
    if(isEmpty(segfile)){
      message("please provide a valid segment file. you can use run_conumee.")
      message(paste0("no segfile for " , ID))
    }else{
      
      seg<-read.delim(paste0(segfolder,segfile),header=T,sep="\t")
      genes<-strsplit(as.character(with(samp,get(cn))),split = ";")[[1]]
      if(isEmpty(genes)){df<-data.frame(ID=NULL,Int=NULL,X=NULL,Var=NULL)}else{
        segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean)
        seggr <- toGRanges(segb)
        int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
        seggr.matched <- seggr[queryHits(int)];
        mcols(seggr.matched) <- cbind.data.frame(
          mcols(seggr.matched),
          mcols(dd[dd$name%in%genes][subjectHits(int)]));
        
        int <- seggr.matched
        message("int finish")
        X=mean(int$log2r,na.rm=T)
        df<-data.frame(
          ID=ID,
          Int=baseline,
          X=X,
          Var=Var)
      }
      
    }
  }
  return(df)
}

conumee_kc<-function(ID,segfolder="analysis/CONUMEE/",samp,dd,
                     K=list(Amp10=3.666485,Amp=1.495666,Gains=1.096756,HetLoss=-1.887569,HomDel=-5.124848)){
  library(dplyr)
  library(regioneR)
  #message(samp$Sample_Name)
  # Load complete log2r of whole array:
  log2file_list <- list.files(paste0(segfolder,"log2/"))
  log2file <- log2file_list[startsWith(log2file_list,ID)][1]
  Data<-list()
  if(isEmpty(log2file)){
    
    message(paste0("no log2file for " , ID))
  }else{
    
    Log2file<-paste0(segfolder,"log2/",log2file)
    message("Log2file: ", Log2file)
    Log2<-fread(Log2file,col.names = c("probeid","log2r"))
    
    # Calculate CNCall:
    baseline <- mean(Log2$log2r, na.rm=T)
    purity<-names(samp)[names(samp) %ilike% "impute" & names(samp) %ilike% "purity"]
    p<-with(samp,get(purity))
    
    Var <- p*sd(Log2$log2r,na.rm=T)

    K$Diploid<-unname(unlist(K[names(which.max(K[K<0]))])+0.00001)
    Tresholds<-sapply(K, function(x) baseline + unname(x) * Var)
    
    # Load segmented data and make Granges:
    segfile_list <- list.files(segfolder)
    segfile <- segfile_list[startsWith(segfile_list,paste0("Segments_",substr(ID,1,16)))][1]
    
    if(isEmpty(segfile)){
      message(paste0("no segfile for " , ID))
    }else{
      
      seg<-read.delim(paste0(segfolder,segfile),header=T,sep="\t")

      for (i in 1:nrow(seg)){
        rest <- seg$seg.mean[i]-Tresholds
        CNA <- names(which.min(rest[rest>0]))
        if(isEmpty(CNA)){ # In the case all Thresholds are positive and seg.mean is below
          CNA<-names(which.max(rest))
          
        }
        # print(i)
        # print(CNA)
        seg$CNCall[i]=CNA
      }
      
      segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean, CNCall=seg$CNCall)
      seggr <- toGRanges(segb)
      int_all<-suppressWarnings(findOverlaps(seggr,dd))
      seggr.matched_all <- seggr[queryHits(int_all)];
      mcols(seggr.matched_all) <- cbind.data.frame(
        mcols(seggr.matched_all),
        mcols(dd[dd$name%in%dd$name][subjectHits(int_all)]));
      geneCall<-table(seggr.matched_all$CNCall)
      for (cn in c(names(K),"GAINS")){
        if(!(cn %in% names(samp))){
          message(cn, " not found in sample_sheet column names.")
          next
        }
        genes<-strsplit(as.character(with(samp,get(cn))),split = ";")[[1]]
        #print(genes)
        int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
        seggr.matched <- seggr[queryHits(int)];
        mcols(seggr.matched) <- cbind.data.frame(
          mcols(seggr.matched),
          mcols(dd[dd$name%in%genes][subjectHits(int)]));
        TP=0
        #apply(res,1,function(x),)
        int<-seggr.matched
        int <- seggr.matched[!duplicated(seggr.matched$name)]

        if (cn == "GAINS")CN<-names(K[K>0])else  if (cn == "LOSS")CN<-names(K[K<0])else CN<-cn
        if(length(int) >=1){
          for(i in 1:length(int)){
           
            res <- int[i,]
            
            dat <- filter(segb, chr%in%as.character(res@seqnames@values) & start <= res@ranges@start & end >= (res@ranges@start + res@ranges@width-1) & CNCall%in%CN)
            if(nrow(dat)!=0)TP=TP+1;
          }
        }
        
        FP<-sum(geneCall[CN])-TP
        Data[[cn]][[ID]]=list(TP=TP,
                              FP=FP,
                              NumGenes=length(unique(seggr.matched$name)),
                              TPR=TP/length(unique(seggr.matched$name)),
                              FPR=FP/(length(dd$name)-length(unique(seggr.matched$name)))
        )
        message(samp$Sample_Name, ": Success")
      }
    }
  }
  
  return(Data)
}

conumee_st<-function(treshold=0.3,ID,segfolder="analysis/CONUMEE/",samp,dd,gains_col="ASCAT_GAINS",losses_col="ASCAT_LOSS"){
  
  message("if")
  # Load complete log2r of whole array:
  
  Data<-list()
  
  Tresholds<-c(GAINS=treshold,LOSS=-treshold)
  
  # Load segmented data and make Granges:
  segfile_list <- list.files(segfolder)
  segfile <- segfile_list[startsWith(segfile_list,paste0("Segments_",substr(ID,1,16)))][1]
  
  if(isEmpty(segfile)){
    message(paste0("no segfile for " , ID))
  }else{
    
    seg<-read.delim(paste0(segfolder,segfile),header=T,sep="\t")
    
    seg$CNCall="Diploid"
    for(i in 1:nrow(seg)){
      if(seg$seg.mean[i] >= Tresholds["GAINS"] )seg$CNCall[i] <- "GAINS"
      
      if(seg$seg.mean[i] <= Tresholds["LOSS"] )seg$CNCall[i] <- "LOSS"
    }
    
    segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean, CNCall=seg$CNCall)
    seggr <- toGRanges(segb)
    int_all<-suppressWarnings(findOverlaps(seggr,dd))
    seggr.matched_all <- seggr[queryHits(int_all)];
    mcols(seggr.matched_all) <- cbind.data.frame(
      mcols(seggr.matched_all),
      mcols(dd[subjectHits(int_all)]));
    geneCall<-table(seggr.matched_all$CNCall)
    # message(geneCall)
    # message(samp[,13:14])
    cnstates<-c(gains_col,losses_col)
    for (cn in cnstates){
      genes<-strsplit(as.character(with(samp,get(cn))),split = ";")[[1]]
      #print(genes)
      int <- suppressWarnings(findOverlaps(seggr,dd[dd$name%in%genes]))
      seggr.matched <- seggr[queryHits(int)];
      mcols(seggr.matched) <- cbind.data.frame(
        mcols(seggr.matched),
        mcols(dd[dd$name%in%genes][subjectHits(int)]));
      TP=0
      #apply(res,1,function(x),)
      int<-seggr.matched
      int <- seggr.matched[!duplicated(seggr.matched$name)]
      
      
      if(length(int) >=1){
        for(i in 1:length(int)){
          #print(i)
          res <- int[i,]
          
          dat <- filter(segb, chr%in%as.character(res@seqnames@values) & start <= res@ranges@start & end >= (res@ranges@start + res@ranges@width-1) & CNCall%in%cn)
          #print(dat)
          if(nrow(dat)!=0)TP=TP+1;
        }
      }
      
      FP<-sum(geneCall[cn])-TP
      Data[[cn]][[ID]]=list(TP=TP,
                            FP=FP,
                            NumGenes=length(unique(seggr.matched$name)),
                            TPR=TP/length(unique(seggr.matched$name)),
                            FPR=FP/(length(dd$name)-length(unique(seggr.matched$name)))
      )
    }
  }
  
  return(Data)
}

make_Kc<-function(ss,cncols,fname){
  K_list<-list()
  cl <- makeCluster(detectCores()-1,outfile="")
  registerDoParallel(cl)
  for (cnstate in cncols){#c("ASCAT_Gain","ASCAT_AMP","ASCAT_AMP10")){
    Kc_file<-paste0(fname,"_",cnstate,".rds")
    if(file.exists(Kc_file)){K<-readRDS(Kc_file)}else{
      scnas_file<-paste0("scnas_",fname,"_",cnstate,".rds")
      #if(file.exists(scnas_file)){readRDS(scnas_file)}
      
      
      all_scnas <-
        foreach(sname=ss$Sample_Name,
                .combine='rbind', 
                .inorder=FALSE,
                .errorhandling = "pass",
                .packages =c("data.table","dplyr","regioneR")
        ) %dopar% {
          scna(cn = cnstate, ID=sname,samp=ss[Sample_Name==sname,][1,],dd=dd)
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
  writexl::write_xlsx(Ks, paste0(name,"_Ks.xlsx"))
}

get_genes<-function(ID,segfolder="analysis/CONUMEE/",samp,dd,
                    K=list(Amp10=3.666485,Amp=1.495666,Gains=1.096756,HetLoss=-1.887569,HomDel=-5.124848)){
  library(dplyr)
  library(GenomicRanges)
  message("if")
  # Load complete log2r of whole array:
  log2file_list <- list.files(paste0(segfolder,"log2/"))
  log2file <- log2file_list[startsWith(log2file_list,ID)][1]
  Data<-data.frame()
  if(isEmpty(log2file)){
    message()
    message(paste0("no log2file for " , ID))
  }else{
    
    Log2file<-paste0(segfolder,"log2/",log2file)
    message("Log2file: ", Log2file)
    Log2<-fread(Log2file,col.names = c("probeid","log2r"))
    
    # Calculate CNCall:
    baseline <- mean(Log2$log2r, na.rm=T)
    purity<-names(samp)[names(samp) %ilike% "impute" & names(samp) %ilike% "purity"]
    p<-unname(with(samp,get(purity)))
    
    Var <- p*sd(Log2$log2r,na.rm=T)
    
    K$Diploid<-unname(unlist(K[names(which.max(K[K<0]))])+0.00001)
    Tresholds<-sapply(K, function(x) baseline + x * Var)
    #write.table(Tresholds,paste0("tre_",i,".csv"))
    # Load segmented data and make Granges:
    segfile_list <- list.files(segfolder)
    segfile <- segfile_list[startsWith(segfile_list,paste0("Segments_",substr(ID,1,16)))][1]
    
    if(isEmpty(segfile)){
      message(paste0("no segfile for " , ID))
    }else{
      
      seg<-read.delim(paste0(segfolder,segfile),header=T,sep="\t")
      
      for (i in 1:nrow(seg)){
        rest <- seg$seg.mean[i]-Tresholds
        CNA <- names(which.min(rest[rest>0]))
        if(isEmpty(CNA)){ # In the case all Thresholds are positive and seg.mean is below
          CNA<-names(which.max(rest))
          
        }
        print(i)
        print(CNA)
        seg$CNCall[i]=CNA
      }
      
      segb <- data.frame(chr=seg$chrom, start=seg$loc.start, end=seg$loc.end, log2r= seg$seg.mean, CNCall=seg$CNCall)
      seggr <- toGRanges(segb)
      int_all<-suppressWarnings(findOverlaps(seggr,dd))
      seggr.matched_all <- seggr[queryHits(int_all)];
      mcols(seggr.matched_all) <- cbind.data.frame(
        mcols(seggr.matched_all),
        mcols(dd[dd$name%in%dd$name][subjectHits(int_all)]));
      geneCall<-table(seggr.matched_all$CNCall)
      Data<-data.frame(gene_name=seggr.matched_all$name,
                       CNCall= seggr.matched_all$CNCall)
    }
  }
  return(Data)
}