# load Libraries:
library(data.table)
setDTthreads(0)
library(readxl)
# load functions:
source("functions.R")

#Load datasets:
excel<-read_excel("analysis/CUP/CUP results_PBlecua.xlsx")

# genes ranges:
library(readxl)
excel<-read_excel("analysis/CUP/CUP results_PBlecua.xlsx",sheet = 2)
genes_list<-excel[-length(excel),1]
genes_list[[1]]
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genome <- TxDb.Hsapiens.UCSC.hg19.knownGene
egid <-
  AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                        genes_list[[1]],
                        c("ENTREZID"),
                        "SYMBOL")

# the plec gene
plec_gene = genes(genome)[which(genes(genome)$gene_id %in% egid$ENTREZID),]
egid_genes<-egid[match(plec_gene$gene_id,egid$ENTREZID),]
plec_gene$name<-egid_genes$SYMBOL
saveRDS(plec_gene,"CUP_genes.rds")
dd<-readRDS("CUP_genes.rds")
# CONUMEE Calibration with 442 TCG samples (Training) vs BLOOD, 
# Rfpurify purity and Izar code: K_c cnostants:
K_list<-list()
K_list$HomDel<-(-5.124848)
K_list$HetLoss<-(-1.887569)
K_list$Gains<-  1.096756
K_list$Amp10<- 3.666485
K_list$Amp<- 1.495666
library(stringr)
files_list<-list.files("~/shares/BACKUP/raw/Arrays/",recursive = T)
f<-str_split(files_list,"_Grn.idat")%>% unlist
f2<-basename(files_list)
f3<-str_split(f2,"_Grn.idat")
f4<-unlist(f3)
#Find a single filepath:
f[which(f4 %in% excel[excel$CUP_ID=="CUP-123","Sentrix_ID"]  )]
#All paths
filenames<-f[sapply(excel$Sentrix_ID, function(x)which(f4 ==x) )]
excel$filenames<-filenames
setDT(excel)

ss<-data.table(
  Sample_Name=excel$CUP_ID,
  filenames=paste0("/home/idevillasante/shares/BACKUP/raw/Arrays/",excel$filenames),
  Cancer=excel$"CUP",
  Purity_Impute_RFPurity.Absolute.=NULL,
  ASCAT_Homdel= excel$HomDel,
  ASCAT_HetLoss= excel$HetLoss,
  ASCAT_diploid= "",
  ASCAT_Gain= excel$Gain,
  ASCAT_AMP= excel$Amp1,
  ASCAT_AMP10 = excel$Amp2,
  Sample_Plate =NA,
  Sample_Group ="Unknown Primary",
  Pool_ID=NA,
  Project=NA,
  Sample_Well =NA
  
)

if(!file.exists("CUP2_intensities.fst")){
  intensities <- pre_process(ss,subf = "CUP2")
  intensities$probeid<-rownames(intensities)
  fst::write_fst(intensities,"CUP2_intensities.fst")
}

library(data.table)
setDTthreads(0)
library(fst)
library(conumee)
library(dplyr)
#Load annotation file
anno<-readRDS("anno.rds")

#Load controls file
controls<-readRDS("controls.rds")

library(doParallel)
library(parallel)


##################### train:
all_samples<-ss
cl <- makeCluster(4,outfile="")
registerDoParallel(cl)
start_time <- Sys.time()
query_list<-foreach(k=1:NROW(all_samples),
                    .combine="rbind",
                    .errorhandling = "pass",
                    .packages =c("conumee","fst","dplyr")
) %dopar%
  {

    run_conumee(anno=anno, controls=controls, ss = all_samples,k=k,
                fst.file = "CUP2_intensities.fst"  )
  }


stopCluster(cl)
intlist<-list()
ss<-as.data.frame(ss)
setDT(ss)
setkey(ss,"Sample_Name")
for (i in ss$Sample_Name){
  #i<-ss$Sample_Name[21]
  ID<-i
  
  # Adaptative thresholding:
  genes <- get_genes(ID=ID,segfolder="analysis/CONUMEE/",samp= ss[i,],dd)
  message("done CONUME_KC")
  intlist[[i]]<-genes
}
library(tidyverse)
a<-reduce(intlist,full_join,by="gene_name")
setDT(a)
names(a)<-c("gene_symbol",names(intlist))
b<-a[,-1]
for (i in names(K_list)){
  colname<-paste0(i,"_count")
  a[,(colname):=rowSums(b==i,na.rm=T)]
}

a[,Amp_Amp10_count:=Amp_count+Amp10_count]
setorder(a,-Amp_Amp10_count,-Amp10_count)
writexl::write_xlsx(a,"CUP_genes_cnstate.xlsx")
library(readxl)
a<-read_xlsx("curated_CUP_genes_cnstate.xlsx")
setDT(a)
setkey(a,"gene_symbol")
#Añadir pos y chrome band:
setDT(dd)
setkey(dd,"name")

library(GenomicRanges)
library(rtracklayer)
library(biovizBase) #needed for stain information
library(IRanges)
library(R.utils)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

query <- rtracklayer::ucscTableQuery("hg18", "cytoBandIdeo")
table1 <- rtracklayer::getTable(query)
table1$Strand <- c("*")
dd <- GRanges(table1$chrom,
                     IRanges(table1$chromStart, table1$chromEnd),
                     table1$Strand,
                     table1$name, table1$gieStain)

seggr<-readRDS("CUP_genes.rds")
int <- suppressWarnings(findOverlaps(seggr,dd))
seggr.matched <- seggr[queryHits(int)];
mcols(seggr.matched) <- cbind.data.frame(
  mcols(seggr.matched),
  mcols(dd[subjectHits(int)]));
ddextend<-as.data.table(seggr.matched)
setkey(ddextend,"name")
cytoband<-ddextend[unique(a$gene_symbol),.SD,.SDcols=c("name","table1.name")]
cyt<-cytoband[match(unique(cytoband$name),cytoband$name),]
a[cyt$name,cytoband:=cyt$table1.name]
a[cyt$name,cytoband:=cyt$table1.name]
full_a<-cbind(dd[a$gene_symbol,c(1:3)],a$cytoband,a)
full_a<-full_a[,-223]
curated<-full_a[full_a$Amp_Amp10_count>=5,]
writexl::write_xlsx(curated,"Data/updated_CUP_genes_cnstate.xlsx")
