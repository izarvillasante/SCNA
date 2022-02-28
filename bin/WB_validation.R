library(data.table)
setDTthreads(0)
library(fst)
library(conumee)
library(ChAMPdata)
library(dplyr)

source("Data/functions.R")

ss_normal <- readRDS("raw/normal_ss.rds")
#Select 3 arrays from each cancer type:
ss_control<-ss_normal[,.SD[1:10],Cancer][!(is.na(Sample_Name)),]
ss_control[,sum(table(unique(Sample_Name))),by=Cancer]
ss_control$basename<-basename(ss_control$filenames)
for (i in 1:nrow(ss_control)){

  file<- ss_control$basename[i]
  name<-substr(file,1,nchar(file)-9)
  Project <- paste0("TCGA-",ss_control$Cancer[i])
  ID <- ss_control$Sample_Name[i]
  folder<-"/data/20211210CUP/3/"
  filename<-paste0(folder,Project,"/",ID,"/",name)
  ss_control$filenames[i]<-filename
}
ss_control$basename<-NULL
ss_control<-unique(ss_control)
mset_file<-"~/Documents/Projects/20211210CUP/analysis/ChAMP/BW_normal/mSetSqn.rds"
if(file.exists(mset_file)){control_mset<-readRDS(mset_file)}else{
  intensities <- pre_process(ss_control,subf = "BW_normal")
  intensities$probeid<-rownames(intensities)
  fst::write_fst(intensities,"BW_normal_intensities.fst")
  control_mset <- readRDS("~/Documents/Projects/20211210CUP/analysis/ChAMP/BW_normal/mSetSqn.rds")
}
#M_control<-control_mset@assays@data@listData$M
# genes ranges:
CancerGenes<-readxl::read_xlsx("Supplementary_table_S3.xlsx")%>%as.data.frame()
dd <- toGRanges(CancerGenes)

#### WB:

blood_mset<-readRDS("raw/blood_mset.rds")

### Normals:
#control_mset <- readRDS("~/Documents/Projects/20211210CUP/analysis/ChAMP/BW_normal/mSetSqn.rds")

### Common probes:
cgcommon<-intersect(rownames(control_mset),rownames(blood_mset))
control_mset<-control_mset[cgcommon,]
blood_mset<-blood_mset[cgcommon,]

betas<-getBeta(control_mset)
colnames(betas)<-colData(control_mset)$Cancer
heat_betas<-data.frame(type="normal")
props=data.frame(props = t(table(names(heat_betas)[-1])))
writexl::write_xlsx(props,"cancer_props_in_normalTCGA.xlsx")
for(i in 1:length(colnames(betas))){
  gr<-granges(control_mset)
  gr$betas<-betas[,i]
  int_all<-suppressWarnings(findOverlaps(gr,dd))
  gr.matched_all <- gr[queryHits(int_all)];
  mcols(gr.matched_all) <- cbind.data.frame(
    mcols(gr.matched_all),
    mcols(dd[dd$name%in%dd$name][subjectHits(int_all)])
  )
  heat_betas<-cbind(heat_betas,gr.matched_all$betas)
  
}
colnames(heat_betas)[-1]<-colData(control_mset)$Cancer

my_heatmap<-pheatmap::pheatmap(heat_betas[-1],main = "Beta values on Normal solid Tissue.")
pheatmap::pheatmap(heat_betas[-1],main = "Beta values on Normal solid Tissue.")
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "BW_normal_heatmap.png")


betas_blood<-getBeta(blood_mset)
colnames(betas_blood)<-rownames(colData(blood_mset))
heat_betas_BW<-data.frame(type="Whole_Blood")
for(i in 1:length(colnames(betas_blood))){
  gr<-granges(blood_mset)
  gr$betas<-betas_blood[,i]
  int_all<-suppressWarnings(findOverlaps(gr,dd))
  gr.matched_all <- gr[queryHits(int_all)];
  mcols(gr.matched_all) <- cbind.data.frame(
    mcols(gr.matched_all),
    mcols(dd[dd$name%in%dd$name][subjectHits(int_all)])
  )
  heat_betas_BW<-cbind(heat_betas_BW,gr.matched_all$betas)
  
}
colnames(heat_betas_BW)[-1]<-rownames(colData(blood_mset))

my_heatmap<-pheatmap::pheatmap(heat_betas_BW[-1],main = "Beta values on Whole Blood Tissue.")
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "BW_Whole_Blood_heatmap.png")

heat_data<-cbind(heat_betas[-1],heat_betas_BW[-1])

saveRDS(heat_data,"heat_data.rds")
heat_data<-readRDS("heat_data.rds")
type<-data.frame(type=c(rep("TCGA normal",75),rep("whole blood",96)))
rownames(type)<-colnames(heat_data)
my_heatmap<-pheatmap::pheatmap(heat_data,
                               main = 'DNA methylation TCGA normal and whole blood',
                               annotation_col = type,
                               cluster_cols = T,show_colnames = F,
                               show_rownames = F,
                               
                               )

save_pheatmap_png(my_heatmap, "BW_normal_heatmap.png")

cellType<-c(rep("WB",ncol(heat_betas_BW)-1),rep("N",ncol(heat_betas)-1))%>%factor 
individual<-factor(c(rep("WB",ncol(heat_betas_BW)-1),colData(control_mset)$Cancer))

combined<-combineArrays(blood_mset,control_mset)
rn<-rownames(mcols(gr.matched_all))
combined<-combined[rn,]
rownames(colData(combined))[97:171]<-colnames(heat_betas)[-1]
colData(combined)$Sample_Group<-cellType
colData(combined)$Sample_Name[1:96]<-rownames(colData(combined)[1:96,])
targets<-colData(combined)
design <- model.matrix(~0+cellType, data=targets)

library(limma)

mVals <-getM(combined)
fit <-lmFit(mVals,design)
fit2<-eBayes(fit)

DMPs<-topTable(fit2,num=Inf)

d<-DMPs[rownames(DMPs) %in% rn,]
summary(decideTests(fit2))

test<-decideTests(fit2)

class(test)
test<-as.data.frame(test)
test$probeid<-rownames(test)
setDT(test)
test[,same:=ifelse(cellTypeN==cellTypeWB,T,F)]
table(test[test$probeid%in%rn,"same"])/1307*100
intersect(test$cellTypeN,test$cellTypeWB)



m2<-model.matrix(~0+cellType+individual, data=targets)
m2<-m2[,-length(colnames(m2))]
colnames(m2)<-c(levels(cellType),levels(individual)[-c(1,length(levels(individual)))])
fit_m2<-lmFit(mVals,m2)
contMatrix<-makeContrasts(WB-N,levels=m2)
fit2_m2 <- contrasts.fit(fit_m2, contMatrix)
fit2_m2 <- eBayes(fit2_m2)
summary(decideTests(fit2_m2))


###Minfi
dmp<-dmpFinder(dat=getM(combined),pheno = cellType)
NROW(dmp[dmp$qval<0.01 & rownames(dmp) %in% rn,])/1307
dmp_beta<-dmpFinder(dat=getBeta(combined),pheno = cellType)
NROW(dmp_beta[dmp_beta$qval<0.001 & rownames(dmp_beta) %in% rn,])/1307

m<-as.data.frame(getBeta(combined))
m$rns<-rownames(m)
setDT(m)
m[,blood_mean:=rowMeans(.SD),.SDcols=colnames(m)[1:97]]
m[,normal_mean:=rowMeans(.SD),.SDcols=colnames(m)[98:171]]
m[,diff_mean:=abs(blood_mean -normal_mean)]

setkey(m,rns)
m[as.numeric(rownames(dmp)),qval:=dmp$qval]
# Results:
# With minfi 0.33
dif_033<-m[qval<0.001 & diff_mean>=0.33 & rns %in% rn,]
100-NROW(dif_02)/NROW(m)*100

# With minfi 0.2
dif_02<-m[qval<0.001 & diff_mean>=0.2 & rns %in% rn,]
100-NROW(dif_02)/1304*100
# With minfi 0.1 
dif_01<-m[qval<0.001 & diff_mean>=0.1 & rns %in% rn,]
100-NROW(dif_01)/1304*100
# With limma decideTests:
table(test[test$probeid%in%rn,"same"])/1307*100



actionable<-c("MYC","KRAS","CCND1","FGF3","FGF4","CCNE2","PIK3CA","MET","BCL6",
              "FGFR1","RICTOR","EGFR","PIK3CB","CCNE1","ERBB2")

library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genome <- TxDb.Hsapiens.UCSC.hg19.knownGene

egid <-
  AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                        actionable,
                        c("ENTREZID"),
                        "SYMBOL")


#select(genome,keys="100",columns=columns(genome),keytype = "GENEID")
# the plec gene
plec_gene = genes(genome)[which(genes(genome)$gene_id %in% egid$ENTREZID),]
egid_actionable<-egid[match(plec_gene$gene_id,egid$ENTREZID),]
plec_gene$name<-egid_actionable$SYMBOL

# Get beta values of the 15 actionable genes for each sample: 
action_betas<-NULL
all_betas<-getBeta(combined)
for(i in 1:length(colnames(all_betas))){
  gr<-granges(combined)
  gr$betas<-all_betas[,i]
  int_all<-suppressWarnings(findOverlaps(gr,plec_gene))
  gr.matched_all <- gr[queryHits(int_all)];
  mcols(gr.matched_all) <- cbind.data.frame(
    mcols(gr.matched_all),
    mcols(plec_gene[subjectHits(int_all)])
  )
  action_betas<-cbind(action_betas,gr.matched_all$betas)
  
}

colnames(action_betas)<-colnames(all_betas)
df_action<-as.data.frame(action_betas)
###Minfi
# Calculate p-values of dmps with minfi::dmpfinder:
dmp_beta<-dmpFinder(dat=action_betas,pheno = cellType)
NROW(dmp_beta[dmp_beta$qval<0.001 & rownames(dmp_beta) %in% rn,])/1307

# Create a variable with diff of mean meth value between blood and normal: 
setDT(df_action)
df_action[,blood_mean:=rowMeans(.SD),.SDcols=colnames(df_action)[1:96]]
df_action[,normal_mean:=rowMeans(.SD),.SDcols=colnames(df_action)[97:171]]
df_action[,diff_mean:=abs(blood_mean -normal_mean)]
# Get rownames:
# df_action$rns<-names(gr.matched_all)
# setkey(df_action,rns)
df_action[as.numeric(rownames(dmp_beta)),qval:=dmp_beta$qval]

# Results:
# With minfi 0.2
dif_02<-df_action[qval<0.001 & diff_mean>=0.2 & rns %in% rn,]
100-NROW(dif_02)/177*100
# With minfi 0.1 
dif_01<-df_action[qval<0.001 & diff_mean>=0.1 & rns %in% rn,]
100-NROW(dif_01)/177*100
# With limma decideTests:
table(test[test$probeid%in%rn,"same"])/1307*100
