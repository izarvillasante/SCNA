## Code to prepare Genes list:

# From all the genes in this study a refined list of candidate genes
# has been manually selected as described in the paper:
library(data.table)
library(readr)
library(dplyr)
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
data("all_CN_ASCAT")
file="Data/CancerGenes.csv"
bed <- read.csv(file, stringsAsFactors = F)
bed <- bed[!bed$chr%in%c("chrX","chrY"),]
bed <- bed[!duplicated(bed$name),]
common <- intersect(names(all_CN_ASCAT) , bed$name)
genes <-all_CN_ASCAT[,.SD,.SDcols=c("Project","barcodes",common)]
CancerGenes<-bed[bed$name %in% common,]
usethis::use_data(CancerGenes,overwrite = T)


writexl::write_xlsx(CancerGenes,"Supplementary_table_S3.xlsx")
#Generate genomic ranges df for full genes list of genes with ASCAT calls:
all_genes<-names(all_CN_ASCAT)[-c(1:3)]                          #gene names
genome <- TxDb.Hsapiens.UCSC.hg19.knownGene                      #ref genome
egid <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                              all_genes, c("ENTREZID"),"SYMBOL") #Annotation

plec_gene = genes(genome)[which(genes(genome)$gene_id %in% egid$ENTREZID),]
egid_all_genes<-egid[match(plec_gene$gene_id,egid$ENTREZID),]
plec_gene$name<-egid_all_genes$SYMBOL
#Remove sex chrs:
plec_gene <- plec_gene[!seqnames(plec_gene)%in%c("chrX","chrY")]
rg<-ranges(plec_gene)
AllGenes<-data.frame(chr=seqnames(plec_gene)[!seqnames(plec_gene)%in%c("chrX","chrY")],
                     start=start(rg),
                     end=end(rg),
                     strand=strand(plec_gene),
                     name=plec_gene$name,
                     id=plec_gene$gene_id

                     )
usethis::use_data(AllGenes,overwrite = T)
