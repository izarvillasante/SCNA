## Upload Copy Number Polymorphism file from Broad Inst. for their exclusion

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
file <- "../20211210CUP/raw/Data/CNV_Germline_GSITIC2_BROAD_SNP6.merged.151117.hg19.CNV.txt"
cnvbroad <- read.delim(file,header=T,sep="\t")

## Add 'chr' tag in front of chr numbers
cnvbroad <- data.frame(chr=paste("chr",cnvbroad$Chromosome,sep = ""), start=cnvbroad$Start, end=cnvbroad$End)

## make Granges object:
cnvbroadgr <- regioneR::toGRanges(cnvbroad)
Exclude <- cnvbroadgr

#Remove chr23 (X) and chr24 (Y)
Exclude <- GenomeInfoDb::keepSeqlevels(Exclude, paste("chr",c(1:22), sep=""), pruning.mode="coarse")
seqnames(Exclude)

## Bed file CancerGenes most amplified/most deleted in cancer autosomal chromosomes
data("CancerGenes")
detail_region <- regioneR::toGRanges(CancerGenes)
seqlevels(detail_region)
anno <- conumee::CNV.create_anno(array_type = "450k",chrXY = F,
                        exclude_regions = Exclude, detail_regions = detail_region)
#saveRDS(anno, "./data-raw/anno.rds")
#usethis::use_data(anno, overwrite = TRUE,internal = TRUE)
