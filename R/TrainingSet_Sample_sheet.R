#' Sample sheet for the 442 TCGA samples used as training set.
#'
#' A dataset containing the IDs & paths to the 450k microarrays used
#' as training set to infer the reference KCN as described in the paper.
#' Reference ASCAT data needed to estimate each KCN is extracted from a previous
#' work https://github.com/riazn/biallelic_hr.git.
#'
#' @format A data frame with 442 rows and 10 variables:
#' \describe{
#'   \item{Sample_Name}{Unique identifier of the sample. TCGA barcode, 28 chars}
#'   \item{filenames}{location of the idat files for each sample}
#'   \item{Cancer}{Primary Cancer Site}
#'   \item{Purity_Impute_RFPurify(Absolute)}{Imputed purity from the array}
#'   \item{Purity_ASCAT}{ Purity ref. value Purity predicted from SNPs array.}
#'   \item{ASCAT_Amp}{Amplified genes (5 to 10 copies) called by ASCAT}
#'   \item{ASCAT_Amp10}{Largely amplified genes (> 10 copies) called by ASCAT}
#'   \item{ASCAT_Gain}{Gains (3 & 4 copies) called by ASCAT}
#'   \item{ASCAT_HetLoss}{Missing 1 allele (1 copy) called by ASCAT}
#'   \item{ASCAT_HomDel}{Deleted gene in both alleles (0 copies)by ASCAT}
#'
#'
#' }
#' @source \url{https://github.com/riazn/biallelic_hr.git}
#' @usage data(TrainingSet_Sample_sheet)
#'
"TrainingSet_Sample_sheet"
