#' Sample sheet for the 442 TCGA samples used as training set.
#'
#' A dataset containing the IDs & paths to the 450k microarrays used
#' as training set to infer the reference KCN as described in the paper.
#' Reference ASCAT data needed to estimate each KCN is extracted from a previous
#' work https://github.com/riazn/biallelic_hr.git.
#'
#' @format A data frame with 442 rows and 10 variables:
#' \describe{
#'   \item{Project}{Cancer project name ex: TCGA-cancer}
#'   \item{barcodes}{Unique identifier of the sample. TCGA barcode, 12 chars}
#'   \item{ASCAT.path}{Path to origin file with CN calls}
#'   \item{gene_name}{Copy number state called by ASCAT for gene_name}
#' }
#' @source \url{https://github.com/riazn/biallelic_hr.git}
#' @usage data(all_CN_ASCAT)
#'
"TrainingSet_Sample_sheet"
