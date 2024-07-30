#' GC content and lengths of 2000 human genes
#'
#' GC content and lengths of 2000 human Ensembl genes. For each gene exons
#' are collapsed to get disjoint meta-exons and GC content is calculated
#' for each meta-exon. A exon length-weighted GC content is calculated for
#' each gene using the [calc_gene_gc()] function.
#'
#' @format
#' A data frame with 2000 rows and 2 columns:
#' \describe{
#'   \item{gc_content}{metagene-level GC content, values are between 0 and 1}
#'   \item{width}{length of metagenes in base pair(bp)}
#' }
#' @source Homo_sapiens.GRCh38.110.gtf.gz
#' @usage
#' data(gene_GC)
#'
"gene_GC"


#' GC content and lengths of 2000 intergenic regions
#'
#' GC content and lengths of 2000 human intergenic regions calculated using the
#' [calc_region_gc()] function.
#'
#' @format
#' A data frame with 2000 rows and 2 columns:
#' \describe{
#'   \item{gc_content}{metagene-level GC content, values are between 0 and 1}
#'   \item{width}{length of metagenes in base pair(bp)}
#' }
#' @source Homo_sapiens.GRCh38.110.gtf.gz
#' @usage
#' data(intergenic_GC)
"intergenic_GC"


#' GC content and lengths of 2000 intergenic regions
#'
#' GC content and lengths of 2000 human intergenic regions calculated using the
#' [calc_region_gc()] function.
#'
#' @format
#   A list containing the 7 or 8 (plant only) elements.
#'  "intergenic_region", "intronic_region", "rRNA", "mitochondrion", "gtf",
#'  "chloroplast" (plant only).
#'  For genomic features, gene, exon, intron, rRNA, mitochondrion, and
#'  chloroplast, only the `stat` elements of the [Rsubread::featureCounts()]
#'  output is stored.
#'  For intergenic region, the `stat`, `annotation` and `counts` elements
#'  are stored; for GTF-based summarization, the `stat` and `counts` elements
#'  are stored.
#' @usage
#' data(feature_counts_list)
#'
"feature_counts_list"

#' GC content and lengths of 2000 intergenic regions
#'
#' GC content and lengths of 2000 human intergenic regions calculated using the
#' [calc_region_gc()] function.
#'
#' @format
#' A list of three elements:
#' \describe{
#'   \item{abundance}{A numeric matrix containing abundance (TPM) for each gene
#'                    of each sample}
#'   \item{counts}{A numeric matrix containing read count (fraction) for each
#'                 gene of each sample}
#'   \item{length}{A numeric matrix containing length (bp) for each gene of
#'                 each sample}
#' }
#' @usage
#' data(salmon_quant)
"salmon_quant"
