% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_gDNA.R
\name{plot_sample_corr}
\alias{plot_sample_corr}
\title{Visualize sample correlation}
\usage{
plot_sample_corr(SummarizedCounts = NULL)
}
\arguments{
\item{SummarizedCounts}{An object of \link{SummarizedCounts}..}
}
\value{
NULL. A plot with pairwise scatter plots and Pearson correlation
coefficients is generated. When the sample size is big, save it as tiff
file.
}
\description{
Generate a panel of plots based on a count table, with the diagonal
showing the sample names, the lower triangle showing smoothed scatterplots
for gene expression of pairwise samples, and the upper triangle showing
Pearson correlation coefficients of gene expression of pairwise samples.
}
\examples{
lib_strand <- 0
col_data_f <- system.file("extdata", "example.colData.txt",
                         package = "CleanUpRNAseq")
col_data <- read.delim(col_data_f, as.is = TRUE)
## create fake bam files
tmp_dir <- tempdir()
bamfiles <- gsub(".+/", "", col_data$BAM_file)
null <- lapply(file.path(tmp_dir, bamfiles), file.create)
## create fake quant.sf files
quant_sf <- file.path(tmp_dir, gsub(".srt.bam$",
                                    "quant.sf",
                                    bamfiles))
null <- lapply(quant_sf, file.create)
col_data$BAM_file <- file.path(tmp_dir, bamfiles)
col_data$salmon_quant_file <- quant_sf

## pretend this is stranded RA=NA-seq data
col_data$salmon_quant_file_opposite_strand <- quant_sf

sc <- create_summarizedcounts(lib_strand, col_data)

data("feature_counts_list")
data("salmon_quant")

sc$set_feature_counts(feature_counts_list)
sc$set_salmon_quant(salmon_quant)
sc$set_salmon_quant_opposite(salmon_quant)
plot_sample_corr(SummarizedCounts = sc)

}
