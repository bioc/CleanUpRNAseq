#' @include check_gDNA.R
NULL

#' Global correction for DNA contamination
#'
#' Correct for DNA contamination in RNA-seq data using the 2.5 times of median
#' counts per base of intergenic regions with at least one count.
#'
#' @inheritParams plot_assignment_stat
#' @param lambda A positive number specifying how many times of the median read
#'   coverage of non-zero count intergenic regions used as an estimate of DNA
#'   contamination. The default of *lambda* is 1,but it could be adjusted based
#'   on the gene-level count distributions of the resulting corrected count
#'   table output by the [plot_read_distr()] function. A value between
#'   1 and 3 can be tried. Ideally the distributions of all samples from a
#'   given condition should be very similar.
#'
#' @return A matrix containing an RNA-seq count table corrected for DNA
#'   contamination, with rows for genes and columns for samples.
#' @export
#' @examples
#' lib_strand <- 0
#' col_data_f <- system.file("extdata", "example.colData.txt",
#'                          package = "CleanUpRNAseq")
#' col_data <- read.delim(col_data_f, as.is = TRUE)
#' ## create fake bam files
#' tmp_dir <- tempdir()
#' bamfiles <- gsub(".+/", "", col_data$BAM_file)
#' null <- lapply(file.path(tmp_dir, bamfiles), file.create)
#' ## create fake quant.sf files
#' quant_sf <- file.path(tmp_dir, gsub(".srt.bam$",
#'                                     "quant.sf",
#'                                     bamfiles))
#' null <- lapply(quant_sf, file.create)
#' col_data$BAM_file <- file.path(tmp_dir, bamfiles)
#' col_data$salmon_quant_file <- quant_sf
#'
#' ## pretend this is stranded RA=NA-seq data
#' col_data$salmon_quant_file_opposite_strand <- quant_sf
#'
#' sc <- create_summarizedcounts(lib_strand, col_data)
#'
#' data("feature_counts_list")
#' data("salmon_quant")
#'
#' sc$set_feature_counts(feature_counts_list)
#' sc$set_salmon_quant(salmon_quant)
#' sc$set_salmon_quant_opposite(salmon_quant)
#' corrected_counts <- correct_global(SummarizedCounts = sc)
#'
correct_global <- function(SummarizedCounts = NULL,
                              lambda = 1) {
    stopifnot(is(SummarizedCounts, "SummarizedCounts"))
    salmon_res <- SummarizedCounts$get_salmon_quant()

    if (lambda < 1) {
        stop("lambda should not be less than 1.")
    }

    ## intergenic read count per base
    intergenic_counts <-
        as.data.frame(SummarizedCounts$get_ir_counts())
    stopifnot(colnames(intergenic_counts) == colnames(salmon_res$counts))
    intergenic_length <-
        SummarizedCounts$get_ir_anno()[, 6, drop = FALSE]
    intergenic_fpb <- mapply(function(.x, .y) {
        .x / .y
    }, intergenic_counts, intergenic_length, SIMPLIFY = FALSE)

    intergenic_fpb <- as.data.frame(do.call("cbind", intergenic_fpb))
    colnames(intergenic_fpb) <- colnames(intergenic_counts)

    ## remove rows from intergenic_fpb, where all samples have 0
    intergenic_fpb <- intergenic_fpb[rowSums(intergenic_fpb) != 0, ]
    median_fpb <- vapply(intergenic_fpb, function(x) {
        10^(median(log10(x + 1))) - 1
    }, numeric(1))

    gene_len <- as.data.frame(salmon_res$length)
    dna_contamination_count <-
        do.call(
            cbind,
            mapply(function(len, cov) {
                len * cov
            }, gene_len, lambda * median_fpb, SIMPLIFY = FALSE)
        )

    colnames(dna_contamination_count) <- colnames(intergenic_counts)
    dna_contamination_count <-
        dna_contamination_count[, colnames(salmon_res$counts)]
    counts <- as.matrix(salmon_res$counts - dna_contamination_count)

    # set negative values to 0
    counts <- ifelse(counts < 0, 0, counts)
    counts <- round(counts)
    SummarizedCounts$set_global_correction(counts)
    counts
}

#' Correct for DNA contamination in consideration of GC-bias
#'
#'
#' @param intergenic_counts A data frame or matrix containing counts assigned to
#'   each intergenic region of each sample, such as the sublist, *counts*  from
#'   sublist named *intergenic_region* of an output from the [summarize_reads()]
#'   function.
#' @param intergenic_gc A data frame or matrix with two columns: gc_content in
#'   proportion between 0 and 1, and width in basepair, containing
#'   GC-content and lengths of each intergenic region, such as an output from
#'   the [calc_region_gc()] function.
#' @param plot A logical(1) vector, specifying whether to output a panel of
#'   scatter plot showing a bi-variate distribution of GC content and count per
#'   base of intergenic regions in each sample. Default is TRUE. A loess
#'   regression line is displayed in each subplot.
#'
#' @return A matrix containing estimated DNA contamination in the count per base
#'   of intergenic region binned by GC content, with a bin size of 5% of each
#'   sample, 20 bins. Column names are sample names, rownames are GC content
#'   bins in the form of '(0,0.05]'.
#'
#' @importFrom stats quantile loess predict
#' @import ggplot2
#' @noRd
gc_bin_contamination <- function(intergenic_counts,
                                 intergenic_gc,
                                 plot = TRUE) {
    estimate_contamination <- function(.count, .sample_name, .gc) {
        ## remove zero count intergenic entries
        .rm <- .count == 0
        .count <- .count[!.rm]
        .gc <- .gc[!.rm, ]

        count_per_base <- .count / .gc$width

        # remove first and last quantile
        keep <-
            # count_per_base >= quantile(count_per_base, probs = 0.01) &
            count_per_base <= quantile(count_per_base, probs = 0.99)

        .gc <- .gc[keep, ]
        count_per_base <- count_per_base[keep]

        ## loess fitting
        l <- loess(count_per_base ~ .gc$gc_content,
                   family = "symmetric",
                   span = 0.25
        )

        data <- cbind(
            .gc, count_per_base,
            predict(l, newdata = .gc$gc_content)
        )
        colnames(data) <- c("gc_content", "width", "CPB", "predict")
        data$sample_name <- .sample_name

        temp_data <- data[order(data$gc_content), ]
        temp_data <- split(temp_data, f = cut(temp_data$gc_content,
                                              breaks = seq(0, 1, 0.05)
        ))
        perc5_bin_cpb <- vapply(
            temp_data,
            function(.x) {
                quantile(.x$CPB, prob = 0.5)
            },
            numeric(1)
        )
        perc5_bin_cpb[10:20] <- c(perc5_bin_cpb[9:1], 0, 0)
        list(contamination = perc5_bin_cpb, plot_data = data)
    }

    contamination <- mapply(
        function(.intergenic_count,
                 .sample_name,
                 .intergenic_gc) {
            contamination_level <- estimate_contamination(.intergenic_count,
                                                          .sample_name,
                                                          .gc = .intergenic_gc
            )
        },
        intergenic_counts,
        colnames(intergenic_counts),
        MoreArgs = list(.intergenic_gc = intergenic_gc),
        SIMPLIFY = FALSE
    )
    plot_data <- do.call("rbind", lapply(contamination, "[[", 2))
    plot_data$sample_name <- factor(plot_data$sample_name,
                                    levels = unique(plot_data$sample_name)
    )
    gc_bin_contamination <-
        do.call(cbind, lapply(contamination, "[[", 1))
    rownames(gc_bin_contamination) <- gsub(
        "\\.\\d+%", "",
        rownames(gc_bin_contamination)
    )

    p <- ggplot(plot_data, aes(x = gc_content, y = CPB)) +
        geom_point(
            color = "gray",
            alpha = 0.3,
            size = 0.3
        ) +
        geom_line(aes(x = gc_content, y = predict),
                  color = "blue"
        ) +
        xlab("GC%") +
        ylab("Count per base") +
        facet_wrap(. ~ sample_name, ncol = 3)
    if (plot) {
        plot(p)
    }
    gc_bin_contamination
}


#' Correct DNA contamination considering GC-bias effect
#'
#' Correct DNA contamination considering GC-bias effect on fragment
#' amplification. Intergenic regions are binned based on their GC content
#' ranging from 0 to 100%, with a bin size of 5%. Per gene DNA contamination
#' is estimated as the product of count per base in a GC content matching bin
#' of intergenic regions and the total collapsed exons of a gene and is
#' subtracted away from the gene count matrix.
#'
#'
#' @inheritParams plot_assignment_stat
#' @param gene_gc A data frame or matrix with two columns: gc_content in
#'   proportion between 0 and 1, and width in basepair, containing
#'   GC-content and total exon lengths of each gene. An output of the
#'   [calc_gene_gc()] function.
#' @param intergenic_gc A data frame or matrix with two columns: gc_content in
#'   proportion between 0 and 1, and width in basepair, containing
#'   GC-content and lengths of each intergenic region, such as an output from
#'   the [calc_region_gc()] function.
#' @param plot A logical(1) vector, specifying whether to output a panel of
#'   scatter plot showing a bi-variate distribution of GC content and count per
#'   base of intergenic regions in each sample. Default is TRUE.
#' @return A data frame containing corrected count for each gene (row) of each
#'   sample (column).
#' @export
#'
#' @examples
#' lib_strand <- 0
#' col_data_f <- system.file("extdata", "example.colData.txt",
#'                          package = "CleanUpRNAseq")
#' col_data <- read.delim(col_data_f, as.is = TRUE)
#' ## create fake bam files
#' tmp_dir <- tempdir()
#' bamfiles <- gsub(".+/", "", col_data$BAM_file)
#' null <- lapply(file.path(tmp_dir, bamfiles), file.create)
#' ## create fake quant.sf files
#' quant_sf <- file.path(tmp_dir, gsub(".srt.bam$",
#'                                     "quant.sf",
#'                                     bamfiles))
#' null <- lapply(quant_sf, file.create)
#' col_data$BAM_file <- file.path(tmp_dir, bamfiles)
#' col_data$salmon_quant_file <- quant_sf
#'
#' ## pretend this is stranded RA=NA-seq data
#' col_data$salmon_quant_file_opposite_strand <- quant_sf
#'
#' sc <- create_summarizedcounts(lib_strand, col_data)
#'
#' data("feature_counts_list")
#' data("salmon_quant")
#'
#' sc$set_feature_counts(feature_counts_list)
#' sc$set_salmon_quant(salmon_quant)
#' sc$set_salmon_quant_opposite(salmon_quant)
#'
#' data("gene_GC")
#' data("intergenic_GC")
#' gc_bias_corrected_counts <-
#'     correct_GC(
#'         SummarizedCounts = sc,
#'         gene_gc = gene_GC,
#'         intergenic_gc = intergenic_GC,
#'         plot = FALSE
#'     )
#'
correct_GC <- function(SummarizedCounts = NULL,
                       gene_gc = NULL,
                       intergenic_gc = NULL,
                       plot = FALSE) {
    stopifnot(is(SummarizedCounts, "SummarizedCounts"))
    salmon_res <- SummarizedCounts$get_salmon_quant()
    intergenic_counts <- SummarizedCounts$get_ir_counts()

    stopifnot(nrow(intergenic_counts) > 1, nrow(salmon_res) > 1)

    if (!is.data.frame(gene_gc) && !is.matrix(gene_gc)) {
        stop("intergenic_counts must be a data.frame or matrix.")
    }

    if (!is.data.frame(intergenic_gc) && !is.matrix(intergenic_gc)) {
        stop("intergenic_counts must be a data.frame or matrix.")
    }

    if (!setequal(rownames(salmon_res$counts), rownames(gene_gc))) {
        warning(
            "Gene names in the Salmon count table is not exactly the same as\n",
            "those of the gene_gc data frame. The intersection of gene ",
            "names from both\n",
            "sources will be used."
        )

        common_genes <-
            intersect(rownames(salmon_res$counts), rownames(gene_gc))
        if (length(common_genes) < 10) {
            stop(
                "There are only",
                length(common_genes),
                "genes shared between the ",
                "Salmon count table and the gene_gc data frame. Please double ",
                "check the rownames."
            )
        }
        salmon_res <- lapply(salmon_res[seq_len(3)], function(.x) {
            .x <- .x[common_genes, ]
        })

        gene_gc <- gene_gc[common_genes, ]
    }

    ## intergenic gc and counts
    if (any(!rownames(intergenic_counts) %in% rownames(intergenic_gc)) ||
        any(!rownames(intergenic_gc) %in% rownames(intergenic_counts))) {
        stop("rownames of intergenic_gc DO NOT match those of",
             " intergenic_counts!")
    }

    intergenic_counts <-
        as.data.frame(intergenic_counts[rownames(intergenic_gc), ])

    gc_bin_cpb <- gc_bin_contamination(intergenic_counts,
                                       intergenic_gc,
                                       plot = plot
    )

    gc_bin_cpb <- as.data.frame(gc_bin_cpb)
    gc_bin_cpb$bin <- rownames(gc_bin_cpb)
    rownames(gc_bin_cpb) <- NULL

    gene_gc <- gene_gc[order(gene_gc$gc_content), ]
    gene_gc$bin <- cut(gene_gc$gc_content, breaks = seq(0, 1, 0.05))
    gene_gc$width <- NULL
    gene_gc$GeneID <- rownames(gene_gc)

    gene_gc_sample <- merge(gene_gc, gc_bin_cpb,
                            by = "bin",
                            all.x = TRUE,
                            sort = FALSE
    )
    rownames(gene_gc_sample) <- gene_gc_sample$GeneID
    gene_gc_sample$GeneID <- NULL
    gene_gc_sample <- gene_gc_sample[rownames(salmon_res$length), ]
    gene_contamination <-
        gene_gc_sample[, -c(1, 2)] * salmon_res$length
    corrected_counts <- salmon_res$counts - gene_contamination
    corrected_counts[corrected_counts < 0] <- 0
    corrected_counts <- as.matrix(corrected_counts)
    SummarizedCounts$set_gc_correction(corrected_counts)
    corrected_counts
}


#' Correct for gDNA contamination in stranded libraries
#'
#' Correct for gDNA contamination in stranded libraries based on Salmon
#' quantitation using the real and opposite library strandedness information
#'
#' @inheritParams plot_assignment_stat
#'
#' @return A list of of matrices of gene-level abundance, counts, and length.
#'   See [tximport::tximport()]. The count matrix is corrected for DNA
#'   contamination and rounded into integers.
#'   \describe{
#'   \item{abundance}{A numeric matrix containing corrected abundance (TPM) for
#'                    each gene of each sample}
#'   \item{counts}{An *integer* matrix containing *corrected* read count for
#'                 each gene of each sample}
#'   \item{length}{A numeric matrix containing length (bp) for each gene of
#'                 each sample}
#' }
#' @export
#' @importFrom tximport tximport
#' @importFrom ensembldb transcripts
#'
#' @examples
#' lib_strand <- 1
#' col_data_f <- system.file("extdata", "example.colData.txt",
#'                          package = "CleanUpRNAseq")
#' col_data <- read.delim(col_data_f, as.is = TRUE)
#' ## create fake bam files
#' tmp_dir <- tempdir()
#' bamfiles <- gsub(".+/", "", col_data$BAM_file)
#' null <- lapply(file.path(tmp_dir, bamfiles), file.create)
#' ## create fake quant.sf files
#' quant_sf <- file.path(tmp_dir, gsub(".srt.bam$",
#'                                     "quant.sf",
#'                                     bamfiles))
#' null <- lapply(quant_sf, file.create)
#' col_data$BAM_file <- file.path(tmp_dir, bamfiles)
#' col_data$salmon_quant_file <- quant_sf
#'
#' ## pretend this is stranded RA=NA-seq data
#' col_data$salmon_quant_file_opposite_strand <- quant_sf
#'
#' sc <- create_summarizedcounts(lib_strand, col_data)
#'
#' data("feature_counts_list")
#' data("salmon_quant")
#'
#' sc$set_feature_counts(feature_counts_list)
#' sc$set_salmon_quant(salmon_quant)
#' sc$set_salmon_quant_opposite(salmon_quant)
#' stranded_correction <- correct_stranded(SummarizedCounts = sc)
#'

correct_stranded <-
    function(SummarizedCounts = NULL) {
        stopifnot(is(SummarizedCounts, "SummarizedCounts"))

        lib_strand <- SummarizedCounts$lib_strand
        stopifnot(lib_strand %in% c(1,2))

        salmon_strand <- SummarizedCounts$get_salmon_quant()
        salmon_opposite_strand <- SummarizedCounts$get_salmon_quant_opposite()
        stopifnot(nrow(salmon_strand$counts) > 1,
                  nrow(salmon_opposite_strand$counts) > 1)

        corrected_salmon <- mapply(
            function(.strand, .reverse_strand) {
                val <- .strand - .reverse_strand
                val[val < 0] <- 0
                val
            },
            salmon_strand[c("abundance", "counts")],
            salmon_opposite_strand[c("abundance", "counts")],
            SIMPLIFY = FALSE
        )

        corrected_salmon$length <- salmon_strand$length
        SummarizedCounts$set_stranded_correction(corrected_salmon)
        corrected_salmon
    }

#' Correct gene expression using a linear model
#'
#' Correct gene expression with IR% as a co-variate in linear model in the
#' limma framework.
#'
#' @inheritParams plot_assignment_stat
#' @param design A design matrix with rows corresponding to samples and
#'   columns to coefficients to be estimated. See [limma::lmFit()] and
#'   Law et al. 2020, F1000Research, doi: 10.12688/f1000research.27893.1.
#'
#' @importFrom limma voom lmFit
#'
#' @return A numeric matrix of normalized expression values on the log2 scale,
#'   corrected for gDNA contamination, and `batch` (if any).
#' @export
#'
#' @examples
#' lib_strand <- 0
#' col_data_f <- system.file("extdata", "example.colData.txt",
#'                          package = "CleanUpRNAseq")
#' col_data <- read.delim(col_data_f, as.is = TRUE)
#' ## create fake bam files
#' tmp_dir <- tempdir()
#' bamfiles <- gsub(".+/", "", col_data$BAM_file)
#' null <- lapply(file.path(tmp_dir, bamfiles), file.create)
#' ## create fake quant.sf files
#' quant_sf <- file.path(tmp_dir, gsub(".srt.bam$",
#'                                     "quant.sf",
#'                                     bamfiles))
#' null <- lapply(quant_sf, file.create)
#' col_data$BAM_file <- file.path(tmp_dir, bamfiles)
#' col_data$salmon_quant_file <- quant_sf
#'
#' ## pretend this is stranded RA=NA-seq data
#' col_data$salmon_quant_file_opposite_strand <- quant_sf
#'
#' sc <- create_summarizedcounts(lib_strand, col_data)
#'
#' data("feature_counts_list")
#' data("salmon_quant")
#'
#' sc$set_feature_counts(feature_counts_list)
#' sc$set_salmon_quant(salmon_quant)
#' sc$set_salmon_quant_opposite(salmon_quant)
#' assigned_per_region <- get_region_stat(SummarizedCounts = sc)
#' design = model.matrix(~ group + batch +IR_rate, data = sc$col_data)
#' ir_corrected <- correct_IR(sc, design)
#'
#'
correct_IR <- function(SummarizedCounts = NULL,
                       design = NULL) {
    stopifnot(is(SummarizedCounts, "SummarizedCounts"))
    counts <- SummarizedCounts$get_salmon_counts()
    stopifnot(nrow(counts) > 1)

    col_data <- SummarizedCounts$col_data

    ## adjusted expression using voom()
    v <- voom(counts, design = design)
    fit <- lmFit(v, design = design)
    group_num <- nlevels(as.factor(col_data$group))
    if (ncol(design) > group_num) {
        col_covariate <- seq_len(ncol(design))[-seq_len(group_num)]
        adj_exp <- v$E - fit$coefficients[, col_covariate] %*%
            t(design[, col_covariate])
    } else {
        adj_exp <- v$E
    }
    SummarizedCounts$set_ir_correction(as.matrix(adj_exp))
    adj_exp
}
