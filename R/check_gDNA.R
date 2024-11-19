#' @include R6_class.R  misc.R
NULL


#' Visualize assignment statistics of reads/fragments by featureCounts
#'
#' @param SummarizedCounts An object of [SummarizedCounts]..
#'
#' @return A ggplot object, showing percentages and number of fragments in each
#'   assignment category as determined by [Rsubread::featureCounts()] based on
#'   a GTF.
#'
#' @export
#' @importFrom reshape2 melt
#' @import ggplot2
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
#' plot_assignment_stat(SummarizedCounts = sc)
#'
#'
plot_assignment_stat <- function(SummarizedCounts = NULL) {
    stopifnot(is(SummarizedCounts, "SummarizedCounts"))
    assignment_stat <- SummarizedCounts$get_gtf_stat()
    stopifnot(is.data.frame(assignment_stat),
              nrow(assignment_stat) >= 1)
    assignment_stat <-
        assignment_stat[rowSums(assignment_stat[, -1]) != 0, ]
    assignment_stat_pct <- assignment_stat
    assignment_stat_pct[, seq_len(ncol(assignment_stat_pct))[-1]] <-
        mapply(
            function(.x, .y) {
                .x / .y * 100
            },
            assignment_stat[, -1],
            colSums(assignment_stat[, -1]),
            SIMPLIFY = FALSE
        )

    assignment_stat[, seq_len(ncol(assignment_stat_pct))[-1]] <-
        assignment_stat[, seq_len(ncol(assignment_stat_pct))[-1]] / 10^6

    assignment_count_long <- melt(
        assignment_stat,
        id.vars = "Status",
        variable.name = "Sample",
        value.name = "Statistics"
    )
    assignment_count_long$Stat <- "Count"
    assignment_pct_long <- melt(
        assignment_stat_pct,
        id.vars = "Status",
        variable.name = "Sample",
        value.name = "Statistics"
    )
    assignment_pct_long$Stat <- "Percentage"
    assignment_long <-
        rbind(assignment_count_long, assignment_pct_long)

    assignment_long$Sample <- factor(assignment_long$Sample,
                                     levels = unique(assignment_long$Sample)
    )

    p <- ggplot(assignment_long, aes(
        x = Sample, y = Statistics,
        fill = Status
    )) +
        geom_bar(stat = "identity", color = "white") +
        ylab("Reads") +
        facet_wrap(~Stat, scales = "free_y") +
        theme(
            text = element_text(size = 8),
            axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0,
                size = 8
            ),
            legend.key.size = unit(0.5, "cm"),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6)
        )
    p
}

#' Calculate read distribution over different types of genomic features
#'
#' Calculate read distribution over different types of genomic features: genes,
#' exons, introns, intergenic regions,rRNA regions, and organelle genome(s).
#' @inheritParams plot_assignment_stat
#'
#' @return A data frame described as below.
#'
#' @export
#' @importFrom graphics pairs par smoothScatter strwidth text
#' @importFrom grDevices dev.off pdf
#' @importFrom stats cor dist median prcomp
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
#' assigned_per_region <- get_region_stat(SummarizedCounts = sc)
#' assigned_per_region

get_region_stat <- function(SummarizedCounts = NULL) {
    stopifnot(is(SummarizedCounts, "SummarizedCounts"))
    featurecounts_list <- SummarizedCounts$get_feature_counts()
    col_data <- SummarizedCounts$col_data

    if (is.null(featurecounts_list)) {
        stop("feature_counts field has not been populated")
    }

    kept <- names(featurecounts_list) != "gtf"
    assigned_per_region <- mapply(
        function(.x, .name) {
            if (.name %in% c("gtf", "intergenic_region"))
            {
                assigned <- .x$stat[.x$stat$Status == "Assigned", ,
                                    drop = FALSE]
            } else {
                assigned <- .x[.x$Status == "Assigned", , drop = FALSE]
            }
            assigned <- assigned[, colnames(assigned) != "Status",
                                 drop = FALSE]
            assigned <- as.data.frame(t(assigned))
            colnames(assigned) <- "assigned_count"
            assigned$sample_name <- rownames(assigned)
            rownames(assigned) <- NULL
            assigned$region_type <- .name
            assigned
        },
        featurecounts_list[kept],
        names(featurecounts_list)[kept],
        SIMPLIFY = FALSE
    )
    assigned_per_region <- do.call("rbind", assigned_per_region)
    gtf_stat <- featurecounts_list$gtf$stat
    gtf_stat <- gtf_stat[!grepl("Unmapped", gtf_stat$Status,
                                ignore.case = TRUE),
                         colnames(gtf_stat) != "Status",
                         drop = FALSE]
    total_mapped_frags <- colSums(gtf_stat)
    total_mapped_frags <- rep(total_mapped_frags, sum(kept))
    assigned_per_region$assigned_percent <-
        assigned_per_region$assigned_count / total_mapped_frags * 100

    levels <- names(featurecounts_list)[kept]
    labels <- gsub("_", " ", levels)
    labels <- paste0(toupper(gsub("(^.).+", "\\1", labels)),
                     gsub("^.(.+)", "\\1", labels))
    labels <- ifelse(labels == "RRNA", "rRNA", labels)
    assigned_per_region$region_type <-
        factor(assigned_per_region$region_type,
               levels = levels,
               labels = labels
        )
    if (!all(sort(as.character(col_data$sample_name)) ==
             sort(unique(assigned_per_region$sample_name)))) {
        stop(
            "sample names in col_data are not consistent with those in ",
            "featurecounts_list"
        )
    }

    assigned_per_region <- merge(assigned_per_region,
                                 col_data[, c("sample_name", "group")],
                                 by = "sample_name",
                                 all.x = TRUE,
                                 sort = FALSE
    )

    assigned_per_region <-
        assigned_per_region[order(
            assigned_per_region$group,
            assigned_per_region$sample_name
        ), ]
    assigned_per_region$group <- factor(assigned_per_region$group,
                                        levels = unique(assigned_per_region$group)
    )
    assigned_per_region$sample_name <-
        factor(assigned_per_region$sample_name,
               levels = unique(assigned_per_region$sample_name)
        )
    IR_rate <- get_IR_rate(assigned_per_region = assigned_per_region)
    SummarizedCounts$add_ir_rate(IR_rate)
    assigned_per_region
}



#' Visualize read distribution among different genomic regions
#'
#' Generate ggplot plots showing percentages of fragments assigned to different
#' type of genomic features: genes, exons, introns, intergenic regions,
#' rRNA regions, and organelle genome(s).
#'
#' @param assigned_per_region A data frame, outputted by [get_region_stat()].
#'
#' @return A ggplot object, showing percentages of fragments assigned to
#'   different genomic features, such as genic regions, intergenic regions,
#'   exonic regions, intronic regions, rRNA genes, mitochondrial genome,
#'   chloroplast genome (only for plants)
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
#' assigned_per_region <- get_region_stat(SummarizedCounts = sc)
#' p <- plot_read_distr(assigned_per_region)
#' p
#'
#'
plot_read_distr <- function(assigned_per_region = NULL) {
    p <- ggplot(
        assigned_per_region,
        aes(
            x = sample_name,
            y = assigned_percent,
            color = group
        )
    ) +
        geom_point() +
        xlab("Sample") +
        ylab("Assigned reads (%)") +
        facet_wrap(~region_type) +
        guides(color = guide_legend(title = "Group")) +
        theme(
            text = element_text(size = 8),
            axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0,
                size = 8
            )
        )
}

#' Access IR rate
#'
#' Get the percentages of reads mapping to intergenic region per sample.
#'
#' @inheritParams plot_read_distr
#'
#' @return A data frame containing the percentages of reads mapping to
#'   intergenic regions, which can be used for the "IR%" method for correcting
#'   gDNA contamination.
#'
#' @noRd
#'
#'
get_IR_rate <- function(assigned_per_region = NULL) {
    IR_rate <-
        assigned_per_region[assigned_per_region$region_type ==
                                "Intergenic region", ]
    rownames(IR_rate) <- IR_rate$sample_name
    IR_rate <- IR_rate[, "assigned_percent", drop = FALSE]
    colnames(IR_rate) <- "IR_rate"
    IR_rate
}


#' Visualize sample correlation
#'
#' Generate a panel of plots based on a count table, with the diagonal
#' showing the sample names, the lower triangle showing smoothed scatterplots
#' for gene expression of pairwise samples, and the upper triangle showing
#' Pearson correlation coefficients of gene expression of pairwise samples.
#'
#' @inheritParams plot_assignment_stat
#' @return NULL. A plot with pairwise scatter plots and Pearson correlation
#'   coefficients is generated. When the sample size is big, save it as tiff
#'   file.
#'
#' @export
#' @importFrom methods is
#' @importFrom graphics smoothScatter
#' @importFrom KernSmooth bkde2D
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
#' plot_sample_corr(SummarizedCounts = sc)
#'

plot_sample_corr <- function(SummarizedCounts = NULL) {
    stopifnot(is(SummarizedCounts, "SummarizedCounts"))
    counts <- SummarizedCounts$get_salmon_counts()
    stopifnot(nrow(counts) > 1)

    ## scatter plot and correlation plot
    counts <- counts[rowSums(counts) > 0, ]
    counts <- log2(counts + 1)
    counts <- counts[, order(colnames(counts))]

    panel.cor <-
        function(x,
                 y,
                 digits = 2,
                 prefix = "",
                 cex.cor,
                 ...) {
            usr <- par("usr")
            on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            Cor <- cor(x, y) # Remove abs function if desired
            txt <-
                paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
            if (missing(cex.cor)) {
                cex.cor <- 0.4 / strwidth(txt)
            }
            text(0.5, 0.5, txt,
                 cex = 0.5 + cex.cor * Cor
            )
        }

    panel.smth <- function(x, y) {
        smoothScatter(x, y, add = TRUE)
    }

    pairs(
        counts,
        lower.panel = panel.smth,
        upper.panel = panel.cor,
        font.labels = 1,
        row1attop = TRUE,
        gap = 0.2,
        log = ""
    )
}


#' Visualize expression distribution
#'
#' Compare expression distribution by boxplot, density plot and empirical
#' cumulative distribution plot.
#'
#' @inheritParams plot_assignment_stat
#' @param normalization A character(1) vector, specifying a between-sample
#'   normalization methods: DESeq2's  median of ratios method, smooth quantile
#'   normalization method [qsmooth::qsmooth()], or `none`.
#' @return A list of 3 ggplot objects.
#' \describe{
#'   \item{box_plot}{boxplots showing DESeq2-normalized gene-level count
#'                  distribution, on a log scale}
#'   \item{density_plot}{density plots showing DESeq2-normalized gene-level
#'                       count distribution, on a log scale}
#'   \item{ecd_plot}{plots showing empricial cumulative distribution of
#'                  fraction  of genes with CPM greater than or equal to a
#'                  given CPM on a log scale}
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom edgeR filterByExpr cpm
#' @importFrom qsmooth qsmooth qsmoothData
#' @importFrom reshape2 melt
#' @examplesIf require("patchwork")
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
#' wrap_plots(plot_expr_distr(SummarizedCounts = sc,
#'                 normalization = "DESeq2"), ncol = 1)

plot_expr_distr <-
    function(SummarizedCounts = NULL,
             normalization = c("DESeq2", "qsmooth", "none")) {
        stopifnot(is(SummarizedCounts, "SummarizedCounts"))
        counts <- SummarizedCounts$get_salmon_counts()
        col_data <- SummarizedCounts$col_data

        stopifnot(nrow(counts) > 1)
        if (!is.matrix(counts)) {
            counts <- as.matrix(counts)
        }
        normalization <- match.arg(normalization)

        if (!setequal(colnames(counts), as.character(col_data$sample_name))) {
            stop(
                "Column names of the raw count matrix DO NOT match ",
                "sample names in the col_data"
            )
        } else {
            counts <- counts[, col_data$sample_name]
        }
        col_data <-
            col_data[order(col_data$group, col_data$sample_name), ]
        col_data$group <- factor(col_data$group,
                                 levels = unique(col_data$group)
        )
        col_data$sample_name <- factor(col_data$sample_name,
                                       levels = col_data$sample_name
        )

        keep <- filterByExpr(counts, group = col_data$group)
        counts <- counts[keep, ]

        if (normalization == "qsmooth") {
            # library size normalization
            lib_sizes <- colSums(counts)
            geometric_mean <- exp(mean(log(lib_sizes)))
            size_factors <- lib_sizes / geometric_mean
            counts <- sweep(counts, 2, size_factors, FUN = "/")
            counts_qs <- qsmooth(
                object = counts,
                group_factor = col_data$group
            )
            counts_1 <- qsmoothData(counts_qs)
        } else if (normalization == "DESeq2") {
            ## use DESeq2 normalization method instead
            dds <-
                DESeqDataSetFromMatrix(
                    countData = as.matrix(round(counts)),
                    colData = col_data,
                    design = ~ 0 + group
                )
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersions(dds,
                                       fitType = "local",
                                       maxit = 1000
            )

            counts_1 <- counts(dds, normalized = TRUE)
        } else {
            counts_1 <- counts
        }

        counts_1_log <- as.data.frame(log2(counts_1 + 1))
        counts_1_log$GeneID <- rownames(counts_1_log)
        counts_1_log_long <- melt(
            counts_1_log,
            id.vars = "GeneID",
            variable.name = "Sample",
            value.name = "log2Count"
        )
        counts_1_log_long <- merge(
            counts_1_log_long,
            col_data[, c("sample_name", "group")],
            by.x = "Sample",
            by.y = "sample_name",
            all.x = TRUE,
            sort = FALSE
        )
        counts_1_log_long <-
            counts_1_log_long[order(
                counts_1_log_long$group,
                counts_1_log_long$Sample
            ), ]
        counts_1_log_long$group <-
            factor(counts_1_log_long$group,
                   levels = unique(counts_1_log_long$group)
            )
        counts_1_log_long$Sample <-
            factor(counts_1_log_long$Sample,
                   levels = unique(counts_1_log_long$Sample)
            )
        p_box <- ggplot(
            counts_1_log_long,
            aes(
                x = Sample,
                y = log2Count,
                fill = group
            )
        ) +
            geom_boxplot() +
            ylab(expression(log[2](counts + 1))) +
            xlab("Sample") +
            guides(fill = guide_legend(title = "Group")) +
            ggtitle("Boxplot-normalized gene expression") +
            theme(
                text = element_text(size = 8),
                axis.text.x = element_text(
                    angle = 90,
                    vjust = 0.5,
                    hjust = 0,
                    size = 8
                ),
                plot.title = element_text(hjust = 0.5),
                legend.key.size = unit(0.5, "cm"),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 6)
            )
        if (length(unique(col_data$sample_name)) > 12) {
            p_density <- ggplot(
                counts_1_log_long,
                aes(
                    x = log2Count,
                    group = Sample,
                    color = group
                )
            ) +
                guides(color = guide_legend(title = "Group"))
        } else {
            p_density <- ggplot(
                counts_1_log_long,
                aes(
                    x = log2Count,
                    group = Sample,
                    color = Sample
                )
            ) +
                guides(color = guide_legend(title = "Sample"))
        }

        p_density <- p_density +
            geom_density(linewidth = 0.8) +
            ylab("Density") +
            xlab(expression(log[2](counts + 1))) +
            ggtitle("Density plot-normalized counts") +
            theme(
                text = element_text(size = 8),
                legend.key.size = unit(0.5, "cm"),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 6)
            )

        cpms_0 <- as.data.frame(cpm(counts))
        cpms_log <- log2(cpms_0 + 1 / 10 * min(cpms_0[cpms_0 != 0]))
        cpms_log$GeneID <- rownames(cpms_log)
        cpms_log_long <- melt(
            cpms_log,
            id.vars = "GeneID",
            variable.name = "Sample",
            value.name = "log2CPM"
        )
        cpms_log_long <- merge(
            cpms_log_long,
            col_data[, c("sample_name", "group")],
            by.x = "Sample",
            by.y = "sample_name",
            all.x = TRUE,
            sort = FALSE
        )
        cpms_log_long <- cpms_log_long[order(
            cpms_log_long$group,
            cpms_log_long$Sample
        ), ]
        cpms_log_long$group <- factor(cpms_log_long$group,
                                      levels = unique(cpms_log_long$group)
        )
        cpms_log_long$Sample <- factor(cpms_log_long$Sample,
                                       levels = unique(cpms_log_long$Sample)
        )

        if (length(unique(col_data$sample_name)) > 12) {
            p_ecdf <- ggplot(
                cpms_log_long,
                aes(
                    x = log2CPM,
                    group = Sample,
                    colour = group
                )
            ) +
                guides(color = guide_legend(title = "Group"))
        } else {
            p_ecdf <- ggplot(
                cpms_log_long,
                aes(
                    x = log2CPM,
                    group = Sample,
                    colour = Sample
                )
            ) +
                guides(color = guide_legend(title = "Sample"))
        }
        p_ecdf <- p_ecdf +
            stat_ecdf(geom = "step", linewidth = 0.8) +
            xlab(bquote(log[2] * CPM)) +
            ylab("Proportion") +
            ggtitle("Empirical cumulative distributions") +
            theme(
                text = element_text(size = 8),
                legend.key.size = unit(0.5, "cm"),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 6)
            )

        list(
            box_pot = p_box,
            density_plot = p_density,
            ecd_plot = p_ecdf
        )
    }

#' Check the percentage of genes with counts greater than minimal CPM
#'
#' @inheritParams plot_assignment_stat
#' @param min_cpm A numeric(1), minimal CPM threshold.
#' @param min_tpm A numeric(1), minimal TPM threshold.
#' @return A ggplot object if `counts` is not `NULL`, showing percentages of
#'   genes with counts above the user-specified minimal CPM (count per million)
#'   in each sample. Or a ggplot object of two panels if both `counts` and
#'   `abundance` are not `NULL`, showing percentages of genes  with counts
#'   above the user-specified minimal CPM and minimal TPM (transcript per
#'   million) in each sample.
#' @details
#'  The axis title contains unicode, so please output the plot in the svg format
#'  using the svglite package, instead of [grDevices::pdf()], for high-
#'  resolution plot.
#' @return A ggplot object.
#'
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
#' plot_gene_content(
#'     SummarizedCounts = sc,
#'     min_cpm = 1,
#'     min_tpm =1
#' )
#'

plot_gene_content <-
    function(SummarizedCounts = NULL,
             min_cpm = 1,
             min_tpm = 1) {
        stopifnot(is(SummarizedCounts, "SummarizedCounts"))
        counts <- SummarizedCounts$get_salmon_counts()
        abundance <- SummarizedCounts$get_salmon_abundance()
        col_data <- SummarizedCounts$col_data

        stopifnot(nrow(counts) > 1, nrow(abundance) > 1)

        if (is.data.frame(counts)) {
            counts <- as.matrix(counts)
        }

        if (!is.numeric(min_cpm) || min_cpm <= 0) {
            stop("min_cpm must be a single positive number")
        }
        if (!is.numeric(min_tpm) || min_tpm <= 0) {
            stop("min_tpm must be a single positive number")
        }

        ## filter all zero genes
        # counts <- counts[rowSums(counts)> 0, ]
        if ((!is.null(counts) &&
             any(!colnames(counts) %in% col_data$sample_name)) ||
            (!is.null(abundance) &&
             any(!colnames(abundance) %in% col_data$sample_name))) {
            stop("colnames of expression data DON'T match those in col_data")
        }

        if (!is.null(counts)) {
            cpms <- as.data.frame(cpm(counts))
            pct_cpm_gt1 <- vapply(
                cpms,
                function(.x) {
                    sum(.x >= min_cpm)
                },
                numeric(1)
            ) / nrow(counts) * 100

            pct_cpm_gt1_df <- data.frame(
                sample_name = names(pct_cpm_gt1),
                percent = pct_cpm_gt1
            )
            pct_cpm_gt1_df <- merge(pct_cpm_gt1_df,
                                    col_data,
                                    by = "sample_name",
                                    all.x = TRUE,
                                    sort = FALSE
            )
            pct_cpm_gt1_df <- pct_cpm_gt1_df[order(
                pct_cpm_gt1_df$group,
                pct_cpm_gt1_df$sample
            ), ]
            pct_cpm_gt1_df$group <- factor(pct_cpm_gt1_df$group,
                                           levels = unique(pct_cpm_gt1_df$group)
            )
            pct_cpm_gt1_df$sample_name <-
                factor(pct_cpm_gt1_df$sample_name,
                       levels = unique(pct_cpm_gt1_df$sample_name)
                )
        }

        if (!is.null(abundance)) {
            if (is.matrix(abundance)) {
                abundance <- as.data.frame(abundance)
            }
            ## filter out all zero abundance
            # abundance <- abundance[rowSums(abundance)> 0, ]
            pct_tpm_gt1 <- vapply(
                abundance,
                function(.x) {
                    sum(.x >=
                            min_tpm)
                },
                numeric(1)
            ) / nrow(abundance) * 100

            pct_tpm_gt1_df <- data.frame(
                sample_name = names(pct_tpm_gt1),
                percent = pct_tpm_gt1
            )

            pct_tpm_gt1_df <- merge(pct_tpm_gt1_df,
                                    col_data,
                                    by = "sample_name",
                                    all.x = TRUE,
                                    sort = FALSE
            )
            pct_tpm_gt1_df <- pct_tpm_gt1_df[order(
                pct_tpm_gt1_df$group,
                pct_tpm_gt1_df$sample
            ), ]
            pct_tpm_gt1_df$group <- factor(pct_tpm_gt1_df$group,
                                           levels = unique(pct_tpm_gt1_df$group)
            )
            pct_tpm_gt1_df$sample_name <-
                factor(pct_tpm_gt1_df$sample_name,
                       levels = unique(pct_tpm_gt1_df$sample_name)
                )
        }

        if (!is.null(counts) && !is.null(abundance)) {
            pct_tpm_gt1_df$type <- "TPM"
            pct_cpm_gt1_df$type <- "CPM"
            cpm_tpm <- rbind(pct_cpm_gt1_df, pct_tpm_gt1_df)

            cpm_tpm$sample_name <-
                factor(cpm_tpm$sample_name,
                       levels = unique(cpm_tpm$sample_name)
                )

            p <- ggplot(cpm_tpm, aes(
                x = sample_name,
                y = percent,
                color = group
            )) +
                geom_point() +
                xlab("Sample") +
                ylab(paste0("%Gene (\U2265 threshold)")) +
                facet_wrap(~type) +
                theme(
                    text = element_text(size = 8),
                    axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 0,
                        size = 6
                    )
                )
        } else if (!is.null(counts)) {
            p <- ggplot(
                pct_cpm_gt1_df,
                aes(
                    x = sample_name,
                    y = percent,
                    color = group
                )
            ) +
                geom_point() +
                xlab("Sample") +
                ylab(paste0("%Gene (\U2265 ", min_cpm, " CPM)")) +
                theme(
                    text = element_text(size = 8),
                    axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 0,
                        size = 6
                    )
                )
        } else {
            p <- ggplot(
                pct_tpm_gt1_df,
                aes(
                    x = sample_name,
                    y = percent,
                    color = group
                )
            ) +
                geom_point() +
                xlab("Sample") +
                ylab(paste0("%Gene (\U2265 ", min_cpm, " TPM)")) +
                theme(
                    text = element_text(size = 8),
                    axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 0,
                        size = 6
                    )
                )
        }
        return(p)
    }

#' Check sample similarity and variation
#'
#' Perform sample-level exploratory analysis of RNA-seq data, generating heatmap
#' showing sample distances and PCA plot showing sample variations. Internally,
#' DESeq2 is used for vst transformation of count data.
#'
#' @inheritParams plot_assignment_stat
#' @param silent A logical(1), specify whether to draw the plot. It is useful
#'   to set it to FALSE useful when using the gtable output.
#' @return A list of a ggplot object and a [gtable::gtable()] object.
#' \describe{
#'   \item{pca}{A *ggplot* object containing the PCA score plot showing sample
#'              similarity}
#'   \item{heatmap}{A *gtable* object containing the heatmap showing
#'                  pairwise sample distances}
#' }
#' @export
#' @importFrom edgeR cpm
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom SummarizedExperiment assay
#' @import DESeq2 ggplot2
#' @importFrom grDevices rainbow
#' @importFrom stats sd
#'
#' @examplesIf require("patchwork") && require("ggplotify")
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
#' p<- plot_pca_heatmap(SummarizedCounts = sc,
#'                      silent = TRUE)
#' wrap_plots(p[["pca"]], as.ggplot(p[["heatmap"]]), ncol = 1)
plot_pca_heatmap <-
    function(SummarizedCounts = NULL,
             silent = TRUE) {
        stopifnot(is(SummarizedCounts, "SummarizedCounts"))
        counts <- SummarizedCounts$get_salmon_counts()
        col_data <- SummarizedCounts$col_data
        stopifnot(nrow(counts) > 1)

        keep <- filterByExpr(counts, group = col_data$group)
        counts <- counts[keep, ]

        col_data$group <- factor(col_data$group,
                                 levels = unique(col_data$group)
        )

        dds <-
            DESeqDataSetFromMatrix(
                countData = as.matrix(round(counts)),
                colData = col_data,
                design = ~ 0 + group
            )
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds,
                                   fitType = "parametric",
                                   maxit = 1000
        )
        ## exploratory analysis
        vsd <- tryCatch({vst(dds, blind = TRUE)},
                        warning = function(w) {
                            message("Less than 1000 rows in the count table!")
                        },
                        error = function(e) {
                            varianceStabilizingTransformation(dds,
                                                              blind = TRUE)
                        })

        sampleDists <- dist(t(assay(vsd)))

        ## Heatmap showing sample distances
        distancePlot <- function(sampleDists, sampleNames, col_data) {
            sampleDistMatrix <- as.matrix(sampleDists)
            rownames(sampleDistMatrix) <- sampleNames
            colnames(sampleDistMatrix) <- sampleNames
            colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

            # add group annotation
            anno <- col_data[, c("sample_name", "group")]
            rownames(anno) <- anno$sample_name
            anno <- anno[, "group", drop = FALSE]
            colnames(anno) <- "Group"
            n_groups <- nlevels(anno$Group)
            if (n_groups <= 20) {
                distinct_cols <- c(
                    "#e6194b",
                    "#3cb44b",
                    "#ffe119",
                    "#4363d8",
                    "#f58231",
                    "#911eb4",
                    "#46f0f0",
                    "#f032e6",
                    "#bcf60c",
                    "#fabebe",
                    "#008080",
                    "#e6beff",
                    "#9a6324",
                    "#fffac8",
                    "#800000",
                    "#aaffc3",
                    "#808000",
                    "#ffd8b1",
                    "#000075",
                    "#808080",
                    "#000000",
                    "#ffffff"
                )[seq_len(n_groups)]
            } else {
                distinct_cols <- rainbow(n_groups)
            }
            names(distinct_cols) <- levels(anno$Group)
            # define the colours
            anno_col <- list(Group = distinct_cols)

            p <- pheatmap(
                sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                annotation_row = anno,
                annotation_col = anno,
                annotation_colors = anno_col,
                fontsize = 6,
                silent = silent,
                main = "Heatmap showing sample distances",
                color = colors
            )
        }

        p_pheatmap <- distancePlot(
            sampleDists = sampleDists,
            sampleNames = vsd$sample_name,
            col_data = col_data
        )

        ## remove genes with sd = 0 across samples
        vsd_filter <- t(assay(vsd))
        vsd_filter <-
            vsd_filter[, unname(vapply(as.data.frame(vsd_filter), sd,
                                       FUN.VALUE = numeric(1)))!= 0]

        ## PCA plot showing PC1 and PC2 only
        pca <-
            prcomp(vsd_filter,
                   scale = TRUE,
                   center = TRUE,
                   retx = TRUE
            )
        pc12 <- as.data.frame(pca$x[, seq_len(2)])
        #colnames(pc12) <- c("PC1", "PC2")
        pc12 <- cbind(pc12, col_data)
        pc12_var <-
            round(pca$sdev[seq_len(2)]^2 / (sum(pca$sdev^2)) * 100, digits = 2)
        pc12$group <- factor(pc12$group, levels = unique(pc12$group))
        pc12$sample_name <- factor(pc12$sample_name,
                                   levels = unique(pc12$sample_name)
        )
        p_pca <- ggplot(pc12, aes(
            x = PC1,
            y = PC2,
            color = group,
            label = sample_name
        )) +
            geom_text_repel(size = 2.5, show.legend = FALSE) +
            geom_point() +
            xlab(paste0("PC1 (", pc12_var[1], "%)")) +
            guides(color = guide_legend(title = "Group")) +
            ylab(paste0("PC2 (", pc12_var[2], "%)")) +
            ggtitle("PCA score plot") +
            theme(
                plot.title = element_text(hjust = 0.5),
                axis.text = element_text(size = 8),
                axis.title = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 8)
            )

        list(pca = p_pca, heatmap = p_pheatmap)
    }

