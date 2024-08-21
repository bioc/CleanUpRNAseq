#' @title SummarizedCounts Object
#'
#' @description
#'   A class for storing and retrieving summarized RNA-seq data.
#' @section Constructors:
#'   Create an object of [SummarizedCounts] by setting `lib_strand` and
#'   `col_data`:
#'   ```r
#'   x <- SummarizedCounts$new(lib_strand = 0,
#'                             col_data = "colData in a tab-delimited txt file")
#'   ````
#'   Alternatively, users can call the wrapper function
#'   [create_summarizedcounts()] to construct an object of
#'   [SummarizedCounts].
#'
#' @section Setting and getting fields:
#'
#'   The library stranded information, `colData`, and corrected expression
#'   data of a SummarizedCounts object, `x`, can be simply accessed as follows:
#'   '''r
#'   x$lib_strand
#'   x$col_data
#'   x$global_correction
#'   x$gc_correction
#'   x$ir_correction
#'   x$stranded_correction
#'   ```
#'
#'   The field and sub-fields of featureCounts-summarized RNA-seq data can be
#'   set and get as follows. The feature_counts list contains the following
#'   elements: "gene", "exon", "intergenic_region", "intronic_region", "rRNA",
#'   "mitochondrion", "gtf", "chloroplast". For genomic features, gene, exon,
#'   intron, rRNA, mitochondrion, and chloroplast, only the `stat` elements of
#'   the [Rsubread::featureCounts()] oupt is stored; for intergenic region,
#'   the `stat`, `annotation` and `counts` elements are stored; for GTF-based
#'   summarization, the `stat` and `counts` elements are stored.
#'
#'   ```r
#'   ## the whole feature_counts list
#'   x$set_feature_counts(feature_counts_list)
#'   x$get_feature_counts()
#'
#'   ## gene-coding regions: intron + exons
#'   x$get_gene_stat ()
#'
#'   ## exons
#'   x$get_exon_stat()
#'
#'   ## intergenic region summaries
#'   x$get_ir_stat()
#'   x$get_ir_counts()
#'   x$get_ir_anno()
#'
#'   ## intronic regions
#'   x$get_intron_stat()
#'
#'   ## rRNA-coding regions
#'   x$get_rRNA_stat()
#'
#'   ## mitochondrion
#'   x$get_mt_stat()
#'
#'   ## chloroplast
#'   x$get_ct_stat()
#'
#'   ## GTF meta-gene, exons as features
#'   x$get_gtf_stat()
#'   x$get_gtf_counts()
#'   ```
#'
#'   Adding Salmon quantification data (setting the library type information
#'   to the opposite type) to the SummarizedCounts object. For library type
#'   information, see https://salmon.readthedocs.io/en/latest/library_type.html.
#'   ```r
#'   x$set_salmon_quant(tximport_list)
#'   x$get_salmon_quant()
#'   x$get_salmon_counts()
#'   x$get_salmon_abundance()
#'   x$get_salmon_length()
#'   ````
#'
#'   Adding Salmon quantification data (setting the library type information
#'   to the opposite type) to the SummarizedCounts object. For library type
#'   information, see https://salmon.readthedocs.io/en/latest/library_type.html.
#'   ```r
#'   x$set_salmon_quant_opposite(tximport_list)
#'   x$get_salmon_quant_opposite()
#'   x$get_salmon_opposite_counts()
#'   x$get_salmon_opposite_abundance()
#'   ```
#'
#'   Adding corrected expression data to the SummarizedCounts object:
#'   ```r
#'   x$set_global_correction(corrected_expr)
#'   x$set_gc_correction(corrected_expr)
#'   x$set_ir_correction(corrected_expr)
#'   x$set_stranded_correction(salmon_correction)
#'   ```
#' @export
#' @importFrom R6 R6Class
#' @docType class
#' @format An R6 class

SummarizedCounts <- R6::R6Class("SummarizedCounts",
                                public = list(
#' @field lib_strand (`integer(1)`)\cr
#'  Library strandedness, which can be 0, 1, or 2. See
#'  [create_summarizedcounts()].
    lib_strand = 0,

#' @field col_data (`data.frame`)\cr
#'  colData for RNA-seq data, see
#'  [create_summarizedcounts()].
    col_data = NULL,

#' @field feature_counts (`list()`) FeatureCounts summaries for different types of genomic
#'  features.
    feature_counts = NULL,

#' @field salmon_quant (`list()`)\cr
#'  Salmon quant using the actual library strandedness,imported by tximport.
    salmon_quant = NULL,

#' @field salmon_quant_opposite (`list()`)\cr
#'  Salmon quant using strandedness info opposite to actual
#'  library strandedness, imported by tximport.
    salmon_quant_opposite = NULL,

#' @field global_correction (`matrix()`)\cr
#'  Count matrix corrected for gDNA contamination by using
#'  the "Global" method.
    global_correction = NULL,

#' @field gc_correction (`matrix()`)\cr
#'  Count matrix corrected for gDNA contamination by using
#'  the "GC%" method.
    gc_correction = NULL,

#' @field ir_correction (`matrix()`)\cr
#'  Count matrix corrected for gDNA contamination by using
#'  the "IR%" method.
    ir_correction = NULL,

#' @field stranded_correction (`matrix()`)\cr
#'  Count matrix corrected for gDNA contamination by using
#'  the method dedicated to stranded RNA-seq data.
    stranded_correction = NULL,


#' @description
#' Creates a new instance of this [R6][R6::R6Class] class.
#'
#' Note that this object is typically constructed via a wrapper function,
#' [create_summarizedcounts()].
#'
#' @param lib_strand (`integer(1)`) The library's strandedness. It has
#'   three possible values:
#'   - 0: unstranded, the default;
#'   - 1: stranded, read 1 (or single-end read) comes from the forward strand;
#'   - 2: reversely stranded, read 1 (or single-end read) comes from the
#'        reverse strand
#'   For more details, See
#'   https://sailfish.readthedocs.io/en/master/library_type.html.
#'  `colData`
#' @param col_data (`data.frame()`) A data frame with rows corresponding to
#'   samples. For unstranded RNA-seq data, it at least contains the
#'   following columns: `sample_name`, `BAM_file`, `group`,
#'   `salmon_quant_file`, and `batch` if the data were generated in more
#'    than one batches. For stranded RNA-seq
#'   data, an extra column, `salmon_quant_file_opposite_strand`, should be
#'   included.

    initialize = function(lib_strand, col_data) {
        validate(lib_strand, col_data)

        ## convert "group and batch to factor
        if ("batch" %in% colnames(col_data)) {
            col_data$batch <- as.factor(col_data$batch)
        }
        col_data$group <- as.factor(col_data$group)

        self$lib_strand <- lib_strand
        self$col_data <- col_data

        ## featureCounts summaries for different types of genomic features
        self$feature_counts <- list(gene = data.frame(),
                                    exon = data.frame(),
                                    intergenic_region = data.frame(),
                                    intronic_region = list(counts = matrix(),
                                                      annotation = data.frame(),
                                                       stat = data.frame()),
                                    rRNA = data.frame(),
                                    mitochondrion = data.frame(),
                                    chloroplast = data.frame(),
                                    gtf = list(counts = matrix(),
                                               stat = data.frame()))

        ## Salmon quant using the actual library strandedness, imported by
        ## tximport
        self$salmon_quant <- list(counts = matrix(),
                                  abundances = matrix(),
                                  length = matrix())

        ## Salmon quant using strandedness info opposite to actual library
        ## strandedness, imported by tximport
        if (lib_strand != 0) {
            self$salmon_quant_opposite <- list(counts = matrix(),
                                               abundances = matrix(),
                                               length = matrix())
        }

        ## corrected expression
        self$global_correction <- matrix()
        self$gc_correction <- matrix()
        self$ir_correction <- matrix()

        self$stranded_correction <- list(counts = matrix(),
                                         abundances = matrix(),
                                         length = matrix())
    },

#' @description
#'  Add IR% to the col_data by merging data frames
#'
#' @param IR_rate (`data.frame()`)\cr
#'  A data frame with sample_name as rownames, and a single column `IR_rate`,
#'  stroing the percentage of reads mapping to intergenic regions.
#' @return An object of [SummarizedCounts].
#'
    add_ir_rate = function(IR_rate) {
        if (!"IR_rate" %in% colnames(self$col_data)){
            self$col_data <- merge(self$col_data, IR_rate,
                                   by.x = "sample_name",
                                   by.y = "row.names",
                                   sort = FALSE, all.x = TRUE)
            stopifnot(self$col_data$sample_name ==
                          colnames(self$get_salmon_counts()))
            stopifnot(self$col_data$sample_name ==
                          colnames(self$get_gtf_counts()))
        }
        invisible(self)
    },
#' @description
#'   Set feature_counts
#'
#' @param feature_counts_list (`data.frame()`)\cr
#'   The feature_counts list contains the following elements: "gene", "exon",
#'   "intergenic_region", "intronic_region", "rRNA", "mitochondrion", "gtf",
#'   "chloroplast". For genomic features, gene, exon, intron, rRNA,
#'   mitochondrion, and chloroplast, only the `stat` elements of the
#'   [Rsubread::featureCounts()] oupt is stored; for intergenic region,
#'   the `stat`, `annotation` and `counts` elements are stored; for GTF-based
#'   summarization, the `stat` and `counts` elements are stored.
#' @return An object of [SummarizedCounts].
    set_feature_counts = function(feature_counts_list) {
        list_names <- c("gene",
                        "exon",
                        "intergenic_region",
                        "intronic_region",
                        "rRNA",
                        "mitochondrion",
                        "gtf",
                        "chloroplast")
        stopifnot(is.list(feature_counts_list),
                  all(names(feature_counts_list) %in% list_names),
                  is.list(feature_counts_list$intergenic_region),
                  is.list(feature_counts_list$gtf),
                  length(feature_counts_list$intergenic_region) == 3,
                  length(feature_counts_list$gtf) == 2,
                  all(c("counts", "stat") %in%
                          names(feature_counts_list$gtf)),
                  all(c("counts", "annotation", "stat") %in%
                          names(feature_counts_list$intergenic_region))
                  )
        self$feature_counts <- feature_counts_list
        invisible(self)
    },

#' @description
#'   set salmon_quant
#'
#' @param tximport_list (`list()`)\cr
#'   A three-element [tximport()] list aggregated from Salmon quantification
#'   results genrated by using the correct library strandedness information,
#'   containing `counts`, `abundance`, and `length`.
#' @return An object of [SummarizedCounts].
#'
    set_salmon_quant = function(tximport_list) {
        stopifnot(is.list(tximport_list),
                  setequal(names(tximport_list),
                           c("counts", "abundance", "length"))
            )
        self$salmon_quant <-tximport_list
        invisible(self)
    },
#' @description
#'   set salmon_quant_opposite
#'
#' @param tximport_list (`list()`)\cr
#'   A three-element [tximport()] list aggregated from Salmon quantification
#'   results genrated by using the opposite library strandedness information,
#'   containing `counts`, `abundance`, and `length`.
#' @return An object of [SummarizedCounts].
#'
    set_salmon_quant_opposite = function(tximport_list) {
        stopifnot(is.list(tximport_list),
                  setequal(names(tximport_list),
                           c("counts", "abundance", "length"))
        )
        self$salmon_quant_opposite <- tximport_list
        invisible(self)
    },

#' @description
#'   Set express matrix by using the "Global" method
#'
#' @param corrected_expr (`matrix()`)\cr
#'   A count matrix containing expression corrected by the "Global" method.
#' @return An object of [SummarizedCounts].
#'
    set_global_correction = function(corrected_expr){
        stopifnot(is.matrix(corrected_expr),
                  colnames(corrected_expr) == self$col_data$sample_name)
        self$global_correction <- corrected_expr
        invisible(self)
    },

#' @description
#'   Set express matrix by using the "GC%" method
#'
#' @param corrected_expr (`matrix()`)\cr
#'   A count matrix containing expression corrected by the "GC%" method.
#' @return An object of [SummarizedCounts].
#'
    set_gc_correction = function(corrected_expr){
        stopifnot(is.matrix(corrected_expr),
                  colnames(corrected_expr) == self$col_data$sample_name)
        self$gc_correction <- corrected_expr
        invisible(self)
    },

#' @description
#'   Set express matrix by using the "IR%" method
#'
#' @param corrected_expr (`matrix()`)\cr
#'   A matrix containing expression corrected by the "Global" method in the form
#'   of log2(count per million).
#' @return An object of [SummarizedCounts].
#'
    set_ir_correction = function(corrected_expr){
        stopifnot(is.matrix(corrected_expr),
                  colnames(corrected_expr) == self$col_data$sample_name)
        self$ir_correction <- corrected_expr
        invisible(self)
    },

#' @description
#'   Set express matrix by using the method dedicated to stranded RNA-seq data
#'
#' @param corrected_expr (`matrix()`)\cr
#'   A count matrix containing expression corrected by the method dedicated to
#'   stranded RNA-seq data.
#' @return An object of [SummarizedCounts].
#'
    set_stranded_correction = function(corrected_expr) {
        stopifnot(is.list(corrected_expr),
                  length(corrected_expr) == 3,
                  setequal(c("counts", "abundance", "length"),
                           names(corrected_expr)))
        self$stranded_correction <- corrected_expr
    },

#' @description
#'   Get featureCounts output as a whole list
#'
#' @return A list containing the following elements: "gene", "exon",
#'   "intergenic_region", "intronic_region", "rRNA", "mitochondrion", "gtf",
#'   "chloroplast". For genomic features, gene, exon, intron, rRNA,
#'   mitochondrion, and chloroplast, only the `stat` elements of the
#'   [Rsubread::featureCounts()] output is stored; for intergenic region,
#'   the `stat`, `annotation` and `counts` elements are stored; for GTF-based
#'   summarization, the `stat` and `counts` elements are stored.
#'
    get_feature_counts = function() {
        self$feature_counts
    },

#' @description
#'   Get featureCounts `stat` for the `gene-coding regions`
#'
#' @return only the `stat` elements of the [Rsubread::featureCounts()] output
#'  for the `gene-coding regions`.
    get_gene_stat = function(){
        self$feature_counts$gene$stat
    },

#' @description
#'   Get featureCounts `stat` for the `exon regions`
#'
#' @return only the `stat` elements of the [Rsubread::featureCounts()] output
#'  for the `exon regions`
    get_exon_stat = function(){
        self$feature_counts$exon
    },

#' @description
#'   Get featureCounts `stat` for the `intergenic regions`
#'
#' @return only the `stat` elements of the [Rsubread::featureCounts()] output
#'  for the `intergenic regions`.
    get_ir_stat = function(){
        self$feature_counts$intergenic_region$stat
    },

#' @description
#'   Get featureCounts `counts` for the `intergenic regions`
#'
#' @return only the `counts` elements of the [Rsubread::featureCounts()] output
#'  for the `intergenic regions`.
    get_ir_counts = function(){
        self$feature_counts$intergenic_region$counts
    },

#' @description
#'   Get featureCounts `annotation` for the `intergenic regions`
#'
#' @return only the `annotation` elements of the [Rsubread::featureCounts()] output
#'  for the `intergenic regions`.
    get_ir_anno = function(){
        self$feature_counts$intergenic_region$annotation
    },

#' @description
#'   Get featureCounts `stat` for the `intronic regions`
#'
#' @return only the `stat` elements of the [Rsubread::featureCounts()] output
#'  for the `intronic regions`.
    get_intron_stat = function(){
        self$feature_counts$intronic_region
    },

#' @description
#'   Get featureCounts `stat` for the `rRNA-encoding regions`
#'
#' @return only the `stat` elements of the [Rsubread::featureCounts()] output
#'  for the `rRNA-encoding regions`.
    get_rRNA_stat = function(){
        self$feature_counts$rRNA
    },

#' @description
#'   Get featureCounts `stat` for the `mitochondiral genome`
#'
#' @return only the `stat` elements of the [Rsubread::featureCounts()] output
#'  for the `mitochondiral genome`.
    get_mt_stat = function(){
        self$feature_counts$mitochondrion
    },
#' @description
#'   Get featureCounts `stat` for the `chloroplast genome`
#'
#' @return only the `stat` elements of the [Rsubread::featureCounts()] output
#'  for the `chloroplast genome`.
    get_ct_stat = function(){
        self$feature_counts$chloroplast
    },

#' @description
#'   Get featureCounts `stat` for the `GTF metagene`
#'
#' @return only the `stat` elements of the [Rsubread::featureCounts()] output
#'  for the `GTF metagenes`.
    get_gtf_stat = function(){
        self$feature_counts$gtf$stat
    },

#' @description
#'   Get featureCounts `counts` for the `GTF metagene`
#'
#' @return only the `counts` elements of the [Rsubread::featureCounts()] output
#'  for the `GTF metagenes`.
    get_gtf_counts = function(){
        self$feature_counts$gtf$counts
    },

#' @description
#'  Get [tximport()] output as a whole list, including the three matrices
#'  `counts`, `abundance`, and `length`
#'
#' @return A list of three elements, `counts`, `abundance`, and `length` of
#'  the [tximport::tximport()] output.
    get_salmon_quant = function(){
        self$salmon_quant
    },

#' @description
#'  Get [tximport()] output as a whole list, only including the `counts` matrix
#'
#' @return only the  `counts` element of the [tximport::tximport()] output.
    get_salmon_counts = function(){
        self$salmon_quant$counts
    },

#' @description
#'  Get [tximport()] output as a whole list, only including the `abundance`
#'  matrix
#'
#' @return only the  `abundance` element of the [tximport::tximport()] output.
    get_salmon_abundance = function(){
        self$salmon_quant$abundance
    },

#' @description
#'  Get [tximport()] output as a whole list, only including the `length`
#'  matrix
#'
#' @return only the  `length` element of the [tximport::tximport()] output.
    get_salmon_length = function(){
        self$salmon_quant$length
    },

#' @description
#'  Get [tximport()] output as a whole list, including the three matrices
#'  `counts`, `abundance`, and `length`
#'
#' @return A list of three elements, `counts`, `abundance`, and `length` of
#'  the [tximport::tximport()] output.
    get_salmon_quant_opposite = function(){
        self$salmon_quant_opposite
    },

#' @description
#'  Get [tximport()] output as a whole list, only including the `count`
#'  matrix
#'
#' @return only the  `counts` element of the [tximport::tximport()] output.
    get_salmon_opposite_counts = function(){
        self$salmon_quant_opposite$counts
    },

#' @description
#'  Get [tximport()] output as a whole list, only including the `abundance`
#'  matrix
#'
#' @return only the  `abundance` element of the [tximport::tximport()] output.
    get_salmon_opposite_abundance = function(){
        self$salmon_quant_opposite$abundance
    })
)

validate <- function(lib_strand, col_data) {
    stopifnot(is.numeric(lib_strand),
              lib_strand %in% c(0, 1, 2),
              length(lib_strand) == 1,
              nrow(col_data) >= 1)

    col_names <- c("sample_name", "BAM_file", "group",
                   "salmon_quant_file",
                   "salmon_quant_file_opposite_strand",
                   "batch")

    if (lib_strand == 0) {
        stopifnot(is.data.frame(col_data),
                  all(col_names[seq_len(4)] %in% colnames(col_data)))
    } else {
        stopifnot(is.data.frame(col_data),
                  all(col_names[seq_len(5)] %in% colnames(col_data)))
        stopifnot(file.exists(col_data$salmon_quant_file_opposite_strand))
        if (any(duplicated(col_data$salmon_quant_file_opposite_strand))) {
            stop("Some salmon_quant_file_opposite_strand are not unique")
        }
    }

    stopifnot(file.exists(col_data$BAM_file),
              file.exists(col_data$salmon_quant_file))

    if (any(is.na(col_data)) ||
        any(col_data == "") || any(col_data == " ")) {
        stop("Missing value(s) found in the col_data")
    }

    if (any(duplicated(col_data$BAM_file)) ||
        any(duplicated(col_data$salmon_quant_file)) ||
        any(duplicated(col_data$sample_name))) {
        stop("Some BAM files, salmon_quant_file, ",
        "or sample names are not unique")
    }
    stopifnot(nlevels(as.factor(col_data$group)) > 1)

    if ("batch" %in% colnames(col_data) &&
        nlevels(as.factor(col_data$batch)) < 2) {
        stop("Only one batch exists, so delete it from the colData")
    }
}



#' Create an object of SummarizedCounts
#'
#' Set up an analysis by creating an object of [SummarizedCounts], which
#' is used to store summary output from featureCounts and tximport, and store
#' gDNA-corrected expression data. This is a convenience wrapper for
#' `SummarizedCounts$new()`
#'
#' @param lib_strand An integer(1), specifying the library's strandedness. It
#'   has three possible values:
#'   - 0: unstranded, the default;
#'   - 1: stranded, read 1 (or single-end read) comes from the forward strand;
#'   - 2: reversely stranded, read 1 (or single-end read) comes from the
#'        reverse strand
#'   For more details, See
#'   https://sailfish.readthedocs.io/en/master/library_type.html.
#' @param colData A data frame with rows corresponding to samples. For
#'   unstranded RNA-seq data, it at least contains the following columns:
#'   `sample_name`, `BAM_file`, `group`, `salmon_quant_file`, and `batch` if
#'   the data were generated in more than one batches. For stranded RNA-seq
#'   data, an extra column, `salmon_quant_file_opposite_strand`, should be
#'   included.
#'
#' @return An object of [SummarizedCounts], which is used to store
#'   summary output from featureCounts and tximport, and store gDNA-corrected
#'   expression data.
#' @export
#'
#' @examples
#' in_dir <- system.file("extdata", package = "CleanUpRNAseq")
#' BAM_file <- dir(in_dir, ".bam$", full.name = TRUE)
#' salmon_quant_file <- dir(in_dir, ".sf$", full.name = TRUE)
#' sample_name = gsub(".+/(.+?).srt.bam", "\\1", BAM_file)
#' salmon_quant_file_opposite_strand <- salmon_quant_file
#' col_data <- data.frame(sample_name = sample_name,
#'                        BAM_file = BAM_file,
#'                        salmon_quant_file = salmon_quant_file,
#'                        salmon_quant_file_opposite_strand =
#'                            salmon_quant_file_opposite_strand,
#'                        group = c("CD1N", "CD1P"))
#'
#' sc <- create_summarizedcounts(lib_strand = 0, colData = col_data)
#'
create_summarizedcounts <- function(lib_strand = 0,
                                    colData = NULL) {
    sc <- SummarizedCounts$new(lib_strand,
                               colData)
    invisible(sc)
}
