tmp_dir <- tempdir()
# ## download feaureCounts results
count_url <- paste0(
    "https://zenodo.org/records/11458839/files/",
    "read_count_summary.RData?download=1"
)

featureCounts_f <- file.path(tmp_dir, "read_count_summary.RData")
download.file(count_url,
              destfile = featureCounts_f,
              mode = "wb")
## load read_count_summary.RData to get counts_summary
load(featureCounts_f)
counts_summary <- mapply(function(res, .name){
    res$targets <- NULL
    if (!.name %in% c("intergenic_region", "gtf")) {
        res$counts <- NULL
        res$annotation <- NULL
        return(res$stat)
    }
    if (.name == "gtf") {
        res$annotation <- NULL
        res$counts_junction <- NULL
    }
    return(res)
}, counts_summary, names(counts_summary), SIMPLIFY = FALSE)
feature_counts_list$metadata <- NULL


## download salmon quant
salmon_url <-
    paste0("https://zenodo.org/records/11458839/files/",
           "salmon_quant_summary.RData?download=1")
salmonquant_f <- file.path(tmp_dir, "salmon_quant_summary.RData")
download.file(salmon_url, destfile = salmonquant_f,
              mode = "wb")
load(salmonquant_f)
salmon_quant$countsFromAbundance <- NULL

## subset data for intergenic regions
ir_counts <- counts_summary$intergenic_region$counts
ir_counts_nonzero <- ir_counts[rowSums(ir_counts) != 0, ]
ir_counts_2000 <- ir_counts_nonzero[1:2000, ]
ir_anno <- counts_summary$intergenic_region$annotation
ir_anno_nonzero <- ir_anno[rowSums(ir_counts) != 0, ]
ir_anno_2000 <- ir_anno_nonzero[1:2000, ]

all(rownames(ir_counts_2000) == ir_anno_2000$GeneID)
counts_summary$intergenic_region$counts <- ir_counts_2000
counts_summary$intergenic_region$annotation <- ir_anno_2000

## ## subset data for GTF metagenes and Salmon quant
gtf_counts <- counts_summary$gtf$counts
gtf_counts_nonzero <- gtf_counts[rowSums(gtf_counts) != 0, ]

salmon_counts <- salmon_quant$counts
salmon_counts_nonzero <- salmon_counts[rowSums(salmon_counts) != 0, ]

common_genes <- intersect(rownames(gtf_counts_nonzero),
                          rownames(salmon_counts_nonzero))

common_genes_2000 <- common_genes[1:2000]

gtf_counts_nonzero_2000 <-
    gtf_counts_nonzero[rownames(gtf_counts_nonzero) %in%
                           common_genes_2000, ]

salmon_counts_nonzero_2000 <-
    salmon_counts_nonzero[rownames(salmon_counts_nonzero) %in%
                              common_genes_2000, ]
counts_summary$gtf$counts <- gtf_counts_nonzero_2000
salmon_quant$counts <-
    salmon_counts_nonzero_2000[rownames(salmon_counts_nonzero_2000) %in%
                                   common_genes_2000, ]
setequal(rownames(counts_summary$gtf$counts), rownames(salmon_quant$counts))

salmon_quant$abundance <-
    salmon_quant$abundance[rownames(salmon_quant$abundance) %in%
                               common_genes_2000, ]

setequal(rownames(counts_summary$gtf$counts), rownames(salmon_quant$abundance))
salmon_quant$length <-
    salmon_quant$length[rownames(salmon_quant$length) %in%
                            common_genes_2000, ]

setequal(rownames(counts_summary$gtf$counts), rownames(salmon_quant$length))

feature_counts_list <- counts_summary
save(feature_counts_list,
     file = "feature_counts_list.RData",
     compression_level = 9)
save(salmon_quant,
     file = "salmon_quant.RData",
     compression_level = 9)


## GC contents and lengths of intergenic regions
intergenic_gc_url <-
    paste0("https://zenodo.org/records/11458839/files/",
           "GRCh38.intergenic.GC.content.RDS?download=1")
intergenic_gc_f <- file.path(tmp_dir, "GRCh38.intergenic.GC.content.RDS")
download.file(intergenic_gc_url, destfile = intergenic_gc_f,
              mode = "wb")
intergenic_GC <- readRDS(intergenic_gc_f)
intergenic_GC <-
    intergenic_GC[rownames(intergenic_GC) %in%
                      rownames(feature_counts_list$intergenic_region$counts), ]
save(intergenic_GC, file = "intergenic_GC.rda",
     compression_level = 9)


## GC contents of the 2000 genes
gene_gc_url <- paste0("https://zenodo.org/records/11458839/files/",
                      "GRCh38.gene.exon.collapsed.GC.content.RDS?download=1")
gene_gc_f <- file.path(tmp_dir, "GRCh38.gene.exon.collapsed.GC.content.RDS")
download.file(gene_gc_url , destfile = gene_gc_f ,
              mode = "wb")
gene_GC <- readRDS(gene_gc_f )
genenames_2000 <- rownames(salmon_quant$counts)
gene_GC <- gene_GC[rownames(gene_GC) %in% genenames_2000, ]
save(gene_GC, file = "gene_GC.rda",
     compression_level = 9)
