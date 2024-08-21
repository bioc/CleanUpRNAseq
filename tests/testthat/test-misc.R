library("BSgenome.Hsapiens.UCSC.hg38")
ucsc_BSgenome <- BSgenome.Hsapiens.UCSC.hg38
ucsc_seqnames <-
    seqnames(ucsc_BSgenome)[!grepl("_alt|fix|hap\\d+",
                            seqnames(ucsc_BSgenome))]

## Wierd: BSgenome.Hsapiens.UCSC.hg19 has both chrM and chrMT for
## mitochondrial genome. renomve chrMT.

if (all(c("chrM", "chrMT") %in% ucsc_seqnames)) {
    ucsc_seqnames <- ucsc_seqnames[!ucsc_seqnames %in% "chrMT"]
}

ensembl_seqnames <- gsub(
    "^chr", "",
    gsub(
        "chrM$", "MT",
        gsub(
            "v", ".",
            gsub("_random|chr[^_]+_", "", ucsc_seqnames)
        )
    )
)
## for BSgenome.Hsapiens.UCSC.hg19, scaffold seqnames start with lower case,
## should be changed to upper case. For example, change "gl000231" to
## "GL000231". Add a suffix ".1" to these scaffold seqnames.

seqname_alias <- data.frame(ucsc = ucsc_seqnames,
                            ensembl = ensembl_seqnames)
ensembl_BSgenome <- style_BSgenome(
    UCSC_BSgenome = ucsc_BSgenome,
    genome_version = "GRCh38",
    seqname_alias = seqname_alias
)

test_that("style_BSgenome works", {
    expect_s4_class(ensembl_BSgenome, "BSgenome")
    expect_equal(metadata(ensembl_BSgenome)$genome, "GRCh38")
})
