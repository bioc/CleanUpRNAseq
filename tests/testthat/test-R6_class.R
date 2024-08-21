in_dir <- system.file("extdata", package = "CleanUpRNAseq")
BAM_file <- dir(in_dir, ".bam$", full.names = TRUE)
salmon_quant_file <- dir(in_dir, ".sf$", full.names = TRUE)
sample_name = gsub(".+/(.+?).srt.bam", "\\1", BAM_file)
salmon_quant_file_opposite_strand <- salmon_quant_file
group <-  c("CD1N", "CD1P")
col_data <- data.frame(sample_name = sample_name,
                       BAM_file = BAM_file,
                       salmon_quant_file = salmon_quant_file,
                       salmon_quant_file_opposite_strand =
                           salmon_quant_file_opposite_strand,
                       group = group)

test_that("create_summarizedcounts works", {
  expect_invisible(create_summarizedcounts(lib_strand = 0, col_data))
  expect_invisible(create_summarizedcounts(lib_strand = 1, col_data))
  expect_invisible(create_summarizedcounts(lib_strand = 2, col_data))
  expect_error(create_summarizedcounts(lib_strand = 3, col_data))
})

# BAM files not existing
BAM_file_fake <- gsub(".bam$", ".1.bam",
                      dir(in_dir, ".bam$", full.names = TRUE))
col_data2 <- data.frame(sample_name = sample_name,
                        BAM_file = BAM_file_fake,
                        salmon_quant_file = salmon_quant_file,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand,
                        group = group)

salmon_quant_file_fake <- gsub(".sf$", ".1.sf",
                          dir(in_dir, ".sf$", full.names = TRUE))

# quant.sf file not existing
col_data3 <- data.frame(sample_name = sample_name,
                        BAM_file = BAM_file,
                        salmon_quant_file = salmon_quant_file_fake,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand,
                        group = group)

## missing column
col_data4 <- data.frame(BAM_file = BAM_file,
                        salmon_quant_file = salmon_quant_file,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand,
                        group = group)

## duplicated sample name
sample_name_dup <- sample_name[c(1,1)]
col_data5 <- data.frame(sample_name = sample_name_dup,
                        BAM_file = BAM_file,
                        salmon_quant_file = salmon_quant_file,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand,
                        group = group)
BAM_file_dup <- BAM_file[c(1,1)]
col_data6 <- data.frame(sample_name = sample_name,
                        BAM_file = BAM_file_dup,
                        salmon_quant_file = salmon_quant_file,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand,
                        group = group)

salmon_quant_file_dup <- salmon_quant_file[c(1,1)]
col_data7 <- data.frame(sample_name = sample_name,
                        BAM_file = BAM_file,
                        salmon_quant_file = salmon_quant_file_dup,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand,
                        group = group)

salmon_quant_file_opposite_strand_dup <-
    salmon_quant_file_opposite_strand[c(1,1)]
col_data8 <- data.frame(sample_name = sample_name,
                        BAM_file = BAM_file,
                        salmon_quant_file = salmon_quant_file,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand_dup,
                        group = group)

## single group or batch
group_single <-group[c(1,1)]
col_data9 <- data.frame(sample_name = sample_name,
                        BAM_file = BAM_file,
                        salmon_quant_file = salmon_quant_file,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand,
                        group = group_single)

## missing value in colData
sample_name_missing1 <- c(sample_name[1], "")
col_data10 <- data.frame(sample_name = sample_name_missing1,
                        BAM_file = BAM_file,
                        salmon_quant_file = salmon_quant_file,
                        salmon_quant_file_opposite_strand =
                            salmon_quant_file_opposite_strand,
                        group = group)
sample_name_missing2 <- c(sample_name[1], " ")
col_data11 <- data.frame(sample_name = sample_name_missing2,
                         BAM_file = BAM_file,
                         salmon_quant_file = salmon_quant_file,
                         salmon_quant_file_opposite_strand =
                             salmon_quant_file_opposite_strand,
                         group = group)
sample_name_missing3 <- c(sample_name[1], NA)
col_data12 <- data.frame(sample_name = sample_name_missing3,
                         BAM_file = BAM_file,
                         salmon_quant_file = salmon_quant_file,
                         salmon_quant_file_opposite_strand =
                             salmon_quant_file_opposite_strand,
                         group = group)

test_that("create_summarizedcounts error", {
    expect_error(create_summarizedcounts(lib_strand = c(0,1), col_data))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data2))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data3))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data4))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data5))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data6))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data7))
    expect_no_error(create_summarizedcounts(lib_strand = 0, col_data8))
    expect_error(create_summarizedcounts(lib_strand = 1, col_data8))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data9))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data10))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data11))
    expect_error(create_summarizedcounts(lib_strand = 0, col_data12))
})
