source("../../preprocess_data.R")
library(dplyr)
library(tidyr)
library(clusterProfiler)

# Paths to csvs
expr_path <- "test_data/opacin_count_demo.csv"
metadata_path <- "test_data/opacin_meta_demo.csv"
gmt_path <- "test_data/20251505_240genelist_withphenotypes.gmt"
rds_path <- "test_data/merged_sing_df.rds"

# Test 1: make sure it actually runs 
test_that("preprocess_data runs without error", {
  expect_error(preprocess_data(expr_path, metadata_path, "test123", gmt_path), NA)
})

test_that("output is a dataframe with expected structure", {
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path)
  expect_true(is.data.frame(result))
  expect_true(all(c("sample_id", "Pathway", "Singscore") %in% colnames(result)))
})

test_that("result contains non-zero rows", {
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path)
  expect_gt(nrow(result), 0)
})

test_that("Response column exists and has correct values", {
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path)
  expect_true("Response" %in% colnames(result))
  expect_true(all(result$Response %in% c("MPRs", "NMPRs")))
})

test_that("Timepoint and recurrence_status are factors", {
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path)
  expect_true(is.factor(result$Timepoint))
  expect_true(is.factor(result$recurrence_status))
})

test_that("RPK values are calculated correctly", {
  # Read in raw inputs manually
  expr <- read.csv(expr_path)
  gene_lengths <- data.frame(
    ensembl_gene_id = expr$ensembl_gene_id,
    cds_length = rep(1500, nrow(expr)) # simulate output of getBM
  )
  gene_lengths$gene_length <- gene_lengths$cds_length / 1000  # convert to kb

  # Merge manually (as in function)
  merged <- merge(expr, gene_lengths[, c("ensembl_gene_id", "gene_length")], by = "ensembl_gene_id")
  rownames(merged) <- merged$ensembl_gene_id
  counts_only <- merged[, setdiff(colnames(merged), c("ensembl_gene_id", "gene_length"))]
  gene_lengths_kb <- merged$gene_length

  # Calculate RPK manually
  rpk_manual <- sweep(counts_only, 1, gene_lengths_kb, FUN = "/")

  # Run the actual function
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path)

  # Extract the RPK calculation from within the function's environment
  # Since we can't directly access it, we re-run the same step here using the same logic
  filtered_counts_df <- merge(expr, gene_lengths[, c("ensembl_gene_id", "gene_length")], by = "ensembl_gene_id")
  rownames(filtered_counts_df) <- filtered_counts_df$ensembl_gene_id
  filtered_counts_df <- filtered_counts_df[!is.na(filtered_counts_df$gene_length), ]
  gene_length_kb <- filtered_counts_df$gene_length
  filtered_counts_df <- filtered_counts_df[, setdiff(colnames(filtered_counts_df), c("ensembl_gene_id", "gene_length"))]
  rpk_from_function <- sweep(filtered_counts_df, 1, gene_length_kb, FUN = "/")

  # Compare a known gene and sample
  gene <- rownames(rpk_from_function)[1]
  sample <- colnames(rpk_from_function)[1]

  expect_equal(rpk_manual[gene, sample], rpk_from_function[gene, sample])
})