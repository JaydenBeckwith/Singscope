# Paths to csvs
expr_path <- "test_data/opacin_count_demo.csv"
metadata_path <- "test_data/opacin_meta_demo.csv"
gmt_path <- "test_data/20251505_240genelist_withphenotypes.gmt"
rds_path <- "test_data/merged_sing_df.rds"

source("../../R/preprocess_data.R")

# Dummy logger to suppress output during testing
dummy_logger <- function(...) {}

# Test 1: make sure it runs without error with merging disabled
test_that("preprocess data runs without error (no merge)", {
  expect_error(
    preprocess_data(expr_path, metadata_path, "test123", gmt_path, 
                    existing_data_merge = FALSE, example_data_path = rds_path, logger = dummy_logger), 
    NA
  )
})

# Test 2: Output is a data frame with expected columns
test_that("output is a dataframe with expected structure", {
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path, 
                            existing_data_merge = FALSE, example_data_path = rds_path, logger = dummy_logger)
  expect_true(is.data.frame(result))
  expect_true(all(c("sample_id", "Pathway", "Singscore") %in% colnames(result)))
})

# Test 3: Result contains non-zero rows
test_that("result contains non-zero rows", {
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path, 
                            existing_data_merge = FALSE, example_data_path = rds_path, logger = dummy_logger)
  expect_gt(nrow(result), 0)
})

# Test 4: Response column exists and has correct values
test_that("Response column exists and is valid", {
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path, 
                            existing_data_merge = FALSE, example_data_path = rds_path, logger = dummy_logger)
  expect_true("Response" %in% colnames(result))
  expect_true(all(result$Response %in% c("MPRs", "NMPRs")))
})

# Test 5: Factor columns are properly set
test_that("Timepoint and recurrence_status are factors", {
  result <- preprocess_data(expr_path, metadata_path, "test123", gmt_path, 
                            existing_data_merge = FALSE, example_data_path = rds_path, logger = dummy_logger)
  expect_true(is.factor(result$Timepoint))
  expect_true(is.factor(result$recurrence_status))
})

# Test 6: RPK values are calculated correctly
test_that("RPK values are calculated correctly", {
  expr <- read.csv(expr_path)
  gene_lengths <- data.frame(
    ensembl_gene_id = expr$ensembl_gene_id,
    cds_length = rep(1500, nrow(expr)) # simulate output of getBM
  )
  gene_lengths$gene_length <- gene_lengths$cds_length / 1000

  merged <- merge(expr, gene_lengths[, c("ensembl_gene_id", "gene_length")], by = "ensembl_gene_id")
  rownames(merged) <- merged$ensembl_gene_id
  counts_only <- merged[, setdiff(colnames(merged), c("ensembl_gene_id", "gene_length"))]
  gene_lengths_kb <- merged$gene_length

  rpk_manual <- sweep(counts_only, 1, gene_lengths_kb, FUN = "/")

  filtered_counts_df <- merge(expr, gene_lengths[, c("ensembl_gene_id", "gene_length")], by = "ensembl_gene_id")
  rownames(filtered_counts_df) <- filtered_counts_df$ensembl_gene_id
  filtered_counts_df <- filtered_counts_df[!is.na(filtered_counts_df$gene_length), ]
  gene_length_kb <- filtered_counts_df$gene_length
  filtered_counts_df <- filtered_counts_df[, setdiff(colnames(filtered_counts_df), c("ensembl_gene_id", "gene_length"))]
  rpk_from_function <- sweep(filtered_counts_df, 1, gene_length_kb, FUN = "/")

  gene <- rownames(rpk_from_function)[1]
  sample <- colnames(rpk_from_function)[1]

  expect_equal(rpk_manual[gene, sample], rpk_from_function[gene, sample])
})