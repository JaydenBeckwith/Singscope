preprocess_data <- function(exprMatrixPath, metadataPath, cohortName, gmtPath, 
                            existing_data_merge = TRUE, example_data_path = NULL, logger = message) {
  
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://may2025.archive.ensembl.org")
  metadata <- read.csv(metadataPath)
  required_cols <- c("sample_id", "study", "recurrence_status", "Gender", "MPRvNMPR", "Timepoint")
  
  missing_cols <- setdiff(required_cols, colnames(metadata))
  
  if (length(missing_cols) > 0) {
    message("The following columns are missing and will be filled with NA: ", paste(missing_cols, collapse = ", "))
    for (col in missing_cols) {
      metadata[[col]] <- NA
    }
  }

  expr_data <- read.csv(exprMatrixPath)
  logger("Starting Ensembl")

  # Retrieve Ensembl IDs and Gene Lengths
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_lengths <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "cds_length"),
    filters = "ensembl_gene_id",
    values = expr_data$ensembl_gene_id,
    mart = mart
  )

  protein_coding <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "biotype",
    values = "protein_coding",
    mart = mart
  )

  logger("got gene lengths")

  gene_lengths <- aggregate(cds_length ~ ensembl_gene_id, data = gene_lengths, max)
  gene_lengths$cds_length <- gene_lengths$cds_length / 1000
  colnames(gene_lengths) <- c("ensembl_gene_id", "gene_length")

  logger("got tpms")
  logger("merging counts")

  filtered_counts_df <- merge(expr_data, gene_lengths, by = "ensembl_gene_id", all.x = TRUE)
  filtered_counts_df <- filtered_counts_df[!is.na(filtered_counts_df$gene_length), ]
  rownames(filtered_counts_df) <- filtered_counts_df$ensembl_gene_id
  filtered_counts_df <- filtered_counts_df %>% dplyr::select(-ensembl_gene_id)
  
  logger("rpk to tpm")
  expr_only <- filtered_counts_df[, setdiff(names(filtered_counts_df), "gene_length")]
  expr_only[] <- lapply(expr_only, as.numeric)
  filtered_counts_df$gene_length <- as.numeric(filtered_counts_df$gene_length)

  rpk <- expr_only / filtered_counts_df$gene_length

  
  logger("passed rpk doing tpms")
  tpm <- sweep(rpk, 2, colSums(rpk), FUN = "/") * 1e6
  
  logger("ensembl to gene")
  ensembl_to_gene <- setNames(protein_coding$hgnc_symbol, protein_coding$ensembl_gene_id)
  tpm$gene_symbol <- ensembl_to_gene[rownames(tpm)]
  tpm <- tpm[!is.na(tpm$gene_symbol) & tpm$gene_symbol != "", ]
  tpm <- tpm[!duplicated(tpm$gene_symbol), ]

  logger("singscores generating..")

  # Singscore analysis
  if (!file.exists(gmtPath)) stop("Invalid GMT path: ", gmtPath)
  PIPdx <- GSEABase::getGmt(gmtPath)

  rownames(tpm) <- tpm$gene_symbol
  tpm$gene_symbol <- NULL

  eranks <- singscore::rankGenes(tpm)
  singscores <- singscore::multiScore(eranks, PIPdx)$Scores %>% as.data.frame()

  transposed_singscores <- t(singscores) %>% as.data.frame()
  transposed_singscores$sample_id <- rownames(transposed_singscores)

  long_singscores <- transposed_singscores %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = "Pathway", values_to = "Singscore")

  new_data <- merge(metadata, long_singscores, by = "sample_id")
  new_data$recurrence_status <- factor(new_data$recurrence_status)
  new_data$Timepoint <- factor(new_data$Timepoint, levels = c("Baseline", "Week 6"))
  new_data$Response <- ifelse(new_data$MPRvNMPR == 1, "MPRs", "NMPRs")

  logger("done")

  if (existing_data_merge) {
  logger(paste("Merging with existing dataset from:", example_data_path))

  if (is.null(example_data_path) || !file.exists(example_data_path)) {
    stop("Example RDS file not found or not specified: ", example_data_path)
  }

  existing_data <- readRDS(example_data_path)

  combined_data <- existing_data %>%
    dplyr::anti_join(new_data, by = "sample_id") %>%
    dplyr::bind_rows(new_data)

  saveRDS(combined_data, example_data_path)

  return(combined_data)
  
} else {
  logger("Skipping merge with example data — using uploaded data only.")
  return(new_data)
}}