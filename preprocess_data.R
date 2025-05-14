# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(biomaRt)
library(dplyr)
library(clusterProfiler)
library(SummarizedExperiment)
library(readr)

# === Load the data ===
singscores <- read.csv("data/neopele_neotrioalone_singscores.csv")
metadata <- read.csv("data/neopele_neotrioalone_meta.csv")
metadata <- metadata %>% select(-X)
metadata

# === Data Transformation ===
# Set the row names to the pathway names
rownames(singscores) <- singscores$X
singscores <- singscores %>% select(-X)

colnames(singscores) <- gsub("^X", "", colnames(singscores))

# === Transpose and prepare for merging ===
transposed_singscores <- t(singscores) %>% as.data.frame()
transposed_singscores$sample_id <- rownames(transposed_singscores)

# === Convert to long format for easier plotting ===
long_singscores <- transposed_singscores %>%
  pivot_longer(cols = -sample_id, names_to = "Pathway", values_to = "Singscore")

# === Merge with Metadata ===
merged_sing_df <- merge(metadata, long_singscores, by = "sample_id")
merged_sing_df$recurrence_status <- factor(merged_sing_df$recurrence_status)

# === Clean up Timepoint names ===
#merged_sing_df$Timepoint <- ifelse(merged_sing_df$Timepoint == "0", "Baseline", "Week 6")
merged_sing_df$Timepoint <- factor(merged_sing_df$Timepoint, levels = c("Baseline", "Week 6"))
unique(merged_sing_df$Timepoint)


merged_sing_df <- merged_sing_df %>%
  mutate(Response = ifelse(MPRvNMPR == 1, "MPRs", "NMPRs"))

merged_sing_df

# Save for the Shiny app
saveRDS(merged_sing_df, "merged_sing_df.rds")


# === Preprocessing Function ===
preprocess_data <- function(exprMatrixPath, metadataPath, cohortName, gmtPath = "data/20241021_188genelist_withphenotypes.gmt") {
  
  # === Load BioMart ===
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://may2025.archive.ensembl.org")
  
  # === 1. Load Metadata ===
  metadata <- read.csv(metadataPath)
  
  # === 2. Validate Metadata Columns ===
  required_cols <- c("sample_id", "study", "recurrence_status", "Gender", "MPRvNMPR", "Timepoint")
  
  # Find missing columns
  missing_cols <- setdiff(required_cols, colnames(metadata))
  
  # If there are missing columns, add them filled with NA
  if (length(missing_cols) > 0) {
    message("The following columns are missing and will be filled with NA: ", paste(missing_cols, collapse = ", "))
    for (col in missing_cols) {
      metadata[[col]] <- NA
    }
  }
  # === 3. Load Gene Expression Data ===
  expr_data <- read.csv(exprMatrixPath)

  message("starting ensembl")
  
  # === 4. Retrieve Ensembl IDs and Gene Lengths ===
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
 
  gene_lengths <- getBM(
    attributes = c("ensembl_gene_id", "cds_length"),
    filters = "ensembl_gene_id",
    values = expr_data$ensembl_gene_id,  
    mart = mart
  )
  
  message("got gene lengths")
  
  # Calculate TPM
  gene_lengths <- aggregate(cds_length ~ ensembl_gene_id, data = gene_lengths, max)
  gene_lengths$cds_length <- gene_lengths$cds_length / 1000
  colnames(gene_lengths) <- c("ensembl_gene_id", "gene_length")
  
  message("got tpms")
  
  message("merging counts")

  # Merge
  filtered_counts_df <- merge(expr_data, gene_lengths, by = "ensembl_gene_id", all.x = TRUE)
  filtered_counts_df <- filtered_counts_df[!is.na(filtered_counts_df$gene_length), ]
  
  message("filtered counts again")
  rownames(filtered_counts_df) <- filtered_counts_df$ensembl_gene_id
  filtered_counts_df <- filtered_counts_df %>% select(-ensembl_gene_id)
  print(head(filtered_counts_df))
  message("rpk to tpm")
  # Compute RPK and TPM
  rpk <- filtered_counts_df[, -c(1, ncol(filtered_counts_df))] / filtered_counts_df$gene_length
  tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6
  
  # Convert Ensembl IDs to HGNC Symbols
  ensembl_to_gene <- setNames(protein_coding$hgnc_symbol, protein_coding$ensembl_gene_id)
  tpm$gene_symbol <- ensembl_to_gene[rownames(tpm)]
  tpm <- tpm[!is.na(tpm$gene_symbol) & tpm$gene_symbol != "", ]
  tpm <- tpm[!duplicated(tpm$gene_symbol), ]
  
  # === 6. Singscore Analysis ===
  PIPdx <- getGmt(gmtPath)
  rownames(tpm) <- tpm$gene_symbol
  tpm$gene_symbol <- NULL
  
  eranks <- rankGenes(tpm)
  singscores <- multiScore(eranks, PIPdx)$Scores %>% as.data.frame()
  
  # === 7. Merge with Metadata ===
  transposed_singscores <- t(singscores) %>% as.data.frame()
  transposed_singscores$sample_id <- rownames(transposed_singscores)
  
  # Convert to long format for merging
  long_singscores <- transposed_singscores %>%
    pivot_longer(cols = -sample_id, names_to = "Pathway", values_to = "Singscore")
  
  # Merge with metadata
  new_data <- merge(metadata, long_singscores, by = "sample_id")
  new_data$recurrence_status <- factor(new_data$recurrence_status)
  new_data$Timepoint <- factor(new_data$Timepoint, levels = c("Baseline", "Week 6"))
  new_data$Response <- ifelse(new_data$MPRvNMPR == 1, "MPRs", "NMPRs")
  
  # === 8. Load Existing Data & Append ===
  existing_data <- readRDS("merged_sing_df.rds")
  # Only add unique samples to avoid duplicates
  combined_data <- existing_data %>%
    anti_join(new_data, by = "sample_id") %>%
    bind_rows(new_data)

  # === 9. Save Back to RDS ===
  saveRDS(combined_data, "merged_sing_df.rds")
  
  return(combined_data)
}
