# Load libraries
library(dplyr)
library(tidyr)
library(readr)

# === Load the data ===
singscores <- read.csv("data/neopele_neotrioalone_singscores.csv")
metadata <- read.csv("data/neopele_neotrioalone_meta.csv")

# === Data Transformation ===
# Set the row names to the pathway names
rownames(singscores) <- singscores$X
singscores <- singscores %>% select(-X)

# ✅ Fix column names to match metadata
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
# ✅ Add MPR and NMPR labels
merged_sing_df <- merged_sing_df %>%
  mutate(Response = ifelse(MPRvNMPR == 1, "MPRs", "NMPRs"))

merged_sing_df

# Save for the Shiny app
saveRDS(merged_sing_df, "merged_sing_df.rds")

