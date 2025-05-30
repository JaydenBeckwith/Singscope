library(testthat)
library(dplyr)

get_mutation_summary <- function(df) {
  df %>%
    filter(!is.na(Mutation)) %>%
    mutate(
      PatientID = sub("_S[0-9]+$", "", sample_id),
      Mutation = as.character(Mutation)
    ) %>%
    distinct(PatientID, Mutation) %>%
    group_by(Mutation) %>%
    summarise(Count = n(), .groups = "drop")
}

test_that("Mutation summary aggregates correctly", {
  # Test input
  test_df <- data.frame(
    sample_id = c("P1_S0", "P1_S1", "P2_S0", "P3_S0", "P4_S0"),
    Mutation = c("BRAFV600E", "BRAFV600E", "WT", "NRAS", NA),
    stringsAsFactors = FALSE
  )

  expected <- data.frame(
    Mutation = c("BRAFV600E", "NRAS", "WT"),
    Count = c(1, 1, 1)
  )

  result <- get_mutation_summary(test_df)

  expect_equal(result[order(result$Mutation), ], expected[order(expected$Mutation), ])
})