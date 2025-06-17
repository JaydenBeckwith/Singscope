source("../../R/survival.R")

#load dummy test data 
test_data <- read.csv("test_data/mock_survival_test_data.csv")
test_data <- as.data.frame(test_data)

test_that("Surv test1: OS survival analysis includes only valid patients", {
    result <- run_survival_analysis(test_data, survival_type = "OS", group_col = "group", time_unit = "days")
    expect_true(all(result$data$eligible_OS))
    expect_true(all(!is.na(result$data$time)))
    expect_true(all(result$data$time >= 0))
})

test_that("Surv test2: EFS survival analysis only assigns time for eligible patients", {
  result <- run_survival_analysis(test_data, survival_type = "EFS", group_col = "group", time_unit = "days")
  df <- result$data
  print(df)
  # Only eligible patients should have a non-NA time value
  expect_true(all(is.na(df$time[!df$eligible_EFS])))
  expect_true(all(!is.na(df$time[df$eligible_EFS])))
  
  # All eligible patients should have time >= 0
  expect_true(all(df$time[df$eligible_EFS] >= 0))
})


test_that("Surv test3: MSS survival analysis only assigns time for eligible patients", {
  result <- run_survival_analysis(test_data, survival_type = "MSS", group_col = "group", time_unit = "days")
  df <- result$data
  print(df)
  # Only eligible patients should have a non-NA time value
  expect_true(all(is.na(df$time[!df$eligible_MSS])))
  expect_true(all(!is.na(df$time[df$eligible_MSS])))

  # All eligible patients should have time >= 0
  expect_true(all(df$time[df$eligible_MSS] >= 0))
})

test_that("Surv test4: EFS survival time is correctly calculated for a known patient", {
  result <- run_survival_analysis(test_data, survival_type = "EFS", group_col = "group", time_unit = "days")
  df <- result$data
  
  # Pick a known patient ID that is eligible for EFS and has known dates
  known_patient <- df[df$patient_id == "PT2", ]
  print(known_patient)
  # Ensure the patient is eligible
  expect_true(known_patient$eligible_EFS)

  # Manually calculate the expected time difference
  expected_time <- as.numeric(as.Date("2020-07-13") - as.Date("2020-01-04"))

  # Check that the calculated time matches
  expect_equal(known_patient$time, expected_time)
})

test_that("survtest 5: Survival analysis runs only if there are at least two groups with non-NA time", {
  df <- test_data
  result <- run_survival_analysis(df, survival_type = "OS", group_col = "group", time_unit = "days")
  
  non_na_groups <- unique(result$data$group[!is.na(result$data$time)])
  
  if (length(non_na_groups) > 1) {
    expect_error(run_survival_analysis(df, survival_type = "OS", group_col = "group", time_unit = "days"), NA)
  } else {
    expect_error(run_survival_analysis(df, survival_type = "OS", group_col = "group", time_unit = "days"),
                 regexp = "There is only 1 group")
  }
})