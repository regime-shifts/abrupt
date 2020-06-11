## Tests for step_trend() & related

## load packages
library("testthat")

context("Test step_trend() & related")

test_that("step_trend() returns a tibble", {
    expect_silent(df <- step_trend(t = 1:20, change_points = c(5, 15),
                                   means = c(1, 2, 3)))
    expect_s3_class(df, "step_trend")
    expect_named(df, c("t", "trend"))
    expect_identical(nrow(df), 20L)
})

test_that("simulate_step_trend() returns a tibble", {
    expect_silent(df <- simulate_step_trend(t = 1:20, change_points = c(5, 15),
                                   means = c(1, 2, 3), seed = 42))
    expect_s3_class(df, "simulate_step_trend")
    expect_named(df, c("t", "trend", "y"))
    expect_identical(nrow(df), 20L)
})
