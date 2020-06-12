## Tests for linear_trend() & related

## load packages
library("testthat")

context("Test linear_trend() & related")

test_that("linear_trend() returns a tibble", {
    expect_silent(df <- linear_trend(t = 1:20))
    expect_s3_class(df, "linear_trend")
    expect_named(df, c("t", "trend"))
    expect_identical(nrow(df), 20L)
})

test_that("simulate_linear_trend() returns a tibble", {
    expect_silent(df <- simulate_linear_trend(t = 1:20, seed = 42))
    expect_s3_class(df, "simulate_linear_trend")
    expect_named(df, c("t", "trend", "y"))
    expect_identical(nrow(df), 20L)
})
