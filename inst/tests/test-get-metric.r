context("get metric")

test_that("get_metric works", {  
  expect_that(get_metric(NULL), is_true())
  expect_that(get_metric(c("ord", "ord")), is_false())
  expect_that(get_metric(c("num", "raw")), is_false())
})
