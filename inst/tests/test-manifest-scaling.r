context("check scaling of factors in data frame")

test_that("check_manifest_scaling works for right scaling", {
  aux =  c(3, 4, 5, 1, 2, 5, 1)
  MV = iris[,aux]
  good_scaling = list(c("num", "num", "nom"),
                      c("raw", "raw"),
                      c("nom", "ord"))
  
  expect_that(check_manifest_scaling(MV, good_scaling), is_true())
  expect_that(get_metric(c("raw", "raw")), is_true())
  expect_that(get_metric(c("num", "raw")), is_true())
})

test_that("check_manifest_scaling throws errors", {  
  aux =  c(3, 4, 5, 1, 2, 5, 1)
  MV = iris[,aux]
  bad_scaling = list(c("num", "num", "ord"),
                     c("raw", "raw"),
                     c("ord", "ord"))
  
  expect_that(check_manifest_scaling(MV, bad_scaling), throws_error())
})