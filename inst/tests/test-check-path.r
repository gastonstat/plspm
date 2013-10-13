context("Checking path matrix")

test_that("check_path works as expected with matrices", {
  some_path = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
  rownames(some_path) = c("LV1", "LV2", "LV3")
  colnames(some_path) = c("LV1", "LV2", "LV3")
  
  expect_that(check_path(some_path), is_identical_to(some_path))
})

test_that("check_path detects bad path matrices", {
  bad1 = as.matrix(0)
  bad2 = matrix(1:12, 4, 3)
  bad3 = matrix(1:9, 3, 3)
  bad4 = matrix(c(0,0,0,2,0,0,0,-1,0), 3, byrow=TRUE)
  
  expect_error(check_path(1:10), "'path_matrix' must be a matrix.")
  expect_error(check_path("string"), "'path_matrix' must be a matrix.")
  expect_error(check_path(bad1), "'path_matrix' must have more than one element.")
  expect_error(check_path(bad2), "'path_matrix' must be a square matrix.")
  expect_error(check_path(bad3), "'path_matrix' must be a lower triangular matrix.")
  expect_error(check_path(bad4), "Elements in 'path_matrix' must be '1' or '0'")
})