context("Checking Data")

test_that("check_data works as expected with matrices", {
  some_matrix = matrix(1:10, 5, 2)
  good_matrix = some_matrix
  rownames(good_matrix) = 1:5
  colnames(good_matrix) = c("MV1", "MV2")

  expect_that(check_data(some_matrix), is_identical_to(good_matrix))
})

test_that("check_data works as expected with data frames", {
  good_df = iris[c(1:5,51:55,101:105),]
  
  expect_that(check_data(good_df), is_identical_to(good_df))
})

test_that("check_data detects bad data", {
  bad_data = list(1:10)
  bad_matrix = matrix(letters[1:15], 5, 3)
  one_row = iris[1,]
  one_column = as.matrix(1:5)
  
  expect_that(check_data(1:10), throws_error())
  expect_that(check_data("string"), throws_error())
  expect_that(check_data(bad_data), throws_error())
  expect_that(check_data(bad_matrix), throws_error())
  expect_that(check_data(one_row), throws_error())
  expect_that(check_data(one_column), throws_error())
})