context("Checking blocks")

test_that("check_scheme detects centroid", {
  expect_that(check_scheme("centroid"), is_identical_to("centroid"))
  expect_that(check_scheme("CENTROID"), is_identical_to("centroid"))
  expect_that(check_scheme("cen"), is_identical_to("centroid"))
  expect_that(check_scheme("CEN"), is_identical_to("centroid"))
  expect_that(check_scheme("c"), is_identical_to("centroid"))
  expect_that(check_scheme("C"), is_identical_to("centroid"))
})


test_that("check_scheme detects factorial", {
  expect_that(check_scheme("factorial"), is_identical_to("factorial"))
  expect_that(check_scheme("FACTORIAL"), is_identical_to("factorial"))
  expect_that(check_scheme("factor"), is_identical_to("factorial"))
  expect_that(check_scheme("FACTOR"), is_identical_to("factorial"))
  expect_that(check_scheme("f"), is_identical_to("factorial"))
  expect_that(check_scheme("F"), is_identical_to("factorial"))
})

test_that("check_scheme detects path", {
  expect_that(check_scheme("path"), is_identical_to("path"))
  expect_that(check_scheme("PATH"), is_identical_to("path"))
  expect_that(check_scheme("p"), is_identical_to("path"))
  expect_that(check_scheme("P"), is_identical_to("path"))
})

test_that("check_scheme detects bad schemes", {
  expect_warning(check_scheme(1:10), 
                 "Invalid 'scheme'. Default 'scheme=centroid' is used.")
  expect_warning(check_scheme("inner"), 
                 "Invalid 'scheme'. Default 'scheme=centroid' is used.")
  expect_warning(check_scheme("centroide"), 
                 "Invalid 'scheme'. Default 'scheme=centroid' is used.")
})