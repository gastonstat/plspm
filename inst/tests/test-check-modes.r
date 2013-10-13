context("Checking modes")

test_that("check_modes works as expected", {
  modA = c("A", "A", "A")
  modB = c("B", "B", "B")
  modnewA = c("newA", "newA", "newA")
  modAB = c("A", "B", "B")
  modAnewA = c("A", "newA", "A")
  modBnewA = c("newA", "B", "B")  
  modABnewA = c("A", "B", "newA")
  block3 = list(1:2, 3:4, 5:6)
  
  expect_that(check_modes(NULL, block3), is_identical_to(modA))
  expect_that(check_modes(modA, block3), is_identical_to(modA))
  expect_that(check_modes(modB, block3), is_identical_to(modB))
  expect_that(check_modes(modnewA, block3), is_identical_to(modnewA))
  expect_that(check_modes(modAB, block3), is_identical_to(modAB))
  expect_that(check_modes(modAnewA, block3), is_identical_to(modAnewA))
  expect_that(check_modes(modABnewA, block3), is_identical_to(modABnewA))
})

test_that("check_mdoes detects bad modes", {
  modA2 = c("A", "A")
  modAAS = c("A", "A", "S")  
  block3 = list(1:2, 3:4, 5:6)
  
  expect_that(check_modes(modA2, block3), gives_warning())
  expect_error(check_modes(modAAS, block3), "Sorry. Don't know how to work with mode: 'S'")
})