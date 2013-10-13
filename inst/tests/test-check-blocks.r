context("Checking blocks")

test_that("check_blocks works as expected", {
  num_block1 = list(1:2, 3:4)
  chr_block1 = list(c("A", "B"), c("C", "D"))
  num_block2 = list(1, c(2,1), 3:5, c(4,6,3))
  chr_block2 = list("A", c("B", "A"), c("C", "D", "E"), c("D","F","C"))
  data = matrix(1:24, 4, 6)
  colnames(data) = c("A", "B", "C", "D", "E", "F")
  
  expect_that(check_blocks(num_block1, data), is_identical_to(num_block1))
  expect_that(check_blocks(num_block2, data), is_identical_to(num_block2))
  expect_that(check_blocks(chr_block1, data), is_identical_to(num_block1))
  expect_that(check_blocks(chr_block2, data), is_equivalent_to(num_block2))
})

test_that("check_blocks detects bad blocks", {
  dupli1 = list(c(1,1), 1:2, 4:6)
  dupli2 = list(c("a", "a"), c("b", "c"))
  mixed = list(c(1,2,3), c("a", "b", "c"))
  bad1 = list(1:2, 3:4, 5:6)
  bad2 = list(c("A", "B"), c("D", "C", "E"))
  data = matrix(1:12, 3, 4)
  colnames(data) = c("A", "B", "C", "D")
  
  expect_error(check_blocks(1:10), "'blocks' must be a list.")
  expect_error(check_blocks("string"), "'blocks' must be a list.")
  expect_error(check_blocks(dupli1), 
               "Invalid 'blocks'. Duplicated variables within a block are not allowed")
  expect_error(check_blocks(dupli2), 
               "Invalid 'blocks'. Duplicated variables within a block are not allowed")
  expect_error(check_blocks(mixed), "All elements in 'blocks' must be of same mode")
  expect_error(check_blocks(bad1, data), 
               "Invalid 'blocks'. Indices outside the number of columns in 'Data'")
  expect_error(check_blocks(bad2, data), 
               "Unrecognized name in 'blocks': 'E'")
})