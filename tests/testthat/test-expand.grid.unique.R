test_that("expand.grid.unique works correctly", {
  x <- c(1, 2, 3)
  y <- c(2, 3, 4)
  
  result1 <- expand.grid.unique(x, y, include.equals = FALSE)
  
  # result must be a dataframe and have only 2 columns 
  expect_true(is.data.frame(result1))
  expect_equal(ncol(result1), 2)
  
  # because include.equals = FALSE, there will be no rows with equal elements
  expect_true(all(result1[,1] != result1[,2]))
  
  # include.equals = TRUE allows equal values
  result2 <- expand.grid.unique(x, y, include.equals = TRUE)
  expect_true(any(result2[,1] == result2[,2]))
})
