test_that("shared_neigh counts common neighbors correctly", {
  
  Ne_list <- list(
    c(2),
    c(1,3),
    c(2)
  )
  
  expect_equal(shared_neigh(1, 2, Ne_list), 0)
  
  expect_equal(shared_neigh(2, 3, Ne_list), 0)
  
  expect_equal(shared_neigh(1, 3, Ne_list), 1)
  
  expect_equal(shared_neigh(2, 2, Ne_list), 2)
})
