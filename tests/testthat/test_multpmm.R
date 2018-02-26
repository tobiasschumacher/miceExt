context("Multivariate PMM")


## set random seed for reproducability
set.seed(4711)

## create data set to test for match_vars error handling
mammal2 <- mammal_data
mammal2[!is.na(mammal2$sws),"pi"] <- 2L

# minor helper variables for error checks
mammal_blocks_fail1 <- c(0,0,0,1,1,0,0,2,2,0,0,0)
mammal_blocks_fail2 <- c(1,0,0,1,1,0,2,2,0,0,0)

## load other test data
mids_nh <- readRDS("mids_nh.rds")
mids_mammal <- readRDS("mids_mammal.rds")
mids_fail1 <- readRDS("mids_fail1.rds")
mids_fail2 <- readRDS("mids_fail2.rds")


## setup result objects of mice.post.matching to compare those to reference data in test

# general functionality, euclidian metric to avoid issues with eigen()/LAPACK
post_nh1 <- mice.post.matching(mids_nh, distmetric = "euclidian")

# apply match_vars and further parameters
post_nh2 <- mice.post.matching(mids_nh, blocks = c("bmi","hyp"), donors = 2L, distmetric = "manhattan", matchtype = 2L, minvar = 0.0002, weights = c(1,5), match_vars = "age")

# apply blocks and weights using with vector notation
post_nh3 <- mice.post.matching(mids_nh, blocks = c(0,1,1,0), weights = c(1,2,3,4), distmetric = "euclidian")


## test input checks
test_that("Test input checks",
          {
            expect_error(mice.post.matching(mammal_data), "Argument 'obj' should be of type mids.")
            expect_error(mice.post.matching(mids_mammal, blocks = "gt"), "Argument 'blocks' contains a tuple of length 1 with imputation method 'pmm'.")
            expect_error(mice.post.matching(mids_mammal, blocks = "ct"), "Argument 'blocks' contains a tuple with invalid column names.")
            expect_error(mice.post.matching(mids_mammal, blocks = c("sws","gt")), "Not all tuples in given columns are either blockwise NA or blockwise non-NA.")
            expect_error(mice.post.matching(mids_mammal, blocks = c(10,20)), "Argument 'blocks' contains a tuple with an out-of-bounds column index.")
            expect_error(mice.post.matching(mids_mammal, blocks = list(c("sws","ps"), c("sws","ps"))), "Argument 'blocks' contains duplicate columns among its elements.")
            expect_error(mice.post.matching(mids_mammal, blocks = mammal_blocks_fail1), "Argument 'blocks' contains a tuple with an out-of-bounds column index.")
            expect_error(mice.post.matching(mids_mammal, blocks = mammal_blocks_fail2), "Argument 'blocks' contains a tuple with a column index that is not in the visit sequence.")
            expect_error(mice.post.matching(mids_fail1), "There are no column tuples with identical missing data patterns and valid imputation methods.")
            expect_error(mice.post.matching(mids_fail2, blocks = c("sws","ps"), match_vars = "pi"), "Column tuple*")
          })


## test functionality
test_that("Test functionality",
          {
            expect_equal_to_reference(post_nh1, file = "post_nh1.rds")
            expect_equal_to_reference(post_nh2, file = "post_nh2.rds")
            expect_equal_to_reference(post_nh3, file = "post_nh3.rds")
          })
