context("Multivariate PMM")


## set random seed for reproducability
set.seed(4711)

## create data set to test for match_vars error handling
mammal2 <- mammal_data
mammal2[!is.na(mammal2$sws),"pi"] <- 2L

## load other test data
mids_nh <- readRDS("mids_nh.rds")
mids_mammal <- readRDS("mids_mammal.rds")
mids_fail1 <- readRDS("mids_fail1.rds")
mids_fail2 <- readRDS("mids_fail2.rds")


## setup result objects of mice.post.matching to compare those to reference data in test

# general functionality, euclidian metric to avoid issues with eigen()/LAPACK
post_nh1 <- mice.post.matching(mids_nh, distmetric = "euclidian")

# apply match_vars and further parameters
post_nh2 <- mice.post.matching(mids_nh, donors = 2L, distmetric = "manhattan", matchtype = 2L, eps = 0.0002, weights_list = c(1,5), match_vars = "age")


## test input checks
test_that("Test input checks",
          {
            expect_error(mice.post.matching(mammal_data), "Argument 'obj' should be of type mids.")
            expect_error(mice.post.matching(mids_mammal, cols = "gt"), "Argument 'cols' contains a tuple of length 1 with imputation method 'pmm'.")
            expect_error(mice.post.matching(mids_mammal, cols = "ct"), "Argument 'cols' contains a tuple with invalid column names.")
            expect_error(mice.post.matching(mids_mammal, cols = c("sws","gt")), "Not all tuples in given columns are either blockwise NA or blockwise non-NA.")
            expect_error(mice.post.matching(mids_mammal, cols = c(10,20)), "Argument 'cols' contains a tuple with out-of-bounds column numbers.")
            expect_error(mice.post.matching(mids_mammal, cols = list(c("sws","ps"), c("sws","ps"))), "Argument 'cols' contains duplicate columns among its elements.")
            expect_error(mice.post.matching(mids_fail1), "There are no column tuples with identical missing data patterns and valid imputation methods.")
            expect_error(mice.post.matching(mids_fail2, cols = c("sws","ps"), match_vars = "pi"), "Column tuple*")
          })


## test functionality
test_that("Test functionality",
          {
            expect_equal_to_reference(post_nh1, file = "post_nh1.rds")
            expect_equal_to_reference(post_nh2, file = "post_nh2.rds")
          })
