context("Factorize")


## load mids objects and parameter lists for tests
mids_boysbin <- readRDS("mids_boysbin.rds")
post_boysbin <- readRDS("post_boysbin.rds")
par_list1 <- readRDS("par_list1.rds")
par_list2 <- readRDS("par_list2.rds")


## run mice.factoize() for test
res_boysbin <- mice.factorize(post_boysbin$midsobj, par_list2)


## test input checks
test_that("Test input checks",
          {
            expect_error(mice.factorize(mids_boysbin, par_list = par_list1), "Argument 'par_list' is not consistent with data of input mids object.")
            expect_error(mice.factorize(mids_boysbin, par_list = par_list2), "The imputed values in the binarized columns are not in proper format.*")
          })


## test functionality
test_that("Test functionality",
          {
            expect_equal_to_reference(res_boysbin, file = "res_boysbin.rds")
          })
