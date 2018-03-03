context("Binarize")


## set up data sets and predictor matrix parameter for tests
pred_matrix <- (1 - diag(1, 9))
pred_matrix[6,7] <- 0
boys_groups1 <- c(0,0,0,0,0,1,1,0,0)
boys_groups2 <- c(0,0,0,0,0,2,1,0,0)
boys_weights <- c(1,1,1,1,1,4,5,1,1)
boys_bin1 <- mice.binarize(boys_data)
boys_bin2 <- mice.binarize(boys_data, include_observed = TRUE, pred_matrix = pred_matrix)
boys_bin3 <- mice.binarize(boys_data, cols = "gen")
boys_bin4 <- mice.binarize(boys_data, blocks = c("gen","phb"), weights = c(2,3))
boys_bin5 <- mice.binarize(boys_data, blocks = boys_groups2)
boys_bin6 <- mice.binarize(boys_data, blocks = boys_groups1, weights = boys_weights)
boys_bin7 <- mice.binarize(boys_data, weights = boys_weights)


## test input checks
test_that("Test input checks",
          {
            expect_error(mice.binarize(mammal_data), "Data doesn't contain any non-binary categorical attributes with unobserved values.")
            expect_error(mice.binarize(mammal_data, cols = "gt"), "Not all columns in argument 'cols' are non-binary factors.")
            expect_error(mice.binarize(mammal_data, cols = "ct"), "Argument 'cols' contains invalid column names.")
            expect_error(mice.binarize(mammal_data, cols = c(10,20)), "Argument 'cols' contains an invalid column index.")
            expect_error(mice.binarize(boys_data, cols = c("gen","gen")), "Argument 'cols' contains duplicates.")
            expect_error(mice.binarize(boys_data, include_ordered = FALSE), "Data doesn't contain any non-binary unordered categorical attributes with unobserved values.")
            expect_error(mice.binarize(boys_data, pred_matrix = diag(1,9)), "Diagonal elements of input predictor matrix have to be zero.")
            expect_error(mice.binarize(boys_data, weights = c(2,3)), "Input argument 'weights' must not be in list format if no blocks have been specified.")
            expect_error(mice.binarize(boys_data, blocks = boys_groups1, weights = c(2,3,4)), "Length of weight tuple 1 does not match length of block 1.")
            expect_error(mice.binarize(boys_data, blocks = c(0,0,0,0,0,2,3,0,0)), "Argument 'blocks' contains invalid group indices.")
          })


## test functionality
test_that("Test functionality",
          {
            expect_equal_to_reference(boys_bin1, file = "boys_bin1.rds")
            expect_equal_to_reference(boys_bin2, file = "boys_bin2.rds")
            expect_equal_to_reference(boys_bin3, file = "boys_bin3.rds")
            expect_equal_to_reference(boys_bin4, file = "boys_bin4.rds")
            expect_equal_to_reference(boys_bin5, file = "boys_bin5.rds")
            expect_equal_to_reference(boys_bin6, file = "boys_bin6.rds")
            expect_equal_to_reference(boys_bin7, file = "boys_bin7.rds")
          })
