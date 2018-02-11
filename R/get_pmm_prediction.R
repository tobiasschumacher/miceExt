###########################################################################################################################################################
# get_pmm_prediction
#
# Internal function that computes univariate predictive means which are later used for multivariate matching
#
# Parameters:
# - y:  Vector to be imputed
# - ry: Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the linear model is fitted.
# - x:  Numeric design matrix with length(y) rows with predictors for y
# - wy: Logical vector of length \code{length(y)} indicating locations in y for which imputations are created.
# - matchtype: Type of matching distance
# - ridge: The ridge penalty used in internal call of mice::norm.draw().
#
# Return value: Vector containing the predictive means y_hat for the missing as well as the observed data
#
# Details:
# Code is pretty much identical to mice.impute.pmm, only that in this case we do not match but exclusively build a column of y_hats that will
# be matched against in multivariate context later. More precisely, we distinguish between missing and observed values in y, build a Bayesian
# linear regression model based on the observed data and its corresponding predictors in x, and then apply this model, considering the given
# matchtype,on the predictors of the missing as well as observed values of y to obtain the values y_hat_obs and y_hat_mis.
#
############################################################################################################################################################

get_pmm_prediction <- function(y, ry, x, wy, matchtype, ridge)
{

  # if y is categoric, it is either binary or a factor that the user missed to binarize,
  # so reduce to integer level -1 to get to binary standard
  ynum <- y
  if (is.factor(y))
    ynum <- as.integer(y) - 1

  # add constant to model
  x <- cbind(1, as.matrix(x))

  # perform linear regression and get parameters
  parm <- mice::norm.draw(ynum, ry, x, ridge = ridge)

  # compute y_hats depending on matching type
  if (matchtype == 0L) {
    y_hat_obs <- x[ry, , drop = FALSE] %*% parm$coef
    y_hat_mis <- x[wy, , drop = FALSE] %*% parm$coef
  }
  if (matchtype == 1L) {
    y_hat_obs <- x[ry, , drop = FALSE] %*% parm$coef
    y_hat_mis <- x[wy, , drop = FALSE] %*% parm$beta
  }
  if (matchtype == 2L) {
    y_hat_obs <- x[ry, , drop = FALSE] %*% parm$beta
    y_hat_mis <- x[wy, , drop = FALSE] %*% parm$beta
  }

  # return y_hats
  y_hat <- ynum
  y_hat[ry] <- y_hat_obs
  y_hat[wy] <- y_hat_mis

  return(y_hat)
}
