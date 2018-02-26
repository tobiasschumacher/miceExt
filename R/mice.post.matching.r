#'Post-processing of Imputed Data by Multivariate Predictive Mean Matching
#'
#'Performs multivariate predictive mean matching (PMM) on a set of columns that
#'have been imputed on by the functionalities of the \code{mice}-package.
#'Also offers functionality to match imputations against observed variables.
#'
#'@param obj \code{mice::mids} object that has been returned by a previous run
#'  of \code{mice()} or \code{mice.mids()} and whose imputations we want to
#'  post-process.
#'@param blocks Column tuple or list of column tuples that multivariate PMM has to
#'  be performed on. Each column tuple has to be represented by an atomic vector
#'  that can contain column names as character strings or column indices in
#'  numerical format. The column list must not contain any duplicates, and each
#'  column in the list has to be included in the \code{visitSequence} of the
#'  given mids object. The imputation method that was applied on each column has
#'  to be either \code{pmm} or \code{norm}. \cr
#'  Univariate tuples are also allowed if they have been imputed on via
#'  \code{mice.impute.norm()}, or if we want to match them against on observed
#'  variable as specified in parameter \code{match_vars}. \cr
#'  The default is \code{blocks = NULL}, in which case this function automatically
#'  looks for tuples that have identical missing value patterns and imputes on
#'  those.
#'@param donors Integer indicating the size of the donor pool among which a draw
#'  in the matching step is made. The default is \code{donors = 5L}. Setting
#'  \code{donors = 1L} always selects the closest match (nearest neighbor), but
#'  is not recommended.
#'@param weights List of numeric vectors that allocates weights to the
#'  elements in each tuple in blocks, giving us the possibility to punish or
#'  mitigate differences in certain columns of the data when computing the
#'  distances in the matching step to determine the donor pool. Hence, this list
#'  has to be the same length as blocks, and each element has to be the same
#'  length as the corresponding column tuple. If we do not want to apply weights
#'  to a single tuple, we can write a single \code{0}, \code{1} or \code{NULL}
#'  in the corresponding spot of the list. The default is \code{weights =
#'  NULL}, meaning that no weights should be applied to any columns at all.
#'@param distmetric Character string that determines which mathematical metric
#'  we want to use to compute the distances between the multivariate
#'  \code{y_obs} and \code{y_mis}. Options are \code{"euclidian"},
#'  \code{"manhattan"}, \code{"mahalanobis"} and \code{"residual"}. The first
#'  three options refer to the distance metrics of the same name, the latter is
#'  a variant of the mahalanobis distance in which we consider the residual
#'  covariance \code{cov(y_hat_obs - y_obs)} of the predicted model. This
#'  distance has been proposed and recommended by Little (1988) and is also the
#'  default. Note that the first two options are faster though, as they do not
#'  involve any matrix computations.
#'@param matchtype Integer indicating the type of matching distance to be used
#'  in PMM. The default choice (\code{matchtype = 1L}) calculates the distance
#'  between the \emph{predicted} values of \code{y_obs} and the \emph{drawn}
#'  values of \code{y_mis} (called type-1 matching). Other choices are
#'  \code{matchtype = 0L} distance between predicted values) and \code{matchtype
#'  = 2L} (distance between drawn values).
#'@param match_vars Vector specifying for each tuple which additional variable
#'  has to be matched against. Can be an integer or character vector, either
#'  specifying column indices or column names. \code{match_vars} must be the
#'  same length as \code{blocks} with \code{0}-elements or \code{""}-elements to
#'  disable this functionality for certain groups. Default is \code{match_vars =
#'  NULL}, which completely disables this functionality.
#'@param ridge The ridge penalty used in an internal call of
#'  \code{mice:::.norm.draw()} to prevent problems with multicollinearity. The
#'  default is \code{ridge = 1e-05}, which means that 0.01 percent of the
#'  diagonal is added to the cross-product. Larger ridges may result in more
#'  biased estimates. For highly noisy data (e.g. many junk variables), set
#'  \code{ridge = 1e-06} or even lower to reduce bias. For highly collinear
#'  data, set \code{ridge = 1e-04} or higher.
#'@param eps The minimum variance that we allow predictors to have when building
#'  a linear model to compute the \code{y_hat}-values. More precisely, this
#'  parameter is used in an internal call of \code{mice:::remove.lindep()}. The
#'  default is \code{eps = 1e-04}.
#'@param maxcor The maximum correlation that we allow predictors to have when
#'  building a linear model to compute the \code{y_hat}-values. More precisely,
#'  this parameter is used in an internal call of
#'  \code{mice:::remove.lindep()}.The default is \code{maxcor = 0.99}.
#'
#'@return List containing the following two elements:
#'\describe{
#'  \item{midsobj}{\code{mice::mids} object that differs from the input object
#'               only in the imputations that have been post-processed, and the
#'               \code{call} and \code{loggedEvents} attributes that have been
#'               updated. In particular, those post-processed imputations are
#'               not affecting the \code{chainMean} or
#'               \code{chainVar}-attributes, and hence, \code{plot()} will not
#'               consider them either.}
#'  \item{blocks}{List of column tuples that multivariate imputation has been
#'              performed on. It is equal to the input parameter \code{blocks} if
#'              it has been specified by the user, otherwise those column tuples
#'              have been determined internally.}
#'}
#'
#'@author Tobias Schumacher, Philipp Gaffert
#'
#'@details The algorithm basically iterates over the \code{m} imputations of the
#'  input \code{mice::mids} object and each column tuple in \code{blocks}, each
#'  time following the following two main steps:
#'\describe{
#'  \item{Prediction}{First, we iterate over the columns of the current tuple,
#'  collecting the complete column \eqn{y} in which we want to impute, along
#'  with the matrix \eqn{x} of its predictors. Among the values of y, we
#'  identify the observed values \eqn{y_{obs}} and the missing values
#'  \code{y_{mis}} along with their corresponting predictors \eqn{x_{obs}} and
#'  \code{x_{mis}}, restricting ourselves to only those values whose predictors
#'  are complete, i.e. don't contain any missing values themselves. Note that
#'  among all the columns in the current tuple, there may be no common
#'  predictors in which case the function breaks. \cr
#'  From the observed values and their predictors, we build a Bayesian linear
#'  regression model and, depending on the specified matching type, use the
#'  model's coefficients and drawn regression weights \eqn{\beta} to compute the
#'  predicted means \eqn{\hat{y}_{obs}} and \eqn{\hat{y}_{mis}} which we then
#'  save in a matrix \eqn{\hat{Y}}.}
#'  \item{Matching}{After iterating over the columns in the current tuple and
#'  building the matrix \eqn{\hat{Y}}, we match the multivariate predictions of
#'  the missing values in \eqn{\hat{Y}} against the predictions of the observed
#'  values. More precisely, for each \eqn{\hat{y}_{mis}}, we perform a
#'  k-nearest-neighbor search among all values \eqn{\hat{y}_{obs}}, where
#'  \eqn{k} is the number of donors specified by the user, and then randomly
#'  sample one element of those \eqn{k} nearest neighbors which is then used to
#'  impute the missing value tuple in \eqn{y}. The distance metric that is used to
#'  determine the nearest neighbors is specified by the user as well, and in
#'  case of Euclidian or Manhattan distance its computation is very
#'  straightforward. If the Mahalanobis distance or the residual distance (as
#'  proposed by Little) has been selected, we first compute the corresponding
#'  covariance matrix and its (pseudo) inverse via its eigen decomposition, and
#'  then use it to transform the values \eqn{\hat{y}_{mis}} and
#'  \eqn{\hat{y}_{obs}} that are fed into the nearest-neighbor search. If
#'  weights have been specified as well, they are also applied before the
#'  kNN-search. \cr
#'  If an additional observed variable to match against has been specified via
#'  \code{match_vars}, the set of predictors and missing values are partitioned
#'  by the values of the external variable first, and then the matching is
#'  performed within each pair of corresponding subsets of that partition.}
#'
#'  Note that the imputed values are only stored as a result and, other than in
#'  the \code{mice}-algorithm, are not used to compute predictive means for any
#'  other missing value. We exclusively use the imputed values that are provided
#'  within the input \code{mids}-object.
#'}
#'
#'@references Little, R.J.A. (1988), Missing data adjustments in large surveys
#'  (with discussion), Journal of Business Economics and Statistics, 6,
#'  287--301.
#'
#'  Van Buuren, S. (2012). Flexible Imputation of Missing Data. CRC/Chapman \&
#'  Hall, Boca Raton, FL.
#'
#'  Van Buuren, S., Groothuis-Oudshoorn, K. (2011). \code{mice}: Multivariate
#'  Imputation by Chained Equations in \code{R}. \emph{Journal of Statistical
#'  Software}, \bold{45}(3), 1-67. \url{http://www.jstatsoft.org/v45/i03/}
#'
#'@examples
#'
#'\dontrun{
#'#------------------------------------------------------------------------------
#'# example on modified 'mammalsleep' data set from mice, that has identical
#'# missing data patterns on the column tuples ('ps','sws') and ('mls','gt')
#'#------------------------------------------------------------------------------
#'
#'# run mice on data set 'mammal_data' and obtain a mids object to post-process
#'mids_mammal <- mice(mammal_data)
#'
#'
#'# run function, as blocks has not been specified, it will automatically detect
#'# the column tuples with identical missing data patterns and then impute on
#'# these
#'post_mammal <- mice.post.matching(mids_mammal)
#'
#'#read out which column tuples have been imputed on
#'post_mammal$blocks
#'
#'#look into imputations within resulting mice::mids object
#'post_mammal$midsobj$imp
#'
#'
#'#------------------------------------------------------------------------------
#'# example on original 'mammalsleep' data set from mice, in which we
#'# want to post-process the imputations in column 'sws' by only imputing values
#'# from rows whose value in 'pi' matches the value of 'pi' in the row we impute
#'# on.
#'#------------------------------------------------------------------------------
#'
#'# run mice on data set 'mammal_data' and obtain a mids object to post-process
#'mids_mammal <- mice(mammalsleep)
#'
#'
#'# run function, specify 'sws' as the column to impute on, and specify 'pi' as
#'# the observed variable to consider in the matching.
#'post_mammal <- mice.post.matching(mids_mammal, blocks = "sws", match_vars = "pi")
#'
#'
#'#look into imputations within resulting mice::mids object
#'post_mammal$midsobj$imp
#'
#'
#'#------------------------------------------------------------------------------
#'# example working on 'boys_data', illustrating some more
#'#   advanced functionalities
#'#------------------------------------------------------------------------------
#'
#'# run mice() first
#'mids_boys <- mice(boys_data)
#'
#'# in boys_data, the columns 'hgt' and 'bmi' have the same NA-pattern,
#'# and they have both been imputed on by PMM within mice()
#'# hence, we can perform our post-processing on those columns:
#'blocks <- c("hgt","bmi")
#'
#'
#'# run post-processing, here with some non-default parameters:
#'# -> a reduced donor pool of only 3 donors is chosen
#'# -> weights: we want to punish bigger differnces in bmi within the
#'#    matching, so we apply weight of 3 to this column while leaving the
#'#    other column as is
#'# -> distmetric: within our matching process, we want to take the variances of
#'#    each column into account by normaizing the column values over their
#'#    variance. This is done by using the Mahalanobis-metric.
#'
#'post_boys <- mice.post.matching(mids_boys,
#'                                blocks = blocks,
#'                                donors = 3L,
#'                                weights = c(1, 3),
#'                                distmetric = "mahalanobis")
#'
#'
#'#------------------------------------------------------------------------------
#'# example that illustrates the combined functionalities of mice.binarize,
#'# mice.factorize and mice.post.matching on the dataset 'boys' from mice, to
#'# impute the factor columns 'gen', 'phb' and 'reg'.
#'#------------------------------------------------------------------------------
#'
#'# binarize all factor columns in boys_data that contain NAs
#'boys_bin <- mice.binarize(boys)
#'
#'# run mice on binarized data, note that we need to use boys_bin$data to grab
#'# the actual binarized data and that we use the output predictor matrix
#'# boys_bin$pred_matrix which is recommended for obtaining better imputation
#'# models
#'mids_boys <- mice(boys_bin$data, predictorMatrix = boys_bin$pred_matrix)
#'
#'# it is very likely that mice imputed multiple ones among one set of dummy
#'# variables, so we need to post-process
#'post_boys <- mice.post.matching(mids_boys, distmetric = "residual")
#'
#'# now we can safely retransform to the original data, with non-binarized
#'# imputations
#'res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)
#'
#'# analyze the distribution of imputed variables, e.g. of the column 'gen',
#'# using the mice version of with()
#'with(res_boys, table(gen))
#'
#'
#'#------------------------------------------------------------------------------
#'# similar example to the previous working on 'boys_data' again, in which the
#'# factor columns  'gen' and 'phb' have the same missing data pattern
#'#------------------------------------------------------------------------------
#'
#'# binarize
#'boys_bin <- mice.binarize(boys_data, blocks = c("gen", "phb"))
#'
#'# run mice on binarized data
#'mids_boys <- mice(boys_bin$data, predictorMatrix = boys_bin$pred_matrix)
#'
#'# again we want to post-process the binarized columns. As the binarized columns
#'# of 'gen' and 'phb' share the same missing data pattern, they would, by default,
#'# be imputed as one big tuple, as this is how they are detected internally.
#'# If we want to post-process the binary columns of 'gen' and 'phb' separately, we
#'# can use the element 'dummy_cols' from 'boys_bin$par_list' which was generated
#'# in the binarization, and feed it into the 'blocks'-argument.
#'post_boys <- mice.post.matching(mids_boys, blocks = boys_bin$par_list$dummy_cols)
#'
#'# retransform to the original format
#'res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)
#'}
#'
#'
#'@seealso \code{\link[mice]{mice}}, \code{\link[mice]{mids-class}}, \code{\link[miceExt]{mice.binarize}}, \code{\link[miceExt]{mice.factorize}}
#'@export
mice.post.matching <- function(obj, blocks = NULL, donors = 5L, weights = NULL, distmetric = "residual", matchtype = 1L, match_vars = NULL, ridge = 1e-05, eps = 1e-04, maxcor = 0.99)
{

  ##  Check whether all input arguments are valid

  # check whether input mids object is valid
  if (!is.mids(obj)) stop("Argument 'obj' should be of type mids.")

	# check whether input set of column tuples is valid
  blocks <- check_blocks(obj, blocks)

  # check whether match_vars are valid
  match_vars <- check_match_vars(obj, blocks, match_vars)

  # check whether input list/vector of weights is valid and convert it to list format if necessary
  weights <- check_weights(weights, blocks, ncol(obj$data))
  
  # check validity of other optional parameters
  optionals <- list(donors = donors, distmetric = distmetric, matchtype = matchtype, ridge = ridge , minvar = minvar, maxcor = maxcor)
  optionals <- check_optionals(optionals)

  # extract possibly converted integral values
  donors <- optionals$donors
  matchtype <- optionals$matchtype


  ## Initialize base parameters

  # get size of data frame
  data <- obj$data
  ncols <- ncol(data)
  nrows <- nrow(data)

  # grab function call
  call <- match.call()


  # initialize more base parameters, e.g. make explicit local copies
  where <- obj$where
  r <- !is.na(data)

  # now perform a deep check on whether given blocks and match_vars work on the given data, and retrieve return setup values
  setup <- check_deep(obj, blocks, match_vars)

  partitions_list <- setup$partitions_list
  complete_R <- setup$complete_R
  complete_W <- setup$complete_W

	# initialize list of resulting imputations
	res_imp <- obj$imp

  # grab event log which is needed in some subfunctions of mice
  loggedEvents <- obj$loggedEvents
  state <- list(it = 0, im = 0, co = 0, dep = "", meth = "",
                log = !is.null(loggedEvents))
  if (is.null(loggedEvents))
    loggedEvents <- data.frame(it = 0, im = 0, co = 0, dep = "",
                               meth = "", out = "")



  ######################################################################################################
	#  CENTRAL ITERATION
  ######################################################################################################

  # iterate over each set of imputations
  for(m in 1:obj$m)
  {
    # build current data set consisting of observed and imputed data from current imputation
    for (j in obj$visitSequence)
    {
      wy <- where[, j]
      ry <- r[, j]

      # copy from obj$imp to avoid using filled in values from multivariate match
      data[(!ry) & wy, j] <- obj$imp[[j]][(!ry)[wy], m]
    }

    # iterate over all column tuples
	  for(i in seq_along(blocks))
	  {
	    # grab current column tuple and correspoding weights
	    tuple <- blocks[[i]]
	    dimweights <- weights[[i]]

	    # if tuple is univariate, just perform standard PMM and skip to next tuple
	    if(length(tuple) == 1)
	    {
	      # get current y column
	      y <- data[, tuple]

	      # get current predictors
	      predictors <- obj$predictorMatrix[tuple, ] == 1
	      x <- data[, predictors, drop = FALSE]
	      x <- expand_factors(x)

	      # filter down to rows that actually have a complete non-NA set of predictors
	      ry <- complete.cases(x) & r[, tuple]
	      wy <- complete.cases(x) & where[, tuple]

	      # remove linear dependent predictors
	      keep <- call_remove_lin_dep(x, y, ry, minvar = minvar, maxcor = maxcor)
	      x <- x[, keep, drop = FALSE]

	      # grab current partition and split matching if necessary
	      match_partitions <- partitions_list[[i]]

	      # perform standard pmm procedure with custom univariate predictions and matching
	      y_hat <- get_pmm_prediction(y, ry, x, wy = wy, matchtype = matchtype, ridge = ridge)
	      res_imp[[tuple]][, m] <- match_univariate(y = as.matrix(y), y_hat = as.matrix(y_hat), ry = ry, wy = wy, donors = donors, match_partitions = match_partitions)

	      next
	    }

	    # build matrix of predictive means
	    y_hat_matrix <- do.call(cbind, lapply(tuple,
	      function(j)
	      {
	        # get current y column
	        y <- data[, j]

	        # get predictor matrix
	        predictors <- obj$predictorMatrix[j, ] == 1
	        x <- data[, predictors, drop = FALSE]
	        x <- expand_factors(x)

	        # filter down to rows that actually have a complete non-NA set of predictors
	        ry <- complete.cases(x) & r[, j]
	        wy <- complete.cases(x) & where[, j]

	        # remove linear dependent predictors
	        keep <- call_remove_lin_dep(x, y, ry, minvar = minvar, maxcor = maxcor)
	        x <- x[, keep, drop = FALSE]

	        # apply current imputation/regression method to obain y_hat
	        return(get_pmm_prediction(y, ry, x, wy = wy, matchtype = matchtype, ridge = ridge))
	      }))

      # prepare multivariate matching, get whole data on current tuple
	    y <- data[, tuple]

      # grab rows of common complete observed and missing data per tuple
	    complete_ry <- complete_R[[i]]
	    complete_wy <- complete_W[[i]]

	    #grab current partition
	    match_partitions <- partitions_list[[i]]

	    # now perform multivariate PMM in current tuple
	    imputes <- match_multivariate(y = y, y_hat = y_hat_matrix, ry = complete_ry, wy = complete_wy, donors = donors, distmetric = distmetric, dimweights = dimweights, match_partitions = match_partitions )

	    # fill in imputed values
	    for (j in seq_along(tuple))
	    {
	      res_imp[[tuple[j]]][, m] <- imputes[,j]
	    }
	  }

  }

  # logging stuff
  if (!state$log)
    loggedEvents <- NULL
  if (state$log)
    row.names(loggedEvents) <- seq_len(nrow(loggedEvents))

	## save, and return
  midsobj <- list(call = call,
                  data = obj$data,
                  where = obj$where,
                  m = obj$m,
                  nmis = obj$nmis,
                  imp = res_imp,
                  method = obj$method,
                  predictorMatrix = obj$predictorMatrix,
                  visitSequence = obj$visitSequence,
                  pad = obj$pad,
                  post = obj$post,
                  seed = obj$seed,
                  iteration = obj$iteration,
                  lastSeedValue = obj$lastSeedValue,
                  chainMean = obj$chainMean,
                  chainVar = obj$chainVar,
                  loggedEvents = loggedEvents)

  oldClass(midsobj) <- "mids"

  return(list(midsobj = midsobj, blocks = blocks))
}
