###########################################################################################################################################################
# check.R
#
# INTERNAL R SCRIPT CONTAINING FUNCTIONS FOR INPUT CHECKS THAT ARE USED AMONG ALL miceExt-FUNCTIONS
#
###########################################################################################################################################################


#==========================================================================================================================================================
# check_blocks
# function dedicated to check whether input list of vector tuples [blocks] in mice.post.matching() is valid
# -> each element of blocks represents a tuple of columns that are contextually intertwined
#
# CHECK CRITERIA:
# - blocks should be either a list of vectors or a single vector, vectors should contain either integral numbers, representing column index,
#   or characters for column names,that actually appear in the data, with no duplicates among them
# - as vectors are contextually intertwined, the row in the mids$where-matrix in those columns shouuld be identical
# - imputation methods on those colum either have to be "pmm" or "norm" to guarantee matching in post-processing that is consistent with previous run of mice
#==========================================================================================================================================================

check_blocks <- function(obj, blocks)
{
  #----------------------------------------------------------------------------------------------------------------
  # First define some helper functions
  #----------------------------------------------------------------------------------------------------------------
  
  
  # helper function that checks whether a single tuple of column names or indices is valid and converts column names to their respective index values
  check_blocks_tuple <- function(tuple)
  {
    #input tuple needs to be atomic
    if(!is.atomic(tuple))
      stop("Argument 'blocks' contains non-atomic element.\n")

    # now check whether input tuple contains valid values w.r.t. names/indices of given data
    #
    # for numeric tuples, cast to integers
    #
    if(is.character(tuple))
    {
      # input tuple shouldn't contain duplicates
      if(anyDuplicated(tuple) > 0)
        stop("Argument 'blocks' contains a tuple with duplicate column names.\n")

      if(!all(tuple %in% varnames))
        stop("Argument 'blocks' contains a tuple with invalid column names.\n")

      # convert character vector to integer vector by mapping column names to their respective index
      tuple <- sapply(tuple, function(i) which(varnames==i), USE.NAMES = FALSE)

    }
    else if(is.numeric(tuple))
    {
      # check whether all column index are finite and not NaN
      if(!all(is.finite(tuple)))
        stop("Argument 'blocks' contains a tuple with NaN or infinite column index.\n")

      # check whether column index are integral
      if(!isTRUE(all.equal(tuple,as.integer(tuple))))
        stop("Argument 'blocks' contains a tuple with a non-integral column index.\n")

      # check whether column indices are valid
      tuple <- as.integer(tuple)
      if(!all((tuple <= nvar) & (tuple > 0L)))
        stop("Argument 'blocks' contains a tuple with an out-of-bounds column index.\n")

      # input tuple shouldn't contain duplicates
      if(anyDuplicated(tuple) > 0)
        stop("Argument 'blocks' contains a tuple with a duplicate column index.\n")

    }
    else
      stop("Argument blocks contains tuple of invalid data type.\n")

    # check whether we are in visit sequence
    if(!all(tuple %in% obj$visitSequence))
      stop("Argument 'blocks' contains a tuple with a column index that is not in the visit sequence.\n")

    ## check whether all imputation methods are valid
    if (!(all(obj$method[tuple] %in% c("pmm","norm","custom"))))
      stop("Argument 'blocks' contains a tuple with an invalid imputation method. The imputation method has to be either 'norm' or 'pmm'.\n")

    if(length(tuple) > 1)
    {
      # now check whether in given column tuple, all rows are either exclusively NA or exclusively non-NA values
      target_matrix <- where[,tuple]
      if(!all(apply(target_matrix, 1, function (row) all(row) | all(!row))))
        stop("Not all tuples in given columns are either blockwise NA or blockwise non-NA.\n")
    }

    ## attach names and return
    names(tuple) <- varnames[tuple]
    return(tuple)
  }
  
  check_blocks_vector_format <- function(groupvec)
  {
    
    # check whether all group numbers are actually finite
    if(!is.numeric(groupvec))
      stop("Argument 'blocks' has to be numeric.\n")
    
    # check whether all column indices are finite and not NaN
    if(!all(is.finite(groupvec) & (groupvec >= 0)))
      stop("Argument 'blocks' contains a group number that is negative, NaN or infinite.\n")
    
    # check whether groupn numbers are integral
    if(!isTRUE(all.equal(unname(groupvec),as.integer(groupvec))))
      stop("Argument 'blocks' contains a non-integral group number.\n")
    
    # check whether group numbers are valid
    groupvec <- as.integer(groupvec)
    n_groups <- max(groupvec)
    if(n_groups < 1)
      stop("Argument 'blocks' must contain a positive group number.\n")
    
    if(!all(1L:max(groupvec) %in% groupvec))
      stop("Argument 'blocks' contains invalid group indices.\n")
    
    
    blocks <- lapply(1L:n_groups, 
      function(j)
      {
        tuple <- which(groupvec == j)
        
        ## check whether we are in visit sequence
        if(!all(tuple %in% obj$visitSequence))
          stop("Argument 'blocks' contains a tuple with a column index that is not in the visit sequence.\n")
        
        ## check whether all imputation methods are valid
        if (!(all(obj$method[tuple] %in% c("pmm","norm","custom"))))
          stop("Argument 'blocks' contains a tuple with an invalid imputation method. The imputation method has to be either 'norm' or 'pmm'.\n")
        
        if(length(tuple) > 1)
        {
          # now check whether in given column tuple, all rows are either exclusively NA or exclusively non-NA values
          target_matrix <- where[,tuple]
          if(!all(apply(target_matrix, 1, function (row) all(row) | all(!row))))
            stop("Not all tuples in given columns are either blockwise NA or blockwise non-NA.\n")
        }
        
        ## attach names and return
        names(tuple) <- varnames[tuple]
        return(tuple)
      })
    
    return(blocks)
  }
  
  
  #----------------------------------------------------------------------------------------------------------------
  # !!! ACTUAL CHECK FUNCTION BEGINS HERE !!!
  #----------------------------------------------------------------------------------------------------------------
  
  # initialize some helper variables
  where <- obj$where
  nvar <- ncol(obj$data)
  varnames <- dimnames(where)[[2]]

  # if blocks is null, look for columns with identical NA distributions
  if(is.null(blocks))
  {
    blocks <- find_blocks(obj)
    if(length(blocks) == 0)
      stop("There are no column tuples with identical missing data patterns and valid imputation methods.\n")

    return(list(blocks = blocks, b_list_format = FALSE))
  }


  # main functionality of check.blocks:
  # if blocks isn't a list, check whether it is either a valid group vector or a "valid" tuple
  # if blocks is a list, check whether all its elements are valid tuples
  # in any way, return a list of valid tuples of indices, not column names
  if(is.atomic(blocks) && length(blocks) == nvar)
  {
    blocks <- check_blocks_vector_format(blocks)
    b_list_format <- FALSE
  }
  else
  {
    # if blocks are atomic, put them into a list
    if(is.atomic(blocks))
      blocks <- list(blocks)
    
    # check every column tuple
    blocks <- lapply(blocks, check_blocks_tuple)
    
    # check whether there are duplicate columns among all tuples
    if(anyDuplicated(unlist(blocks)) > 0)
      stop("Argument 'blocks' contains duplicate columns among its elements.\n")
    
    # set list format bool
    b_list_format <- TRUE
  }

  return(list(blocks = blocks, b_list_format = b_list_format))
}



#==========================================================================================================================================================
# check_weights
# function dedicated to check whether input list of dimension weights [weights] in mice.post.matching() is valid
# -> each element of weights represents a tuple of weights that is related to column tuple of same index in blocks
#
# CHECK CRITERIA:
# - weights should be either a list of vectors or a single vector, vectors should be of the same length as the column tuple of same index,
#   containing only postitive numbers
# - alternatively, an element of the list may be NULL, 0 or 1 which are substitute values that indicates that no weights are to be applied on current
#   column tuple, or weights may be NULL to indicate that no weights should be applied at all
#==========================================================================================================================================================

check_weights <- function(weights, blocks, nvar, b_list_format)
{
  
  #----------------------------------------------------------------------------------------------------------------
  # First define some helper functions
  #----------------------------------------------------------------------------------------------------------------
  
  check_weights_tuple <- function(weights_index)
  {
    # get current element
    curr_weights <- weights[[weights_index]]

    # if current element is null, return
    if(is.null(curr_weights))
      return(NULL)

    # check whether curr_weights are numeric
    if(!is.numeric(curr_weights))
      stop("Argument 'weights' contains non-numeric element.\n")

    # check whteher current value is NULL-substitute and return NULL if so
    if(length(curr_weights) == 1 && curr_weights %in% c(0,1))
      return(NULL)

    # check whether curr_weights vector has same length as corresponding column tuple
    if(length(curr_weights) != length(blocks[[weights_index]]))
      stop("Argument 'weights' contains curr_weights tuple of invalid length.\n")

    # check whether curr_weights are neither NaN nor infinite
    if(!all(is.finite(curr_weights)))
      stop("Argument 'weights' contains a tuple with an element that is either NaN or infinite.\n")

    # check whether all curr_weights are positive
    if(!all(curr_weights > 0))
      stop("Argument 'weights' contains a tuple with a non-positive element.\n")

    return(curr_weights)
  }
  
  convert_weights_vec <- function(weights_vec)
  {
    # check whether weights are numeric
    if(!is.numeric(weights_vec))
      stop("Argument 'weights' is not numeric.\n")
    
    # check whether weights are neither NaN nor infinite
    if(!all(is.finite(weights_vec)))
      stop("Argument 'weights' contains an element that is either NaN or infinite.\n")
    
    # check whether all weights are positive
    if(!all(weights_vec > 0))
      stop("Argument 'weights' contains a tuple with a non-positive element.\n")
    
    if(length(unique(weights_vec[unlist(blocks)])) == 1)
      return(NULL)
    
    weights <- lapply(blocks, 
      function(tuple)
      {
        weights <- weights_vec[tuple]
        if(length(unique(weights)) == 1)
          return(NULL)
        else
          return(weights)
      })
    
    return(weights)
  }

  
  #----------------------------------------------------------------------------------------------------------------
  # !!! ACTUAL CHECK FUNCTION BEGINS HERE !!!
  #----------------------------------------------------------------------------------------------------------------
  
  # check for null, which would be valid
  if(is.null(weights))
    return(weights)

  
  if(is.atomic(weights) && length(weights) == nvar)
    weights <- convert_weights_vec(weights)
  else
  {
    if(!b_list_format)
      stop("Input argument 'weights' must not be in list format if blocks haven't been specified in list format as well.")
    
    # if weights is atomic, make it a list
    if(is.atomic(weights))
      weights <- list(weights)
    
    # check whether weights is of type list
    if(!is.list(weights))
      stop("Argument 'weights' is not atomic or of type list.\n")
    
    # check whether weights list is of same length as blocks
    if(length(weights) != length(blocks))
      stop("The arguments 'weights' and 'blocks' have different lengths.\n")
    
    weights <- lapply(seq_along(weights), check_weights_tuple)
  }
  return(weights)
}



#==========================================================================================================================================================
# check_match_vars
# function dedicated to check whether input argument match_vars of mice.post.matching() is valid
# -> each element of match_vars represents an extra column in the data that is matched against
#
# CHECK CRITERIA:
# - match_vars should contain either integral numbers, representing column indices, or characters for columns names, which should all appear in the data
# - all columns should be either integers or factors to allow a safe split
# - match_vars should be of the same length as blocks, and no element of match_vars should be in the tuple of blocks of the same index
#==========================================================================================================================================================

check_match_vars <- function(obj, blocks, match_vars)
{
  if(is.null(match_vars))
    return(NULL)

  data <- obj$data

  if(is.character(match_vars))
  {
    match_vars[is.null(match_vars)] <- ""

    if(!all(match_vars %in% colnames(data) | match_vars == ""))
      stop("Argument 'match_vars' contains invalid column names.\n")

    # convert character vector to integer vector by mapping column names to their respective index
    match_vars <- unlist(sapply(match_vars,
                         function(i)
                         {
                           if(i != "")
                             return(which(colnames(data)==i))
                           else
                             return(0L)
                         }, USE.NAMES = FALSE))
  }
  else if(is.numeric(match_vars))
  {
    match_vars[is.null(match_vars)] <- 0

    # check whether all column indices are valid
    if(!all(match_vars %in% 1:ncol(data) | match_vars == 0))
      stop("Argument 'match_vars' contains an invalid column index.\n")

    # cast columns to integer
    match_vars <- as.integer(match_vars)

  }
  else
    stop("Argument 'match_vars' is neither numeric nor a character vector.\n")

  if(length(blocks) != length(match_vars))
    stop("Argument 'match_vars' has to be of the same length as argument 'blocks'.\n")

  if(!all(!unlist(lapply(seq_along(match_vars), function(j) match_vars[j] %in% blocks[[j]]))))
    stop("Elements of argument 'match_vars' must not be contained in corresponding element of argument 'blocks'.\n")

  if(!all(unlist(lapply(match_vars, function(j) j == 0 || is.factor(data[,j]) || is.integer(data[,j]) ))))
    stop("Columns specified in argument 'match_vars' must be either factors or integers.\n")

  if(!all(unlist(lapply(match_vars, function(j) all(!is.na(data[,j]))))))
    stop("Columns specified in argument 'match_vars' must not contain any NAs.\n")

  return(match_vars)
}





#==========================================================================================================================================================
# check_optionals
# function dedicated to check whether all optional input arguments of mice.post.matching() next to blocks and weights are valid
#
# CHECK CRITERIA:
# - argument "distmetric" should be a character string in "euclidian", "manhattan", "residual" and "mahalanobis"
# - arguments "donors" and "matchtype" should be integral values, with matchtype between 0 and 2 and donors bigger than 0
# - arguments eps, ridge, maxcor shoupld be numeric values bigger than zero
#==========================================================================================================================================================

check_optionals <- function(optionals, nblocks)
{

  #----------------------------------------------------------------------------------------------------------------
  # First define some helper functions
  #----------------------------------------------------------------------------------------------------------------
  
  check_integral <- function(n,argname)
  {
    # check whether n is numeric
    if(!is.numeric(n))
      stop(paste0("Argument '",argname,"' is not numeric.\n"))

    # check whether n is a single number
    if(length(n) == 1)
    {
      # check whether n is finite and not NaN
      if(!is.finite(n))
        stop(paste0("Argument '",argname,"' is either NaN or infinite.\n"))
      
      # check whether n is integral
      if(!isTRUE(all.equal(n,as.integer(n))))
        stop(paste0("Argument '",argname,"' is not an integer.\n"))
      
      n <- rep(n, nblocks)
    }
    else
    {
      if(length(n) != nblocks)
        stop(paste0("Argument '",argname,"' has to be either a single number or a vector with as many elements as there are blocks.\n"))
      
      # check whether n is finite and not NaN
      if(!is.finite(n))
        stop(paste0("Argument '",argname,"' is either NaN or infinite.\n"))
      
      # check whether n is integral
      if(!isTRUE(all.equal(n,as.integer(n))))
        stop(paste0("Argument '",argname,"' is not an integer.\n"))
    }
    
    return(as.integer(n))
  }

  check_delta <- function(delta, argname)
  {
    # check whether delta is numeric
    if(!is.numeric(delta))
      stop(paste0("Argument '",argname,"' is not numeric.\n"))

    # if delta is a single number, check whether it is valid and expand it
    if(length(delta) == 1)
    {
      # check whether delta is finite and not NaN
      if(!is.finite(delta))
        stop(paste0("Argument '",argname,"' is either NaN or infinite.\n"))
      
      # check whether delta is bigger than 0
      if(delta <= 0)
        stop(paste0("Argument '",argname,"' is smaller than 0.\n"))
      
      delta <- rep(delta, nblocks)
    }
    else
    {
      # if delta is a vector, check whether is has correct length
      if(length(delta) != nblocks)  
        stop(paste0("Argument '",argname,"' has to be either a single number or a vector with as many elements as there are blocks.\n"))
      
      # check whether delta is finite and not NaN
      if(any(!is.finite(delta)))
        stop(paste0("Argument '",argname,"' contains an element that is either NaN or infinite.\n"))
      
      # check whether delta is greater than 0
      if(any(delta <= 0))
        stop(paste0("Argument '",argname,"' contains an element that is smaller than 0.\n"))
    }
    
    return(delta)
  }

  
  #----------------------------------------------------------------------------------------------------------------
  # !!! ACTUAL CHECK FUNCTION BEGINS HERE !!!
  #----------------------------------------------------------------------------------------------------------------
  
  ## check donors

  # check whether donors is integral
  optionals$donors <- check_integral(optionals$donors, "donors")

  # check whether donors is bigger than one
  if(optionals$donors < 1)
    stop("Argument 'donors' is smaller than 1.")


  ## check distmetric

  # check whether distmetric is a character vector
  if(!is.character(optionals$distmetric))
    stop("Argument 'distmetric' is not a character string.\n")

  # check whether distmetric is a single character string
  if(length(optionals$distmetric) == 1)
  {
    # check whether dstfunction is one of the four valid character strings
    if(!(optionals$distmetric %in% c("manhattan", "euclidian", "mahalanobis", "residual")))
      stop("Argument 'distmetric' is invalid. It has to be one of the following: \n \t 'manhattan', 'euclidian', 'mahalanobis', 'residual'.\n")
    
    optionals$distmetric <- rep(optionals$distmetric, nblocks)
  }
  else
  {
    if(length(optionals$distmetric) != nblocks)
      stop("Argument 'distmetric' has to be a single character of a character vector with as many elements as there are blocks.\n")
    
    # check whether dstfunction is one of the four valid character strings
    if(!all(optionals$distmetric %in% c("manhattan", "euclidian", "mahalanobis", "residual")))
      stop("Argument 'distmetric' contains an invalid element. All metrics have to be one of the following: \n \t 'manhattan', 'euclidian', 'mahalanobis', 'residual'.\n")
    
  }
  

  ## check matchtype
  optionals$matchtype <- check_integral(optionals$matchtype, "matchtype")

  # check whether matchtype is between 0 and 2
  if(!all(optionals$matchtype %in% c(0L,1L,2L)))
    stop("Argument 'matchtype' is not an integer between 0 and 2.\n")


  ## check ridge
  optionals$ridge <- check_delta(optionals$ridge, "ridge")

  # check whether ridge is finite and not NaN
  if(any(optionals$ridge) > 1)
    stop("Argument 'ridge' contains an element is bigger than 1.\n")


  ## check minvar
  optionals$minvar <- check_delta(optionals$minvar, "minvar")


  ## check maxcor
  optionals$maxcor <- check_delta(optionals$maxcor ,"maxcor")

  return(optionals)
}



#==========================================================================================================================================================
# check_deep
#
# Function dedicated to perform deeper checks on whether the arguments "blocks" and "match_vars" of mice.post.matching() work with the input data/mids object,
# which also returns partitions of both missing and observed data based on match_vars
#
# More precisely, we need to check for two scenarios:
# 1. When collecting the predictive means y_obs and y_mis for the matching step, we can only use those values whose designated predictors have
#    actually been completely observed or imputed. As we collect those predictors column-wise, we need to intersect between all those predictors, which may
#    result in an empty set of observed or missing values that we want to match on.
# 2. In case that we match against an external variable, we need to make sure that all the values that the external variable takes on the rows of missing
#    data are also taken in the rows of observed data, as otherwise there is nothing to match against. To perform this check, we actually need to compute
#    partitions of observed and missing data based on those external values, which are also needed later in the matching step
#    -> Hence, we return those partitions so we do not need to recompute them later
#
# NOTE: This function introduces some redundancy within the structure mice.post.matching, as identical iterations/computations are performed
#       again in the main function. The whole functionality of this check could also be implemented within the main function, but this would have
#       the consequence that such errors might be discovered at a late stage of the whole computation, possibly causing the user to run the program for
#       a long time before an error is detected.
#==========================================================================================================================================================

check_deep <- function(obj, blocks, match_vars)
{

  data <- obj$data
  nrows <- nrow(data)
  r <- !is.na(data)
  where <- obj$where

  complete_R <- vector(mode="list", length = length(match_vars))
  complete_W <- vector(mode="list", length = length(match_vars))

  if(is.null(match_vars))
    partitions_list <- NULL
  else
    partitions_list <- vector(mode="list", length = length(match_vars))


  # build current data set consisting of observed and imputed data from current imputation
  for (j in obj$visitSequence)
  {
    wy <- where[, j]
    ry <- r[, j]

    # copy from obj$imp to avoid using filled in values from multivariate match
    data[(!ry) & wy, j] <- obj$imp[[j]][(!ry)[wy], 1]
  }

  # iterate over all column tuples
  for(i in seq_along(blocks))
  {
    # grab current column tuple and correspoding weights
    tuple <- blocks[[i]]

    # keep track of in which rows the predictor values are not complete, as we cannot use those in multivariate matching
    complete_ry <- !vector(mode="logical", length = nrows)
    complete_wy <- !vector(mode="logical", length = nrows)

    # impute every column in current tuple for each imputation and collect y_hats
    for(j in tuple)
    {

      # get predictor matrix
      predictors <- obj$predictorMatrix[j, ] == 1
      x <- data[, predictors, drop = FALSE]
      x <- expand_factors(x)
      
      # filter down to rows that actually have a complete non-NA set of predictors
      ry <- complete.cases(x) & r[, j]
      wy <- complete.cases(x) & where[, j]

      # update shared filter
      complete_ry <- complete_ry & ry
      complete_wy <- complete_wy & wy

      # check whether there are no common predictors left
      if(all(!complete_wy) || all(!complete_ry))
      {
        stop(paste0("There are either no common donors or no common predictors in column tuple (", toString(tuple), ").\n"))
      }

    } # end of tuple loop

    complete_R[[i]] <- complete_ry
    complete_W[[i]] <- complete_wy

    # build partition if there is external variable to match against
    match_col <- match_vars[i]
    if(!is.null(match_col) && match_col != 0L)
    {
      match_partitions <- get_partition(data[,match_col], complete_ry, complete_wy)
      if(is.null(match_partitions))
        stop(paste0("Column tuple (", toString(tuple), ") has to be matched against the values in column ", match_col,
                    ", but some of the missing rows in this tuple have no matching values in that column.\n"))

      partitions_list[[i]] <- match_partitions
    }
    else
    {
      # nothing to partition against
      # -> if input tuple has length 1, input method has to be "norm" cause otherwise there is nothing to do
      if(length(tuple) == 1 && obj$method[tuple] == "pmm")
        stop("Argument 'blocks' contains a tuple of length 1 with imputation method 'pmm' and no external column to match against.\n")
    }

  } #end blocks iteration

  return(list(partitions_list = partitions_list, complete_R = complete_R, complete_W = complete_W))
}


#==========================================================================================================================================================
# check_tuples_weights_binarize
# function dedicated to check whether input arguments cols, blocks, weights, include_ordered, include_observed
# of mice.binarize() are valid. Also convert tules and weights into vector format if necessary and returns vector with 
# indices of those facator columns that are to be binarized
# -> cols represents a set of factor columns that have to be binarized
#
# CHECK CRITERIA:
# - each element should be a integer in within range of ncol(data) or a character represeting a column name of the data
# - each column index should correspond to a factor column
# - there should be no duplicates
#==========================================================================================================================================================

check_blocks_weights_binarize <- function(cols, blocks, weights, include_ordered, include_observed, data)
{
  
  #----------------------------------------------------------------------------------------------------------------
  # First define some helper functions
  #----------------------------------------------------------------------------------------------------------------
  
  
  # helper function that checks whether a single tuple of column names or indices is valid and converts column names to their respective index values
  check_tuple <- function(tuple)
  {
    #input tuple needs to be atomic
    if(!is.atomic(tuple))
      stop("Argument 'blocks' contains non-atomic element.\n")
    
    # now check whether input tuple contains valid values w.r.t. names/indices of given data
    #
    # for numeric tuples, cast to integers
    #
    if(is.character(tuple))
    {
      if(!all(tuple %in% varnames))
        stop("Argument 'cols' contains a tuple with invalid column names.\n")
      
      # input tuple shouldn't contain duplicates
      if(anyDuplicated(tuple) > 0)
        stop("Argument 'cols' contains a tuple with duplicate column names.\n")
      
      # convert character vector to integer vector by mapping column names to their respective index
      tuple <- sapply(tuple, function(i) which(varnames==i), USE.NAMES = FALSE)
      
    }
    else if(is.numeric(tuple))
    {
      # check whether all column indices are finite and not NaN
      if(!all(is.finite(tuple)))
        stop("Argument 'cols' contains a tuple with a NaN or infinite column index.\n")
      
      # check whether column indices are integral
      if(!isTRUE(all.equal(tuple,as.integer(tuple))))
        stop("Argument 'cols' contains a tuple with a non-integral column index.\n")
      
      # check whether column indices are valid
      tuple <- as.integer(tuple)
      if(!all((tuple <= nvar) & (tuple > 0)))
        stop("Argument 'cols' contains a tuple with an out-of-bounds column index.\n")
      
      # input tuple shouldn't contain duplicates
      if(anyDuplicated(tuple) > 0)
        stop("Argument 'cols' contains a tuple with duplicate column indices.\n")
      
    }
    else
      stop("Argument cols contains tuple of invalid data type.\n")
   
    if(length(tuple) > 1)
    {
      # now check whether in given column tuple, all rows are either exclusively NA or exclusively non-NA values
      target_matrix <- is.na(data[,tuple])
      if(!all(apply(target_matrix, 1, function (row) all(row) | all(!row))))
        warning("Not all tuples in given columns are either blockwise NA or blockwise non-NA.\n")
    }
    
    ## attach names and return
    names(tuple) <- varnames[tuple]
    return(tuple)
  }
  
  check_blocks_vector_format <- function(groupvec)
  {
    
    # check whether all group vector is actually finite
    if(!is.numeric(groupvec))
      stop("Argument 'blocks' has to be numeric.\n")
    
    # check whether group vector has valid length
    if(length(groupvec) != nvar)
      stop("Argument 'blocks' has invalid length.\n")
    
    # check whether all column index are finite and not NaN
    if(!all(is.finite(groupvec) & groupvec >= 0))
      stop("Argument 'blocks' contains a tuple index that is negative, NaN or infinite.\n")
    
    # check whether column index are integral
    if(!isTRUE(all.equal(groupvec,as.integer(groupvec))))
      stop("Argument 'blocks' contains a non-integral group number.\n")
    
    # check whether group numbers are valid
    groupvec <- as.integer(groupvec)
    n_groups <- max(groupvec)
    if(n_groups < 1)
      stop("Argument 'blocks' must contain a positive group number.\n")
    
    if(!all(1L:max(groupvec) %in% groupvec))
      stop("Argument 'blocks' contains invalid group indices.\n")
    
    src_factor_cols <- unlist(lapply(1L:n_groups, 
      function(j)
      {
        tuple <- which(blocks == j)
        if(length(tuple) > 1)
        {
          target_matrix <- is.na(data[,tuple])
          if(!all(apply(target_matrix, 1, function (row) all(row) | all(!row))))
            warning("Not all column tuples in given set of blocks are either blockwise NA or blockwise non-NA.\n")
        }
        
        factor_inds <- unlist(lapply(tuple, function(j) return(nlevels(data[,j]) > 2)))
        return(tuple[factor_inds])  
      }))
    
    if(length(src_factor_cols) == 0)
      stop("None of the column tuples specifed in argument 'col_groups' contain a non-binary factor variable.")
    
    return(list(groupvec = groupvec, src_factor_cols = src_factor_cols))
  }
  
  
  check_weights <- function(weights_index)
  {
    # get current element
    curr_weights <- weights[[weights_index]]
    
    # if current element is null, return
    if(is.null(curr_weights))
      return(rep(1, length(blocks[[weights_index]])))
    
    # check whether curr_weights are numeric
    if(!is.numeric(curr_weights))
      stop("Argument 'weights' contains non-numeric element.\n")
    
    # check whteher current value is NULL-substitute and return NULL if so
    if(length(curr_weights) == 1 && curr_weights %in% c(0,1))
      return(rep(1, length(blocks[[weights_index]])))
    
    # check whether curr_weights vector has same length as corresponding column tuple
    if(length(curr_weights) != length(blocks[[weights_index]]))
      stop("Argument 'weights' contains curr_weights tuple of invalid length.\n")
    
    # check whether curr_weights are neither NaN nor infinite
    if(!all(is.finite(curr_weights)))
      stop("Argument 'weights' contains a tuple with an element that is either NaN or infinite.\n")
    
    # check whether all curr_weights are positive
    if(!all(curr_weights > 0))
      stop("Argument 'weights' contains a tuple with a non-positive element.\n")
    
    return(curr_weights)
  }
  
  check_weights_vec <- function(weights_vec)
  {
    # check whether weights are numeric
    if(!is.numeric(weights_vec))
      stop("Argument 'weights' is not numeric.\n")
    
    # check whether weights vector has correct length
    if(length(weights_vec) > nvar)
      stop("Input weights vector has too many elements.\n")
    
    # check whether weights are neither NaN nor infinite
    if(!all(is.finite(weights_vec)))
      stop("Argument 'weights' contains an element that is either NaN or infinite.\n")
    
    # check whether all weights are positive
    if(!all(weights_vec > 0))
      stop("Argument 'weights' contains a tuple with a non-positive element.\n")

    return(weights_vec)
  }
  
  
  #----------------------------------------------------------------------------------------------------------------
  # !!! ACTUAL CHECK FUNCTION BEGINS HERE !!!
  #----------------------------------------------------------------------------------------------------------------
  
  nvar <- ncol(data)
  varnames <- colnames(data)
  
  
  if(!is.null(blocks))
  {
    # check column format first
    
    if(!is.list(blocks) && length(blocks) == nvar)
    {
      check_res <- check_blocks_vector_format(blocks)
      res_blocks <- check_res$groupvec
      src_factor_cols <- check_res$src_factor_cols
    }
    else
    {
      
      res_blocks <- rep(0, nvar)
      factor_inds <- rep(FALSE, nvar)
      if(class(blocks) != "list")
        blocks <- list(blocks)
      
      for (i in seq_along(blocks))
      {
        curr_tuple <- check_tuple(blocks[[i]])
        blocks[[i]] <- curr_tuple       ## need to overwrite this here for valid duplicate check at later stage
        factor_inds[curr_tuple] <- unlist(lapply(curr_tuple, function(j) return(nlevels(data[,j]) > 2)))
        res_blocks[curr_tuple] <- i
      }
      
      # check whether there are duplicate columns among all tuples
      if(anyDuplicated(unlist(blocks)) > 0)
        stop("Argument 'blocks' contains duplicate columns among its elements.\n")
      
      
      ## check whether there actually is a non-binary factor in given set of target columns
      ## -> has to be at the end
      if(!any(factor_inds))
        stop("None of the column tuples specifed in argument 'col_groups' contain a non-binary factor variable.")
      
      src_factor_cols <- (1L:nvar)[factor_inds]
    }
  }
  else 
  {
    res_blocks <- NULL
    
    if(!is.null(cols))
    {
      # in this case check whether input tuple contains valid values w.r.t. names/indices of given data
      #
      # for numeric tuples, cast to integers
      #
      if(is.character(cols))
      {
        # input cols shouldn't contain duplicates
        if(anyDuplicated(cols) > 0)
          stop("Argument 'cols' contains duplicates.\n")
        
        if(!all(cols %in% colnames(data)))
          stop("Argument 'cols' contains invalid column names.\n")
        
        # convert character vector to integer vector by mapping column names to their respective index
        cols <- sapply(cols, function(i) which(colnames(data)==i), USE.NAMES = FALSE)
      }
      else if(is.numeric(cols))
      {
        # check whether column index are valid
        if(!all(cols %in% 1:ncol(data)))
          stop("Argument 'cols' contains an invalid column index.\n")
        
        # cast columns to integer
        cols <- as.integer(cols)
        
        # cols shouldn't contain duplicates
        if(anyDuplicated(cols) > 0)
          stop("Argument 'cols' contains a duplicate column index.\n")
        
      }
      else
        stop("Argument 'cols' is neither numeric nor a character vector.\n")
      
      
      if(!all(unlist(lapply(cols, function(j) is.factor(data[,j]) && nlevels(data[,j]) > 2))))
        stop("Not all columns in argument 'cols' are non-binary factors.\n")
      
      src_factor_cols <- cols
    }
    else
    {
      src_factor_ind <- unlist(lapply(1L:nvar, 
        function (j)
        {
          col <- data[,j]
          if(nlevels(col) > 2 && (include_observed || any(is.na(col))) && (include_ordered || !is.ordered(col)) )
            return(TRUE)
          else
            return(FALSE)
        }))
      
      if(!any(src_factor_ind))
      {
        if(include_observed)
        {
          if(include_ordered)
            stop("Data doesn't contain any non-binary categorical attributes.\n")
          else
            stop("Data doesn't contain any non-binary unordered categorical attributes.\n")
        }
        else
        {
          if(include_ordered)
            stop("Data doesn't contain any non-binary categorical attributes with unobserved values.\n")
          else
            stop("Data doesn't contain any non-binary unordered categorical attributes with unobserved values.\n")
        }
        
      }
      src_factor_cols <- (1L:nvar)[src_factor_ind]
    }
  }
  
  
  ### CHECK WEIGHTS HERE
   
  if(is.null(weights))
    res_weights <- rep(1, nvar)
  else
  {
    # if weights are in vector format, delegate to vector format check
    if(is.atomic(weights) && length(weights) == nvar)
    {
      res_weights <- check_weights_vec(weights)
    }
    else
    {
      # if weights are in atomic vector format but not of length ncol(data), we interpret it as single eement weights list
      if(is.atomic(weights))
        weights <- list(weights)
      
      # weights in list format are only allowed, if blocks have been specified in list format as well 
      if(is.null(blocks))
        stop("Input argument 'weights' must not be in list format if no blocks have been specified.\n")
      
      # weights in list format are only allowed, if blocks have been specified in list format as well 
      if(!is.list(blocks))
        stop("Input argument 'weights' must not be in list format if blocks haven't been specified in list format as well.\n")

      # check whether all weights are valid
      check_weights_vec(unlist(weights))
      
      # check whether weights is of type list
      if(!is.list(weights))
        stop("Argument 'weights' is not atomic or of type list.\n")
      
      # check whether weights list is of same length as cols
      if(length(weights) != length(blocks))
        stop("The arguments 'weights' and 'blocks' have different lengths.\n")
      
      res_weights <- rep(1, nvar)
      
      for (j in seq_along(blocks))
      {
        tuple <- blocks[[j]]
        curr_weights <- check_weights(j)
        
        if(length(tuple) != length(curr_weights))
          stop(paste0("Length of weight tuple ",j," does not match length of block ",j,".\n"))
        
        res_weights[tuple] <- curr_weights
      }
      
    }    
  }
  
  setup <- list(blocks = res_blocks, weights = res_weights, src_factor_cols = src_factor_cols)
  return(setup)
}




#==========================================================================================================================================================
# check_pred_matrix
# function dedicated to check whether input argument pred_matrix of mice.binarize() is valid
# -> matrix should be usable as predictorMatrix argument in mice().
#
# CHECK CRITERIA:
# - predictor matrix has to be quadratic with as many columns as the original data
# - all entries should be either 0, 1 or 2
# - diagonal entires should be 0
#==========================================================================================================================================================

check_pred_matrix <- function(pred_matrix, n)
{
  # check for matrix type
  if(!is.matrix(pred_matrix))
    stop("Argument 'pred_matrix' has to be matrix.\n")

  # check matrix size
  if(nrow(pred_matrix) != n || ncol(pred_matrix) != n)
    stop("Input predictor matrix does not have the correct size.\n")

  # check values
  if(!all(pred_matrix %in% c(0,1,2)))
    stop("Input predictor matrix contains invalid values.\n")

  # check diagonal entries
  if(!all(diag(pred_matrix) == 0))
    stop("Diagonal elements of input predictor matrix have to be zero.\n")

}



#==========================================================================================================================================================
# check_par_list
# function dedicated to check whether input argument par:list of mice.factorize is valid
# -> par_list is a list of parameters containing information of the original data as well as some transformation parameters
# -> elements are src_data [the original data], n_src_cols [#cols of originial data], src_factor_cols [factor cols in original data that have been binarized],
#    dummy_cols [indices of binary columns in transformed data], src_levels [levels of transformed factors], and src_names [original names of transformed columns]
# -> we have to check whether par_list fits input data, i.e. whtether all elements in par_list are valid by type, fit to data and are also
#
# CHECK CRITERIA:
# - par_list should be a named list with exactly the elements that are mentioned above
# - each element should have a valid type and fit to data [e.g. all columns in dummy_cols should indded be binary]
# - elements should be consistent with each other [e.g. number of tuples in dummy_cols should equal number of columns in src_factor_cols]
#==========================================================================================================================================================

check_par_list <- function(obj, par_list)
{

  ## first check whether par_list is of the correct format

  #check if par_list is not NULL
  if(is.null(par_list))
    stop("Argument 'par_list' is NULL.\n")

  # check if par_list actually is a list
  if(!is.list(par_list))
    stop("Argument 'par_list' has to be a list.\n")

  # check whether par_list only contains the elements that it is supposed to contain
  if(!identical(sort(names(par_list)), c("dummy_cols", "n_pad_cols", "n_src_cols", "pad_names", "src_data", "src_factor_cols", "src_levels", "src_names" )))
    stop("Argument 'par_list' has to exclusively contain the elements 'src_data', 'n_src_cols', 'n_pad_cols', 'src_factor_cols', 'dummy_cols', 'src_names', 'pad_names', 'src_levels'.\n")

  # create local copy of each element
  src_data <- par_list$src_data
  n_src_cols <- par_list$n_src_cols
  n_pad_cols <- par_list$n_pad_cols
  src_factor_cols <- par_list$src_factor_cols
  dummy_cols <- par_list$dummy_cols
  src_levels <- par_list$src_levels
  src_names <- par_list$src_names
  pad_names <- par_list$pad_names


  ## now make a deeper check whether all elements contain values of valid types

  # begin with src_data
  if(!is.data.frame(par_list$src_data))
    stop("Element 'src_data' in argument 'par_list' has to be a data frame.\n")

  # n_src_cols
  if(!is.numeric(n_src_cols) || length(n_src_cols) != 1 || !is.finite(n_src_cols) || !isTRUE(all.equal(n_src_cols,as.integer(n_src_cols))) || n_src_cols < 2)
    stop("Element 'n_src_cols' of argument 'par_list' is not an integral number bigger than one.\n")

  # src_factor_cols
  if(!is.numeric(src_factor_cols) || !all(is.finite(src_factor_cols)) || !isTRUE(all.equal(src_factor_cols,as.integer(src_factor_cols))) || !all(src_factor_cols > 0))
    stop("Element 'src_factor_cols' of argument 'par_list' is either NULL or contains invalid values.\n")

  # dummy_cols
  if(!is.list(dummy_cols))
    stop("Element 'dummy_cols' of argument 'par_list' is not a list.\n")

  if(!all(unlist(lapply(dummy_cols, is.numeric))))
    stop("Element 'dummy_cols' of argument 'par_list' contains non-numeric elements.\n")

  ul_dummy_cols <- unlist(dummy_cols)
  if(!all(is.finite(ul_dummy_cols)) || !isTRUE(all.equal(ul_dummy_cols,as.integer(ul_dummy_cols))) || !all(ul_dummy_cols > 0))
    stop("Element 'dummy_cols' of argument 'par_list' contains invalid values.\n")


  # src_levels
  if(!is.list(src_levels))
    stop("Element 'src_levels' of argument 'par_list' is not a list.\n")

  if(!all(unlist(lapply(src_levels, is.character))))
    stop("Element 'src_levels' of argument 'par_list' contains non-numeric elements.\n")


  # src_names
  if(!is.character(src_names))
    stop("Element 'src_names' of argument 'par_list' is either NULL or not a character vector.\n")



  ## now check whether values of all elements are consistent with each other and, i.e., with obj$data

  # src_factor_cols must not be bigger than n_src_cols
  if(max(src_factor_cols) > n_src_cols)
    stop("Element 'src_factor_cols' of argument 'par_list' contains an out-of-bounds value.\n")

  # src_factor_sols and dummy_cols have to be same length
  if(length(src_factor_cols) != length(dummy_cols))
    stop("Elements 'src_factor_cols' and 'dummy_cols' of argument 'par_list' are not of the same length.\n")

  # dummy_cols has to be consistent with src_levels
  if(length(dummy_cols) != length(src_levels) || !identical(unlist(lapply(dummy_cols, length)), unlist(lapply(src_levels, length))))
    stop("Elements 'dummy_cols' and 'src_levels' of argument 'par_list' are not consistent with each other.\n")

  # src_factor_cols must not be bigger than n_src_cols
  if(length(src_names) != n_src_cols)
    stop("'par_list$src_names' doesn't have par_list$n_src_cols arguments.\n")

  # check whether we have consistency with src_data
  if(ncol(src_data) != n_src_cols)
    stop("'par_list$src_data' doesn't have par_list$n_src_cols columns.\n")

  if(!identical(names(src_data), src_names))
    stop("'par_list$src_data' doesn't have par_list$n_src_cols columns.\n")

  if(!all(unlist(lapply(seq_along(src_factor_cols), function(j) is.factor(src_data[,src_factor_cols[j]]) && identical(levels(src_data[,src_factor_cols[j]]), src_levels[[j]] )))))
    stop("Elements 'src_factor_cols', 'src_levels' and 'src_data' of argument 'par_list' are not consistent with each other.\n")



  ## finally, check whether values of all elements are consistent with obj$data

  # dummy_cols has to fit to obj$data
  if(n_pad_cols != ncol(obj$data))
    stop("Argument 'par_list' is not consistent with data of input mids object.\n")

  if(!identical(pad_names, names(obj$data)))
    stop("Argument 'par_list' is not consistent with data of input mids object.\n")

  if(!all(as.matrix(obj$data[,ul_dummy_cols]) %in% c(0,1, NA)))
    stop("Element 'dummy_cols' of argument 'par_list' is not consistent with data of input mids object.\n")

  # make sure all colmns in dummy_cols actually are binary and include one non-zero entry at max
  if(!all(unlist(lapply(dummy_cols,
    function(tuple)
    {
      all(apply(obj$data[,tuple], MARGIN = 1,
            function(row)
            {
              row <- row[!is.na(row)]
              return (all(is.na(row)) || (all(row %in% c(0,1)) && sum(row) == 1))
            }))
    }))))
    stop("Not every column tuple in given list of dummmy columnss is in proper binarized format.\n")

}
