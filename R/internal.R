###########################################################################################################################################################
# internal.R
#
# INTERNAL R SCRIPT CONTAINING ADDITIONAL HELPER FUNCTIONS FOR THE miceExt PACKAGE
###########################################################################################################################################################

#----------------------------------------------------------------------------------------------------------------------------------
# find_cols
# helper function to find column tuples with identical missing data patterns (or imputation target patterns) in given data.frame
#
# IDEA:
# - get target matrix of all columns of mids$where-matrix that are in visist sequence of input mids-object
# - initalize 'to_visit'-vector of 'unchecked' column indices of target matrix and iterate over it:
#    -> always use the first column in 'to_visit' as reference column
#       and compare it to all other columns in 'to_visit', saving indices of all column that are equal to it,
#    -> obtain a tuple of indices of identical columns and append it to result list
#    -> delete all indices in appended tuple from to_visit
#    -> break iteration when 'to_visit' vector contains less than two elements
#----------------------------------------------------------------------------------------------------------------------------------

find_blocks <- function(obj)
{
  # determine set of candidate columns, taking visit sequence as well as imputation methods into account
  candidates <- intersect(which(obj$method %in% c("pmm","norm")), obj$visitSequence)
  names(candidates) <- names(obj$method)[candidates]
  n_length <- length(candidates)

  if(n_length < 2)
    return(NULL)

  # initalize target matrix
  target_matrix <- obj$where[,candidates]

  # initialize result list
  blocks <- vector(mode="list", length = n_length)

  # initialize helper variables for iteration
  to_visit <- 1:n_length
  blocks_index <- 1

  #iterate
  while (length(to_visit) > 1)
  {
    # get reference column
    first <- to_visit[1]
    fcol <- target_matrix[,first]

    # get indices of all columns that are equal to current reference column
    # first comparison is always redundant, but should not be slower than two extra vector allcations
    equal_cols_ind <- unlist(lapply(to_visit, function(i) if(identical(fcol, target_matrix[,i])) i else 0))

    # append resulting index tuple to output list
    blocks[[blocks_index]] <- candidates[equal_cols_ind]
    blocks_index <- blocks_index + 1

    #delete indices of identical columns from 'to_visit' vector
    to_visit <- setdiff(to_visit,equal_cols_ind)
  }

  # filter resulting list by only returning tuples of length bigger than one
  return(blocks[lapply(blocks,length) > 1])
}



#----------------------------------------------------------------------------------------------------------------------------------
# get_partition
#
# helper function to partition row sets of observed and missing values of y by the values of another column
# that has been specified by match_vars
#
#----------------------------------------------------------------------------------------------------------------------------------
get_partition <- function(target_col, ry, wy)
{
  # divide target column by rows in which y is observed and in which y is to be imputed
  target_obs <- target_col[ry]
  target_mis <- target_col[wy]

  # get unique values in both split set and make sure that all values in the "missing" rows are also contained in the
  # "observed" rows to ensure a proper matching later
  obs_values <- unique(target_obs)
  mis_values <- unique(target_mis)

  if(!all(mis_values %in% obs_values))
    return(NULL)

  # use lapply to iterate over values of "missing" rows only, save for each value the indices in which they are taken
  obs_partition <- lapply(mis_values, function(k) which(target_obs == k))
  mis_partition <- lapply(mis_values, function(k) which(target_mis == k))

  # return both "partitions" (technically obs_partition might not be a real partition as some values of it could be missing)
  return(list(obs_partition = obs_partition, mis_partition = mis_partition))
}



#----------------------------------------------------------------------------------------------------------------------------------
# expand_factors
#
# helper function that expands matrices by replacing factor columns in place with a dummy coded version of themselves
#
#----------------------------------------------------------------------------------------------------------------------------------
expand_factors <- function(x)
{
  if(all(!unlist(lapply(1:ncol(x), FUN = function(j) is.factor(x[[j]])))))
    return(x)
  else
    return(do.call(cbind, lapply(1:ncol(x), 
      function(j)
      {
        col <- x[,j]
        if(!is.factor(col))
          return(col)
        else
          return(model.matrix(~ col - 1)[,-1])
      })))
}



#----------------------------------------------------------------------------------------------------------------------------------
# inflate_vector
#
# helper function that inflates vectors which is needed for transforming predictor matrices in mice.binarize
#
#----------------------------------------------------------------------------------------------------------------------------------
inflate_vector <- function(vec, by)
{
  unlist(lapply(seq_along(vec), function(j) rep(vec[j],by[j])))
}



#-------------------------------------------------------------------------------------------------------------------------------
# call_remove_lin_dep
# hack to access right frame from lower environment in buried function call
#-------------------------------------------------------------------------------------------------------------------------------
call_remove_lin_dep <- function(x, y, ry, minvar, maxcor)
{
  keep <- remove.lindep(x, y, ry, eps = minvar, maxcor = maxcor)
  return(keep)
}
