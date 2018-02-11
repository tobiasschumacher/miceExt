###########################################################################################################################################################
# match_multivariate
#
# Internal function that matches precomputed predictive means of multivariate missing data against the predictive means of the observed data
#
# Parameters:
# - y:  Matrix that is to be imputed on.
# - y_hat: Matrix containing all predictive means.
# - ry: Logical vector of length length(y) indicating the rows of observed data that are matched against.
# - wy: Logical vector of length length(y) indicating the rows in y for which imputations are created.
# - donors: Size of the donor pool that is computed and drawn from.
# - distmetric: Character string specifying which distance metric is used when computing the donor pool.
# - dimweights: Numerical vector of weights that are applied on the columns of the predictive mean matrices when computing the donor pool.
# - match_partition: List containing partitions of y[ry] and y[wy] that are based on the values of an additional column that has to be matched by which
#                    has been specified externally by the variable match_vars.
#
# Return value: Matrix containing the imputed values for the missing data.
#
# Details:
# Matching procedure as proposed by Little (1988). When computing the distances, we use the RANN package to perform a fast k-nearest-neighbor
# search via kd-trees. If we choose Mahalanobis distance or its analogon that uses the residual covariance (which is recommended by Little) to
# compute the donor pool, we first have to transform the y_hats as we cannot feed quadratic forms into RANN::nn2. To do this, we first compute
# the square root of the (pseudo)inverse of the covariance matrix via its eigendecomposition, taking advantage of its symmetric positive
# semidefinite structure. We then multiply this squareroot onto the y_hats to obtain the desired transformation and feed the results into the
# nearest-neighbor search, after applying weights in case they have been specified. In the end, we randomly sample from the resulting donor pool
# and return the drawn values as the imputations.
#
############################################################################################################################################################

match_multivariate <- function(y, y_hat, ry, wy, donors, distmetric, dimweights, match_partitions)
{

  #if a column is categoric, it has to be binary, so reduce to integer level -1 to get to binary standard
  ynum <- do.call(cbind, lapply(1:ncol(y), FUN =
    function(i)
    {
      if(is.factor(y[,i]))
        as.integer(y[,i]) -1
      else
        y[,i]
    }))

  #initialize y_obs and y_hat_obs/mis
  y_obs <- ynum[ry,]
  y_hat_obs <- y_hat[ry,]
  y_hat_mis <- y_hat[wy,]

  #if number of donors is bigger than or equal to number of observed values, we can sample immediately
  if(donors >= sum(ry))
  {
    message("Specified number of donors is bigger than number of observed values that donors can be taken from.")
    idx <- replicate(sum(wy), sample(1:sum(ry),1))
    return(idx)
  }

  ## if we want to compute distance via a covariance matrix, we have some more work to do:
  ## we want to use RANN for fastest nearest neighbor search which only provides euclidian distance,
  ## hence we have to multiply y with the square root of the inverse to get the same effect
  if(!distmetric %in% c("manhattan", "euclidian"))
  {
    #compute fitting covariance matrix
    if(distmetric == "mahalanobis")
      cov_matrix <- cov(y_obs)
    else
      cov_matrix  <- cov(y_obs - y_hat_obs)

    ## compute inverse of covariance matrix via eigen decomposition
    ## if matrix isn't invertible, we can use and simply compute the pseudo inverse from this decomposition as well

    # get eigen decomposition of covariance matrix
    cov_eigen <- eigen(cov_matrix, symmetric = TRUE)

    # check which eigenvalues are sufficiently bigger than 0 so we can invert them
    inv_ind <- cov_eigen$values > (cov_eigen$values[[1]] * .Machine$double.eps)

    #inverted eigenvalues > eps and get their square root
    ev_inv_sqrt <- rep(0, ncol(y))
    ev_inv_sqrt[inv_ind] <- 1/sqrt(cov_eigen$values[inv_ind])

    # compute distance matrix = (pseudo) inverse of covariance matrix via eigen decomposition
    dist_matrix <- cov_eigen$vectors %*% diag(ev_inv_sqrt) %*% t(cov_eigen$vectors)

    # transform y_hats
    y_hat_obs <-  y_hat_obs %*% dist_matrix
    y_hat_mis <-  y_hat_mis %*% dist_matrix
  }

  ## apply weights
  if(!is.null(dimweights))
  {
    y_hat_obs <- t(t(y_hat_obs) * dimweights) # '*' operator multiplies component-wise periodically in first dimension, so we need to transpose
    y_hat_mis <- t(t(y_hat_mis) * dimweights) # to get column-wise multiplication
  }

  #=================================================================================================================================================
  # DEPRECATED: NORMALIZATION AND JITTER
  # normalize y_hat matrices by median norm of shared rows
  #med_norm <- median(apply(rbind(y_hat_mis, y_hat_obs), MARGIN = 1, function(x) sqrt(sum(x^2))))
  #y_hat_mis <- med_norm * y_hat_mis
  #y_hat_obs <- med_norm * y_hat_obs

  # add random epsilon error on normailized y_hat_mis to simulate a random tie break in RANN
  #eps_matrix <- matrix(runif(nrow(y_hat_mis) * ncol(y_hat_mis), min = -1e-15, max = 1e-15), nrow = nrow(y_hat_mis) , ncol = ncol(y_hat_mis) )
  #y_hat_mis <- y_hat_mis + eps_matrix
  #=================================================================================================================================================


  # check whether we also match by external data
  # -> if yes, match on partitioned data
  if(is.null(match_partitions))
  {
    # now do nearest neighbor search, either by manhattan or euclidian norm
    if(distmetric == "manhattan")
    {
      nearest_neighbors <- RANN.L1::nn2(y_hat_obs, query = y_hat_mis, k = donors)
    }
    else
    {
      nearest_neighbors <- RANN::nn2(y_hat_obs, query = y_hat_mis, k = donors)
    }

    # sample from the donor pool of each missing observation
    idx <- apply(nearest_neighbors$nn.idx, 1, sample, size = 1)
  }
  else
  {
    # prepare index of soon-to-be imputed rows and current partitions of missing and observed data
    idx <- 1:sum(wy)
    obs_partition <- match_partitions$obs_partition
    mis_partition <- match_partitions$mis_partition

    for(i in seq_along(obs_partition))
    {

      # if number of donors is bigger than or equal to number of observed values, we can sample immediately
      if(donors >= length(obs_partition[[i]]))
      {
        idx[mis_partition[[i]]] <- replicate(length(mis_partition[[i]]), sample(obs_partition[[i]],1))
        message("Specified number of donors is bigger than number of observed values in current partition that donors can be taken from.")
      }
      else
      {
        # grab y_hat, y_obs of current partition subset
        y_hat_obs_tmp <- y_hat_obs[obs_partition[[i]],]
        y_hat_mis_tmp <- y_hat_mis[mis_partition[[i]],]

        # now do nearest neighbor search, either by manhattan or euclidian norm
        if(distmetric == "manhattan")
        {
          nearest_neighbors <- RANN.L1::nn2(y_hat_obs_tmp, query = y_hat_mis_tmp, k = donors)
        }
        else
        {
          nearest_neighbors <- RANN::nn2(y_hat_obs_tmp, query = y_hat_mis_tmp, k = donors)
        }

        # sample from the donor pool of each missing observation
        # we want index of drawn match donor among all observed values
        # -> we draw from indices of current partition which is itself a set of indices, so we have to apply drawn index on current subset of partition
        idx[mis_partition[[i]]] <- obs_partition[[i]][apply(nearest_neighbors$nn.idx, 1, sample, size = 1)]
      }

    }

  }

  #return sample index, not value of y, as given y might contain dummy values
  return(y[ry,][idx,])
}




###########################################################################################################################################################
# match_univariate
#
# Univariate analogon to match_multivariate-function above, having the same functionalities, but, by nature of univariate matching, does not have parameters
# for metric or weight.
#
# Parameters:
# - y:  Vector that is to be imputed on.
# - y_hat: Vector containing all predictive means.
# - ry: Logical vector of length length(y) indicating the rows of observed data that are matched against.
# - wy: Logical vector of length length(y) indicating the rows in y for which imputations are created.
# - donors: Size of the donor pool that is computed and drawn from.
# - match_partition: List containing partitions of y[ry] and y[wy] that are based on the values of an additional column that has to be matched by which
#                    has been specified externally by the variable match_vars.
#
# Return value: Vector containing the imputed values for the missing data.
#
############################################################################################################################################################
match_univariate <- function(y, y_hat, ry, wy, donors, match_partitions)
{

  #if y is categoric, it has to be binary, so reduce to integer level -1 to get to binary standard
  ynum <- y
  if (is.factor(y))
    ynum <- as.integer(y) - 1

  #initialize y_obs and y_hat_obs/mis
  y_obs <- y[ry]
  y_hat_obs <- y_hat[ry]
  y_hat_mis <- y_hat[wy]

  #if number of donors is bigger than or equal to number of observed values, we can sample immediately
  if(donors >= sum(ry))
  {
    idx <- replicate(sum(wy), sample(1:sum(ry),1))
    return(idx)
  }


  # normalize y_hat matrices by median norm of shared rows
  med_norm <- median(abs(c(y_hat_mis, y_hat_obs)))
  y_hat_mis <- med_norm * y_hat_mis
  y_hat_obs <- med_norm * y_hat_obs

  # add random epsilon error on normalized y_hat_mis to simulate a random tie break in RANN
  eps_vec <- runif(length(y_hat_mis), min = 1e-15, max = 1e+15)
  y_hat_mis <- y_hat_mis + eps_vec;

  # check whether we also match by external data
  # -> if yes, match on partitioned data
  if(is.null(match_partitions))
  {
    # draw nearest neighbors
    nearest_neighbors <- RANN::nn2(y_hat_obs, query = y_hat_mis, k = donors)

    # sample from the donor pool of each missing observation
    idx <- apply(nearest_neighbors$nn.idx, 1, sample, size = 1)
  }
  else
  {
    idx <- 1:sum(wy)

    obs_partition <- match_partitions$obs_partition
    mis_partition <- match_partitions$mis_partition

    for(i in seq_along(obs_partition))
    {

      # if number of donors is bigger than or equal to number of observed values, we can sample immediately
      if(donors >= length(obs_partition[[i]]))
      {
        idx[mis_partition[[i]]] <- replicate(length(mis_partition[[i]]), sample(obs_partition[[i]],1))
      }
      else
      {
        # grab y_hat, y_obs of current partition subset
        y_hat_obs_tmp <- y_hat_obs[obs_partition[[i]]]
        y_hat_mis_tmp <- y_hat_mis[mis_partition[[i]]]

        # now do nearest neighbor search
        nearest_neighbors <- RANN::nn2(y_hat_obs_tmp, query = y_hat_mis_tmp, k = donors)

        # sample from the donor pool of each missing observation
        # we want index of drawn match donor among all observed values
        # -> we draw from indices of current partition which is itself a set of indices, so we have to apply drawn index on current subset of partition
        idx[mis_partition[[i]]] <- obs_partition[[i]][apply(nearest_neighbors$nn.idx, 1, sample, size = 1)]
      }

    }

  }

  # return sample index, not value of y, as given y might contain dummy values
  return(y_obs[idx])
}
