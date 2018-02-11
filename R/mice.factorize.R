#'Transform Imputations of Binarized Data Into Their Corresponding Factors
#'
#'This function acts as the counterpart to \code{mice.binarize}, as it
#'effectively retransforms imputations of binarized data that \code{mice} has
#'been run on and that has been post-processed via \code{mice.post.matching}
#'after. The post-processing is usually necessary as \code{mice} is very likely
#'to impute multiple ones among the dummy columns belonging to to a single
#'factor entry. The resulting \code{mice::mids} object is not suited for further
#'\code{mice.mids()} iterations or the use of \code{plot}, but works well as
#'input to \code{with()}.
#'
#'@param obj \code{mice::mids} object resulting from a call of
#'  \code{mice.post.matching()} and whose underlying data frame results from a
#'  call of \code{mice::binarize()}.
#'@param par_list List that has been returned in a previous call of
#'  \code{mice::binarize()} next to the underlying data of the argument
#'  \code{obj}.
#'@return A \code{mice::mids} object in which data and imputations have been
#'  retransformed from their respective binarized versions in the input
#'  \code{obj}. As this isn't a proper result of a mice iteration and many of
#'  the attributes of \code{obj} cannot be transformed well, only the slots
#'  \code{data}, \code{nmis}, \code{where} and \code{imp}, which are needed in
#'  \code{with()} are not \code{NULL}. Hence, it does not work as input for
#'  \code{mice.mids()}.
#'
#'@author Tobias Schumacher, Philipp Gaffert
#'
#'@examples
#'
#'
#'\dontrun{
#'#------------------------------------------------------------------------------
#'# this example illustrates the combined functionalities of mice.binarize,
#'# mice.factorize and mice.post.matching on the dataset 'boys' from mice, which
#'# yields different imputations on the factor columns 'gen', 'phb' and 'reg'
#'# than mice() would output
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
#'# now we can safely retransform to the original data, with non-binarized imputations
#'res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)
#'
#'# analyze the distribution of imputed variables, e.g. of the column 'gen',
#'# using the mice version of with()
#'with(res_boys, table(gen))
#'}
#'
#'
#'
#'@seealso \code{\link[miceExt]{mice.binarize}},
#'  \code{\link[miceExt]{mice.post.matching}}, \code{\link[mice]{mice}}
#'@export
mice.factorize <- function(obj, par_list)
{
  ## check whether input is valid
  if (!is.mids(obj))
    stop("Object should be of type mids.")

  check_par_list(obj, par_list)

  # grab function call
  call <- match.call()

  ## make local copies of some frequently used parameters
  src_factor_cols <- par_list$src_factor_cols
  dummy_cols <- par_list$dummy_cols
  n_src_cols <- par_list$n_src_cols

  # get number of columns in padded data frame
  n_padded_cols <- ncol(obj$data)

  # initialize result imputation list, nmis and where
  res_imp <- vector("list", length = n_src_cols)
  res_method <- vector(mode = "character", length = n_src_cols)
  res_nmis <- vector(mode = "numeric", length = n_src_cols)
  res_where <- matrix(FALSE, nrow = nrow(obj$data), ncol = n_src_cols)


  ## first copy imputations of non-categorical variables

  # grab numeric columns in source and padded data
  src_numeric_cols <- setdiff(1:n_src_cols,src_factor_cols)
  pad_data_numeric_cols <- setdiff(1:n_padded_cols, unlist(dummy_cols))

  # copy from mids$imp
  res_imp[src_numeric_cols] <- obj$imp[pad_data_numeric_cols]
  res_method[src_numeric_cols] <- obj$method[pad_data_numeric_cols]
  res_nmis[src_numeric_cols] <- obj$nmis[pad_data_numeric_cols]
  res_where[,src_numeric_cols] <- obj$where[,pad_data_numeric_cols]

  ###
  # MAIN STEP: transform binarized imputations back to categorical, use three convoluted iterations, where the inner two are carried out
  #    with apply/lapply functionalities
  # -> outermost iteration runs on target elements of result imputation list, gathering imputations for each eement
  # -> mid-level iteration runs over number of imputations which form the column of current entry in result imputation list
  # -> innermost iteration runs over all impuations of dummy atrributes of current target attribute, building a binary matrix representing
  #    the encoding of current target attribute, which is then transformed via another call of apply that scans each row of the binary
  #    matrix for the entry that is equal to 1
  ###

  #res_imp[src_factor_cols] <- lapply(seq_along(src_factor_cols),
  for(col_index in seq_along(src_factor_cols))
  {
    ## outer iteration: iterate over index of target column tuple

    # grab tmp values of current target column, corrensponding tuple of dummy columns, and list of levels of target column
    curr_src_index <- src_factor_cols[[col_index]]
    curr_dummy_tuple <- dummy_cols[[col_index]]
    curr_levels <- par_list$src_levels[[col_index]]

    # grab imputation of first dummy variable in current dummy tuple as reference
    # for naming and null-checking purposes
    ref_imp <- obj$imp[[curr_dummy_tuple[1]]]

    # if reference imputation is null, there are no imputed values -> skip to next iteration
    if(is.null(ref_imp))
      next

    # get imputations of current attribute
    curr_imp <- data.frame(lapply(1:obj$m,
      function(m)
      {
        ## mid-level iteration: iterate over imputation index m

        ## build binary matrix representing the encoding of current target attribute of current imputation
        # -> use inner iteration over current dummy columns
        bin_matrix <- do.call(cbind, lapply(obj$imp[curr_dummy_tuple], function(imp) imp[,m]))

        # return transformed factor column
        # -> transformation is obtained by scanning each row in binary matrix for the entry that is "1" via the apply function
        factor(apply(bin_matrix, MARGIN = 1,
          function(row)
          {
            # check whether binary encoding is correct, as mice() might return rows with multiple '1's and therefore requires use of mice.post.matching
            if(!all(row %in% c(0,1)) || sum(row) != 1)
              stop("The imputed values in the binarized columns are not in proper format. Maybe you forgot to run mice.post.matching() on input mids object.\n")

            # transform, column index that equals '1' indicates level to use
            return(curr_levels[which(row == 1)])
          }),
          levels = curr_levels)
      }))

    # copy row and column names from reference imputation
    colnames(curr_imp) <- colnames(ref_imp)
    rownames(curr_imp) <- rownames(ref_imp)

    res_imp[[curr_src_index]] <- curr_imp
    res_method[curr_src_index] <- "pmm"
    res_where[,curr_src_index] <- obj$where[,curr_dummy_tuple[1]]
    res_nmis[curr_src_index] <- obj$nmis[curr_dummy_tuple[1]]
  }

  # set names of elements of resulting imputation list
  names(res_imp) <- par_list$src_names

  # now build result mids object that can be used within mids
  ## save, and return
  midsobj <- list(call = call,
                  data = par_list$src_data,
                  where = res_where,
                  m = obj$m,
                  nmis = res_nmis,
                  imp = res_imp,
                  method = res_method,
                  predictorMatrix = NULL,
                  visitSequence = NULL,
                  post = NULL,
                  seed = obj$seed,
                  iteration = obj$iteration,
                  lastSeedValue = obj$lastSeedValue,
                  chainMean = NULL,
                  chainVar = NULL,
                  loggedEvents = obj$loggedEvents)

  oldClass(midsobj) <- "mids"

  # return result
  return(midsobj)
}
