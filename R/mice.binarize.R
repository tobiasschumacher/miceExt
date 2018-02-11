#'Binarize Factor Columns in Data Frames
#'
#'This function replaces factor columns in data frames in-place by a set of
#'binary columns which represent the so-called one-hot encoding of this factor.
#'More precisely, a column of a factor with \code{n} levels will be transformed
#'into a set of \code{n} binary columns, each representing exactly one category
#'of the original factor. Hence, the value \code{1} occurs in a column if and
#'only if the original factor had the value corresponding to that column. \cr
#'Further, this function also returns a predictor matrix that fits to the
#'binarized data frame and, when used as input parameter in \code{mice()},
#'ensures that binary columns that relate to the same factor do not predict each
#'other.
#'
#'@param data Matrix or data frame that contains factor columns which we want to
#'  convert into an equivalent set of binary columns.
#'@param include_ordered Logical variable indicating whether we also want to
#'  transform ordered factors. Default is \code{TRUE}.
#'@param include_observed Logical variable indicating whether we also want to
#'  transform factor columns in which all values are observed. Default is
#'  \code{FALSE}.
#'@param cols Numerical vector corresponding to the indices or character vector
#'  corresponding to the names of factor columns which we want to transform. By
#'  default, its value is \code{NULL}, indicating that the algorithm
#'  automatically identifies all factor columns that are to be binarized. If
#'  however the user specifies its value, the function exclusively transforms
#'  the specified columns, ignoring the values of the other optional parameters.
#'@param pred_matrix A custom predictor matrix relating to input \code{data},
#'  which will get transformed into the format that fits to the binarized output
#'  data frame. The result of this transformation will be stored in the
#'  \code{pred_matrix}-element of the output and should then be used as the
#'  \code{predictorMatrix} parameter in \code{mice()} to ensure that binary
#'  columns relating to the same factor column in the original data do not
#'  predict each other, yielding cleaner imputation models. If not specified,
#'  the default is the massive imputation predictor matrix.
#'@return List containing the following three elements:
#'\describe{
#'  \item{data}{ The binarized data frame.}
#'  \item{par_list}{A list containing the original data frame as well as some
#'  parameters with further information on the transformation. This list is
#'  needed to retransform the (possibly imputed) data at later stager via the
#'  \code{mice.factorize()} function, and should not be edited by the user
#'  under any circumstance. Next to the original data, the most notable
#'  element of this list would be "dummy_cols", which itself is a list of the
#'  column tuples that correspond to the transformed factor column from the
#'  original data set, and therefore works perfectly as input for
#'  \code{mice.post.matching()} (cf. examples below).}
#'  \item{pred_matrix}{Transformed version of input \code{pred_matrix} that
#'  should be used as the input argument \code{predictorMatrix} of
#'  \code{mice()}.}
#'}
#'@author Tobias Schumacher, Philipp Gaffert
#'
#'@examples
#'
#'#------------------------------------------------------------------------------
#'# first set of examples illustrating basic functionality
#'#------------------------------------------------------------------------------
#'
#'# binarize all factor columns in boys_data that contain NAs
#'boys_bin <- mice.binarize(boys_data)
#'
#'
#'# binarize only column 'gen' in boys_data
#'boys_bin <- mice.binarize(boys_data, cols = c("gen"))
#'
#'# read out binarized data
#'boys_bin$data
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
#'# now we can safely retransform to the original data, with non-binarized
#'# imputations
#'res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)
#'
#'# analyze the distribution of imputed variables, e.g. of the column 'gen',
#'# using the mice version of with()
#'with(res_boys, table(gen))
#'}
#'
#'
#'@seealso \code{\link[miceExt]{mice.factorize}},
#'  \code{\link[miceExt]{mice.post.matching}}, \code{\link[mice]{mice}}
#'@export
mice.binarize <- function(data, include_ordered = TRUE, include_observed = FALSE, cols = NULL, pred_matrix = (1 - diag(1, ncol(data))))
{
  ## check whether input data is valid

  #basic type checks
  if (!(is.matrix(data) || is.data.frame(data)))
    stop("Data should be a matrix or data frame")

  data <- as.data.frame(data) # cast to data frame

  # grab number of columns of source data, check whether we have more than two columns
  n_src_cols <- ncol(data)
  if (n_src_cols < 2)
    stop("Data should contain at least two columns")


  ## check optional parameters:

  # check predictor matrix
  check_pred_matrix(pred_matrix,  ncol(data))

  # if cols has been specified, check whether it is valid and then binarize these cols, ignoring other optionals -
  # otherwise scan whole data frame and binarize w.r.t. other optionals

  if(!is.null(cols))
  {
    src_factor_cols <- check_cols_binarize(cols, data)
    n_factors <- length(src_factor_cols)
    v_padded_column_counts <- rep(1, n_src_cols)
    v_padded_column_counts[src_factor_cols] <- unlist(lapply(src_factor_cols, function(j) length(levels(data[,j]))))
  }
  else
  {
    # get number of columns in binarized data of each attribute, non-categorical, non-binary attributes stay at one column
    v_padded_column_counts <- unlist(lapply(1:n_src_cols,
      function(j)
      {
        col <- data[,j]
        if(is.factor(col) && length(levels(col)) > 2 && (include_observed || !all(!is.na(col)) ))
        {
          if(!include_ordered)
          {
            if(!is.ordered(col))
              return(length(levels(col)))
            else
              return(1)
          }
          else
            return(length(levels(col)))
        }
        else
          return(1)
      }))

    # read out factor cols that are to be binarized
    src_factor_cols <- which(v_padded_column_counts > 1)

    # grab number of factors that have to be transformed
    n_factors = length(src_factor_cols)

    # check whether there actually is something to binarize and abort if not
    if(n_factors == 0)
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
  }


  # initialize result lists
  dummy_cols <- vector(mode = "list", length = n_factors)
  src_levels <- vector(mode = "list", length = n_factors)

  # grab number of columns that padded data has
  n_padmodel_cols <- sum(v_padded_column_counts)

  # initialize frame of padded data
  padded_data <- data.frame(matrix(NA, nrow = nrow(data), ncol = n_padmodel_cols))

  # initialize result predictor matrix
  res_matrix <- matrix(NA, nrow = n_padmodel_cols, ncol = n_padmodel_cols)

  # prepare iteration, initialize indices for padded data and result lists
  pad_index <- 1
  list_index <- 1

  ##################################################################################################
  # MAIN ITERATION:
  # iterate over columns of source data, if column isn't a factor, copy it into padded data frame,
  # else create matrix of dummy columns like in mice::padmodel and insert it into padded data frame
  ##################################################################################################

  for(j in 1:n_src_cols)
  {

    # check whether current column has to be binarized
    if(!(j %in% src_factor_cols))
    {
      ## current column in source data won't be binarized

      # copy current column from source data into padded data frame and adapt column name
      padded_data[,pad_index] <- data[,j]
      names(padded_data)[pad_index] <- names(data)[j]

      # update resulting predictive matrix
      res_matrix[,pad_index] <- inflate_vector(pred_matrix[,j], v_padded_column_counts)

      # increment index of padded data
      pad_index <- pad_index + 1
    }
    else
    {
      ## current column in source data has to be binarized
      # -> binarize source colun like in mice::padmodel

      # grab number of columns in padded model that correspond to current column in source data frame
      n_dummy_cols <- v_padded_column_counts[j]

      # grab current column in source data
      cat.column <- data[!is.na(data[, j]), j]

      # grab indices of corresponding columns in padded data frame
      curr_dummy_indices <- pad_index:(pad_index + n_dummy_cols - 1)

      # create binary matrix and copy it into padded data
      padded_data[!is.na(data[, j]), curr_dummy_indices] <- model.matrix(~ cat.column - 1)

      # name new columns
      names(padded_data)[curr_dummy_indices] <- paste(attr(data, "names")[j], levels(cat.column), sep = ".")

      # update documentation variables
      dummy_cols[[list_index]] <- curr_dummy_indices
      src_levels[[list_index]] <- levels(data[,j])

      # update resulting predictor matrix
      res_matrix[,curr_dummy_indices] <- inflate_vector(pred_matrix[,j], v_padded_column_counts)

      # increment indices
      pad_index <- pad_index + n_dummy_cols
      list_index <- list_index + 1
    }
  }

  pad_names <- names(padded_data)  # name rows and columns of predictor matrix
  rownames(res_matrix) <- pad_names
  colnames(res_matrix) <- pad_names

  # encapsulate all information of factors from source data that is needed to retransform via mice.factorize in one list
  par_list <- list(src_data = data, n_src_cols = n_src_cols, n_pad_cols = ncol(padded_data), src_factor_cols = src_factor_cols, dummy_cols = dummy_cols, src_names = names(data), pad_names = pad_names, src_levels = src_levels)

  # return
  return(list(data = padded_data, par_list = par_list, pred_matrix = res_matrix))
}

