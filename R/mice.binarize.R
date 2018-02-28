#'Binarize Factor Columns in Data Frames
#'
#'This function replaces factor columns in data frames in-place by a set of
#'binary columns which represent the so-called one-hot encoding of this factor.
#'More precisely, a column of a factor with \code{n} levels will be transformed
#'into a set of \code{n} binary columns, each representing exactly one category
#'of the original factor. Hence, the value \code{1} occurs in a column if and
#'only if the original factor had the value corresponding to that column. \cr
#'Further, this function also returns a weights vector and a predictor matrix
#'that fit to the binarized data frame. The weights vector is recommended to be
#'used as input for \code{mice.post.matching()}, as this avoids the effect of
#'overweighing a big set of binary columns relating to one factor against other
#'columns within one imputation block, while the predictor matrix, when used as
#'input parameter in \code{mice()}, ensures that binary columns that relate to
#'the same factor do not predict each other.
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
#'  the specified columns, ignoring the values of the previous optional parameters.
#'@param blocks Optional vector or list of vectors specifying the column tuples
#'  that are to be imputed on blockwise later when running
#'  \code{mice.post.matching()}. If this parameter is specified, the function
#'  will dectect all factor columns withing this argument and binarize these,
#'  while also producing a corresponding expanded block vector that can be found
#'  in the output of this function. In particular, whenever this parameter is
#'  specified, the values of the previous three parameters will be ignored. For
#'  a detailed explanation about the valid formats of this parameter, see
#'  section \strong{\emph{input formats}} below.\cr
#'  The default of this parameter is \code{blocks = NULL}, in which case the
#'  output parameter \code{blocks} will also be \code{NULL}.
#'@param weights Optional numeric vector or list of numeric vectors that
#'  allocates weights to the columns in the data, in particular on those that
#'  are going to be imputed on blockwise. If specified, the input will be
#'  transformed into vector format (if necessary), and columns that correspond
#'  to factor columns that are binarized will also be expanded. Within this
#'  expansion, all weights that belong to factor columns are thereby also
#'  divided by the number of levels of their corresponding factor column. This
#'  is done to ensure a balanced weighting in the later matching process within
#'  \code{mice.post.matching()}, as in the case that a factor column with many
#'  levels is in the same block with factor columns that much fewer levels or
#'  even with single numeric columns, the predictive means of the dummy columns
#'  of this factor would have much more impact on the matching than predictive
#'  means of the other column(s) simply because the latter would be
#'  outnumbered.\cr
#'  The default of this parameter is \code{weights = rep(1, ncol(data))}, which
#'  initially assigns the weight \code{1} to all columns before expanding them
#'  and reducing the weight of binarized columns as explained above. In any way,
#'  it is strongly recommended to use the transformed output \code{weights} as
#'  the \code{weights} parameter in \code{mice.post.matching()}. \cr
#'  Note that in order to avoid any ambiguities, specifying this parameter in
#'  list format is only allowed if blocks have been specified in list format as
#'  well.
#'@param pred_matrix A custom predictor matrix relating to input \code{data},
#'  which will get transformed into the format that fits to the binarized output
#'  data frame. The result of this transformation will be stored in the
#'  \code{pred_matrix} element of the output and should then be used as the
#'  \code{predictorMatrix} parameter in \code{mice()} to ensure that binary
#'  columns relating to the same factor column in the original data do not
#'  predict each other, yielding cleaner imputation models. If not specified,
#'  the default is the massive imputation predictor matrix.
#'  
#'  
#'@return List containing the following five elements:
#'\describe{
#'  \item{\strong{\code{data}}}{ The binarized data frame.}
#'  
#'  \item{\strong{\code{par_list}}}{A list containing the original data frame as
#'  well as some parameters with further information on the transformation. This
#'  list is needed to retransform the (possibly imputed) data at later stager
#'  via the \code{mice.factorize()} function, and should not be edited by the
#'  user under any circumstance. Next to the original data, the most notable
#'  element of this list would be "dummy_cols", which itself is a list of the
#'  column tuples that correspond to the transformed factor column from the
#'  original data set, and therefore works perfectly as input for
#'  \code{mice.post.matching()} (cf. examples below).}
#'  
#'  \item{\strong{\code{blocks}}}{If input parameter \code{blocks} has been
#'  specified, an expanded version of that input is returned in vector format
#'  via this element. It should be used as input parameter \code{blocks} in
#'  \code{mice.post.matching()} later, after imputing the binarized data via
#'  \code{mice()} first.}
#'  
#'  \item{\strong{\code{weights}}}{Transformed version of input parameter
#'  \code{weights} in vector format that should be used as input parameter
#'  \code{weights} in \code{mice.post.matching()} later, after imputing the
#'  binarized data via \code{mice()} first. If input parameter \code{weights}
#'  has not been specified, a default vector is still going to be output.}
#'  
#'  \item{\strong{\code{pred_matrix}}}{Transformed version of input
#'  \code{pred_matrix} that should be used as the input argument
#'  \code{predictorMatrix} of \code{mice()}.}
#'}
#'
#'
#'@section Input Formats: 
#'  Within \code{mice.binarize()} and \code{mice.post.matching()}, there are two
#'  formats that can be used to specify the input parameters \code{blocks} and
#'  \code{weights}, namely the \emph{list format} and the \emph{vector format}.
#'  The basic idea behind the list format is that we exclusively specify
#'  parameters for those column blocks that we want to impute on and summarize
#'  them in a list where each element is a vector that represents one such
#'  block, while in the vector format we use a single vector in which each
#'  element represents one column in the data, and therefore also specify
#'  information for columns that we do not want to impute on. The exact use of
#'  these two formats on both the \code{blocks} and the \code{weights} parameter
#'  will be illustrated in the following.
#'  
#'  \describe{
#'  
#'     \item{\strong{\code{blocks}}}{ 
#'     \strong{1. List Format} \cr
#'     To specify the imputation \code{blocks} using the list format, a list of
#'     atomic vectors has to be passed to the \code{blocks} parameter, and each
#'     vector in this list represents a column block that has to be imputed on.
#'     These vectors can either contain the names of the columns in this block
#'     as character strings or the corresponding column indices in numerical
#'     format. Note that this input list must not contain any duplicate columns
#'     among and within its elements. If there is only a single block to impute
#'     on, a single atomic vector representing this tuple may also be passed to
#'     \code{blocks}.
#'     
#'     \strong{Example:}\cr 
#'     Within this and the following examples of this section, we consider the
#'     \code{boys_data} data set which contains 9 columns, out of which the
#'     column tuples \code{(hgt, bmi)} and \code{(hc, gen, phb)} have identical
#'     missing value patterns, while the columns \code{gen} and \code{phb} are
#'     also categorical. (check \code{?boys_data} for further details on the
#'     data).\cr
#'     If we now wanted to specify these blocks in list format, we would have to
#'     write \cr
#'     \code{blocks = list(c("hgt","bmi"), c("hc", "gen", "phb"))} \cr
#'     or analogously, when using column indices instead, \cr
#'     \code{blocks = list(c(2,4), c(5,6,7))}.
#'     
#'     \strong{2. Vector Format} \cr
#'     If we want to specify the imputation \code{blocks} via the vector format,
#'     a single vector with as many elements as there are columns in the data
#'     has to be used. Each element of this vector assigns a block number to its
#'     corresponding column, and columns that have the same number are imputed
#'     together, while columns that are not to be imputed have to carry the
#'     value \code{0}. All block numbers have to be integral, starting from
#'     \code{1}, and the total number of imputation blocks has to be the maximum
#'     block number.
#'     
#'     \strong{Example:}\cr
#'     Again we want to blockwise impute the column tuples \code{(hgt, bmi)} and
#'     \code{(hc, gen, phb)} from \code{boys_data}. To specify these blocks via
#'     the vector notation, we assign block number \code{1} to the columns of
#'     the first tuple and group number \code{2} to the columns of the second
#'     tuple, while all other columns have block number \code{0}. Hence, we
#'     would pass \cr
#'     \code{blocks = c(0,1,0,1,2,2,2,0,0)} \cr
#'     to \code{mice.post.matching()}.
#'     }
#'     
#'     \item{\strong{\code{weights}}}{ 
#'      \strong{1. List Format} \cr
#'     To specify the imputation \code{weights} using the list format, the
#'     corresponding imputations blocks must have been specified in list format
#'     as well. In this case, the \code{weights} list has to be the same length
#'     as blocks, and each of its elements has a numeric vector of the same
#'     length as its corresponding column block, thereby assigning each of its
#'     columns a (strictly postitive) weight. If we do not want to apply weights
#'     to a single tuple, we can write a single \code{0}, \code{1} or
#'     \code{NULL} in the corresponding spot of the list.
#'     
#'     \strong{Example:}\cr 
#'     In our example, we want to assign the \code{hgt} column a \code{1.5}
#'     times higher weight than \code{bmi}, while not assigning any values to
#'     the second tuple at all. To achieve that, we specify \cr
#'     \code{weights = list(c(3,2), NULL)}.
#'     
#'     \strong{2. Vector Format} \cr
#'     When specifying the imputation \code{weights} via the vector format, once
#'     again a single vector with as many elements as there are columns in the
#'     data has to be used, in which each element assigns a weight to its
#'     corresponding column. Weights of columns that are not imputed on will
#'     have no effect, while in all blocks that weights should not be applied
#'     on, each column should carry the same value. In general, the value
#'     \code{1} should be used as the standard weight value for all columns that
#'     either are not imputed on or that weights are not applied on.
#'     
#'     \strong{Example:}\cr 
#'     We want to assign the same weights to the first tuple as in the previous
#'     example, while not assigning weights the second block again. In this
#'     case, we would use the vector \cr
#'     \code{weights = c(1,3,1,2,1,1,1,1,1)} \cr
#'     in which all the values of the unimputed and unweighted columns are
#'     \code{1}.
#'     }
#'  
#'  }
#'  
#'  Internally, \code{mice.post.matching()} converts both parameters into the
#'  vector format as this works best within the main iteration over all columns
#'  of the data. Hence, the output parameters \code{blocks} and \code{weights}
#'  are also in list format.
#'
#'
#'@author Tobias Schumacher, Philipp Gaffert
#'
#'@examples
#'
#'\dontrun{
#'#------------------------------------------------------------------------------
#'# first set of examples illustrating basic functionality
#'#------------------------------------------------------------------------------
#'
#'# binarize all factor columns in boys_data that contain NAs
#'boys_bin <- mice.binarize(boys_data)
#'
#'# binarize only column 'gen' in boys_data
#'boys_bin <- mice.binarize(boys_data, cols = c("gen"))
#'
#'# binarize all factor columns with the blocks ('hgt','bmi') and ('gen','phb')
#'# to impute on these (binarized) blocks later
#'boys_bin <- mice.binarize(boys_data, blocks = list(c("hgt","bmi"), c("gen", "phb")))
#'
#'# read out binarized data
#'boys_bin$data
#'
#'
#'
#'
#'#------------------------------------------------------------------------------
#'# Examples illustrating the combined usage of blocks and weights, relating to
#'# the examples in the input format section above. As before, we want to impute  
#'# on the column tuples ('hgt','bmi') and ('hc','gen','phb) from boys_data, while  
#'# assigning weights to the first block, in which 'hgt' gets a 1.5 times
#'# higher weight than 'bmi'. The second tuple is not weighted. 
#'#------------------------------------------------------------------------------
#'
#'## Now there are three options to specify the blocks and weights:
#'
#'# First option: specify blocks and weights in list format
#'boys_bin <- mice.binarize(data = boys_data,
#'                                 blocks = list(c("hgt","bmi"), c("hc","gen","phb")),
#'                                 weights = list(c(3,2), NULL))
#'                               
#'# or equivalently, using colums indices:
#'boys_bin <- mice.binarize(data = boys_data,
#'                            blocks = list(c(2,4), c(5,6,7)),
#'                            weights = list(c(3,2), NULL))
#'                            
#'# Second option: specify blocks in list and weights in vector format
#'post_mammal <- mice.binarize(data = boys_data,
#'                                    blocks = c(0,1,0,1,2,2,2,0,0),
#'                                    weights = c(1,3,1,2,1,1,1,1,1))
#'
#'# Third option: specify blocks in list format and weights in vector format
#'post_mammal <- mice.binarize(data = boys_data,
#'                                    blocks = list(c("hgt","bmi"), c("hc","gen", "phb")),
#'                                    weights = c(1,3,1,2,1,1,1,1,1))
#'
#'# check expanded blocks vector
#'boys_bin$blocks
#'
#'# check expanded weights vector
#'boys_bin$weights
#'
#'
#'
#'#------------------------------------------------------------------------------
#'# Example that illustrates the combined functionalities of mice.binarize(),
#'# mice.factorize() and mice.post.matching() on the data set 'boys_data', which
#'# contains the column blocks ('hgt','bmi') and ('hc','gen','phb') that have
#'# identical missing value patterns, and out of which the columns 'gen' and
#'# 'phb' are factors. We are going to impute both tuples blockwise, while
#'# binarizing the factor columns first. Note that we never need to specify any
#'# blocks or columns to binarize, as these are all determined automatically 
#'#------------------------------------------------------------------------------
#'
#'# By default, mice.binarize() expands all factor columns that contain NAs,
#'# so the columns 'gen' and 'phb' are automatically binarized
#'boys_bin <- mice.binarize(boys_data)
#'
#'# Run mice on binarized data, note that we need to use boys_bin$data to grab
#'# the actual binarized data and that we use the output predictor matrix
#'# boys_bin$pred_matrix which is recommended for obtaining better imputation
#'# models
#'mids_boys <- mice(boys_bin$data, predictorMatrix = boys_bin$pred_matrix)
#'
#'# It is very likely that mice imputed multiple ones among one set of dummy
#'# variables, so we need to post-process. As recommended, we also use the output
#'# weights from mice.binarize(), which yield a more balanced weighting on the
#'# column tuple ('hc','gen','phb') within the matching. As in previous examples,
#'# both tuples are automatically discovered and imputed on
#'post_boys <- mice.post.matching(mids_boys, weights = boys_bin$weights)
#'
#'# Now we can safely retransform to the original data, with non-binarized
#'# imputations
#'res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)
#'
#'# Analyze the distribution of imputed variables, e.g. of the column 'gen',
#'# using the mice version of with()
#'with(res_boys, table(gen))
#'
#'
#'
#'#------------------------------------------------------------------------------
#'# Similar example to the previous, that also works on 'boys_data' and
#'# illustrates some more advanced funtionalities of all three functions in miceExt: 
#'# This time we only want to post-process the column block ('gen','phb'), while
#'# weighting the first of these tuples twice as much as the second. Within the
#'# matching, we want to avoid matrix computations by using the euclidian distance
#'# to determine the donor pool, and we want to draw from three donors only.
#'#------------------------------------------------------------------------------
#'
#'# Binarize first, we specify blocks in list format with a single block, so we 
#'# can omit an enclosing list. Similarly, we also specify weights in list format.
#'# Both blocks and weights will be expanded and can be accessed from the output
#'# to use them in mice.post.matching() later
#'boys_bin <- mice.binarize(boys_data, 
#'                          blocks = c("gen", "phb"), 
#'                          weights = c(2,1))
#'
#'# Run mice on binarized data, again use the output predictor matrix from
#'# mice.binarize()
#'mids_boys <- mice(boys_bin$data, predictorMatrix = boys_bin$pred_matrix)
#'
#'# Post-process the binarized columns. We use blocks and weights from the previous
#'# output, and set 'distmetric' and 'donors' as announced in the example
#'# description
#'post_boys <- mice.post.matching(mids_boys,
#'                                blocks = boys_bin$blocks,
#'                                weights = boys_bin$weights,
#'                                distmetric = "euclidian",
#'                                donors = 3L)
#'
#'# Finally, we can retransform to the original format
#'res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)
#'}
#'
#'
#'@seealso \code{\link[miceExt]{mice.factorize}},
#'  \code{\link[miceExt]{mice.post.matching}}, \code{\link[mice]{mice}}
#'@export
mice.binarize <- function(data, include_ordered = TRUE, include_observed = FALSE, cols = NULL, blocks = NULL, weights = rep(1, ncol(data)), pred_matrix = (1 - diag(1, ncol(data))))
{
  ## check whether input data is valid

  # basic type checks
  if (!(is.matrix(data) || is.data.frame(data)))
    stop("Data should be a matrix or data frame")

  data <- as.data.frame(data) # cast to data frame

  # grab number of columns of source data, check whether we have more than two columns
  n_src_cols <- ncol(data)
  if (n_src_cols < 2)
    stop("Data should contain at least two columns")


  ## check optional parameters:
  
  # check whether input cols, blocks, weights and data work together, convert them into vector format if necessary
  # and return set of factor columns to convert
  setup <- check_blocks_weights_binarize(cols, blocks, weights, include_ordered, include_observed, data)
  
  # check predictor matrix
  check_pred_matrix(pred_matrix,  ncol(data))


  ## initialize main variables
  src_factor_cols <- setup$src_factor_cols
  n_factors = length(src_factor_cols)
  v_padded_column_counts <- rep(1, n_src_cols)
  v_padded_column_counts[src_factor_cols] <- unlist(lapply(src_factor_cols, function(j) nlevels(data[,j])))
  
  # initialize result lists
  dummy_cols <- vector(mode = "list", length = n_factors)
  src_levels <- vector(mode = "list", length = n_factors)

  # grab number of columns that padded data has
  n_padmodel_cols <- sum(v_padded_column_counts)

  # initialize frame of padded data
  padded_data <- data.frame(matrix(NA, nrow = nrow(data), ncol = n_padmodel_cols))

  # initialize result predictor matrix
  res_matrix <- matrix(NA, nrow = n_padmodel_cols, ncol = n_padmodel_cols)
  
  # initialize base and result column weights vector
  base_tuples <- setup$blocks
  if(is.null(base_tuples))
    res_blocks <- NULL
  else
    res_blocks <- rep(0, n_padmodel_cols)
  
  # initialize base and result column weights vector
  base_weights <- setup$weights
  res_weights <- rep(1, n_padmodel_cols)

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
      
      # update column tuples vector
      res_blocks[pad_index] <- base_tuples[j]
      
      # update result weights vector
      res_weights[pad_index] <- base_weights[j]

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
      
      # update column tuples vector
      res_blocks[curr_dummy_indices] <- base_tuples[j]
      
      # update column weights vector
      res_weights[curr_dummy_indices] <- base_weights[j]/n_dummy_cols

      # increment indices
      pad_index <- pad_index + n_dummy_cols
      list_index <- list_index + 1
    }
  }

  # grab names of binarized data
  pad_names <- names(padded_data)
  
  # name rows and columns of predictor matrix
  rownames(res_matrix) <- pad_names
  colnames(res_matrix) <- pad_names
  
  # if not NULL, name values of blocks vector
  if(!is.null(res_blocks))
    names(res_blocks) <- pad_names
  
  # name values of res_weights vector
  names(res_weights) <- pad_names

  # encapsulate all information of factors from source data that is needed to retransform via mice.factorize() in one list
  par_list <- list(src_data = data, n_src_cols = n_src_cols, n_pad_cols = ncol(padded_data), src_factor_cols = src_factor_cols, dummy_cols = dummy_cols, src_names = names(data), pad_names = pad_names, src_levels = src_levels)

  # return
  return(list(data = padded_data, par_list = par_list, blocks = res_blocks, weights = res_weights, pred_matrix = res_matrix))
}

