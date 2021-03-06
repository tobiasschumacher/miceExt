% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mice.post.matching.r
\name{mice.post.matching}
\alias{mice.post.matching}
\title{Post-processing of Imputed Data by Multivariate Predictive Mean Matching}
\usage{
mice.post.matching(obj, blocks = NULL, donors = 5L, weights = NULL,
  distmetric = "residual", matchtype = 1L, match_vars = NULL,
  ridge = 1e-05, minvar = 1e-04, maxcor = 0.99)
}
\arguments{
\item{obj}{\code{mice::mids} object that has been returned by a previous run
of \code{mice()} or \code{mice.mids()} and whose imputations we want to
post-process.}

\item{blocks}{Vector or list of vectors specifying the column tuples that are
to be imputed on blockwise. For a detailed explanation about the valid
formats of this parameter, see section \strong{\emph{input formats}} below.
Each element of each specified block has to be included in the
\code{visitSequence} of the given \code{mids} object, and the imputation
method that was applied on each of these columns has to be either \code{pmm}
or \code{norm}. \cr Univariate blocks are also allowed if they have been
imputed on via \code{mice.impute.norm()}, or if we want to match them
against an observed variable that is specified in the parameter
\code{match_vars}. \cr The default is \code{blocks = NULL}, which tells this
function to automatically determine and impute on column blocks that have
identical missing value patterns.}

\item{donors}{Integer or integer vector indicating the size of the donor pool
among which a draw in the matching step is made. If only a single number of
donors is specified, it is applied on all blocks, otherwise this parameter
needs to be a vector with as many elements as there blocks, specifying for
each of these blocks how many donors are to be drawn from. \cr
The default is \code{donors = 5L}. Setting \code{donors = 1L} always selects
the closest match (nearest neighbor), but is not recommended.}

\item{weights}{Numeric vector or list of numeric vectors that allocates
weights to the columns of each block, giving us the possibility to punish or
mitigate differences in certain columns of the data when computing the
distances in the matching step to determine the donor pool. Further details
on the valid formats of this parameter can be found below in section
\strong{\emph{input formats}}. \cr
The default is \code{weights = NULL}, meaning that no weights should be
applied to any column at all. \cr
Note that in order to avoid any ambiguities, specifying this parameter in
list format is only allowed if blocks have been explicitly specified and NOT
automatically determined.}

\item{distmetric}{Character string or character vector that determines which
mathematical metric we want to use to compute the distances between the
multivariate \code{y_obs} and \code{y_mis}. If only a single
metric is specified, it is applied on all blocks, otherwise this parameter
needs to be a vector with as many elements as there blocks, specifying for
each of these blocks which metric there is to use. \cr
The options are \code{"euclidian"}, \code{"manhattan"}, \code{"mahalanobis"}
and \code{"residual"}. The first three options refer to the distance metrics
of the same name, the latter is a variant of the mahalanobis distance in
which we consider the residual covariance \code{cov(y_hat_obs - y_obs)} of
the predicted model. This distance has been proposed and recommended by
Little (1988) and is also the default. Note that the first two options are
faster though, as they do not involve any matrix computations.}

\item{matchtype}{Integer of integer vector indicating the type of matching
distance to be used in PMM. If only a single matchtype is specified, it is
applied on all blocks, otherwise this parameter needs to be a vector with as
many elements as there blocks, specifying for each of these blocks which
matchtype has to be used. \cr
The default choice (\code{matchtype = 1L}) calculates the distance between
the \emph{predicted} values of \code{y_obs} and the \emph{drawn} values of
\code{y_mis} (called type-1 matching). Other choices are \code{matchtype =
0L} distance between predicted values) and \code{matchtype = 2L} (distance
between drawn values).}

\item{match_vars}{Vector specifying for each tuple which additional variable
has to be matched against. Can be an integer or character vector, either
specifying column indices or column names. \code{match_vars} must be the
same length as \code{blocks} with \code{0}-elements or \code{""}-elements to
disable this functionality for certain groups. Default is \code{match_vars =
NULL}, which completely disables this functionality.}

\item{ridge}{The ridge penalty used in an internal call of
\code{mice::.norm.draw()} to prevent problems with multicollinearity. Can be
a single number or a numeric vector. If only a single ridge value is
specified, it is applied on all blocks, otherwise this parameter needs to be
a vector with as many elements as there blocks, specifying for each of these
blocks which ridge value has to be used. \cr
The default is \code{ridge = 1e-05}, which means that 0.01 percent of the
diagonal is added to the cross-product. Larger ridges may result in more
biased estimates. For highly noisy data (e.g. many junk variables), set
\code{ridge = 1e-06} or even lower to reduce bias. For highly collinear
data, set \code{ridge = 1e-04} or higher.}

\item{minvar}{The minimum variance that we require predictors to have when
building a linear model to compute the \code{y_hat}-values. Can be a single
number or a numeric vector. If only a single value is specified, it is
applied on all blocks, otherwise this parameter needs to be a vector with as
many elements as there blocks, specifying for each of these blocks which
minimum variance will b allowed. The default is \code{minvar = 1e-04}.}

\item{maxcor}{The maximum correlation that we allow predictors to have when
building a linear model to compute the \code{y_hat}-values.  Can be a single
number or a numeric vector. If only a single value is specified, it is
applied on all blocks, otherwise this parameter needs to be a vector with as
many elements as there blocks, specifying for each of these blocks which
maximum variance will be allowed. The default is \code{maxcor = 0.99}.}
}
\value{
List containing the following two elements:
\describe{

 \item{\strong{\code{midsobj}}}{\code{mice::mids} object that differs from the input object
 only in the imputations that have been post-processed, and the \code{call}
 and \code{loggedEvents} attributes that have been updated. In particular,
 those post-processed imputations are not affecting the \code{chainMean} or
 \code{chainVar}-attributes, and hence, \code{plot()} will not consider them
 either.}
 
 \item{\strong{\code{blocks}}}{Set of column blocks in list format that multivariate
 imputation has been performed on. It is equivalent to the input parameter
 \code{blocks} if it has been specified by the user, otherwise those column
 tuples have been determined internally.}
}
}
\description{
Performs multivariate predictive mean matching (PMM) on a set of columns that
have been imputed on by the functionalities of the \code{mice} package.
Also offers a functionality to match imputations against observed variables.
}
\section{Input Formats}{
 
 Within \code{mice.post.matching()} and \code{mice.binarize()}, there are two
 formats that can be used to specify the input parameters \code{blocks} and
 \code{weights}, namely the \emph{list format} and the \emph{vector format}.
 The basic idea behind the list format is that we exclusively specify
 parameters for those column blocks that we want to impute on and summarize
 them in a list where each element is a vector that represents one such
 block, while in the vector format we use a single vector in which each
 element represents one column in the data, and therefore also specify
 information for columns that we do not want to impute on. The exact use of
 these two formats on both the \code{blocks} and the \code{weights} parameter
 will be illustrated in the following.
 
 \describe{
 
   \item{\strong{\code{blocks}}}{
   \strong{1. List Format} \cr
   To specify the imputation \code{blocks} using the list format, a list of
   atomic vectors has to be passed to the \code{blocks} parameter, and each
   vector in this list represents a column block that has to be imputed on.
   These vectors can either contain the names of the columns in this block as
   character strings or the corresponding column indices in numerical format.
   Note that this input list must not contain any duplicate columns among and
   within its elements. If there is only a single block to impute on, a
   single atomic vector representing this tuple may also be passed to
   \code{blocks}.
   
   \bold{Example:}\cr 
   Within this and the following examples of this section, we consider the
   \code{mammal_data} data set which contains 11 columns, out of which the
   column tuples \code{(sws, ps)} and \code{(mls, gt)} have identical missing
   value patterns (check \code{?mammal_data} for further details on the
   data).\cr
   If we now wanted to specify these blocks in list format, we would have to
   write \cr
   \code{blocks = list(c("sws","ps"), c("mls", "gt"))} \cr
   or analogously, when using column indices instead, \cr
    \code{blocks = list(c(4,5), c(7,8))}.

   \strong{2. Vector Format} \cr
   If we want to specify the imputation \code{blocks} via the vector format,
   a single vector with as many elements as there are columns in the data has
   to be used. Each element of this vector assigns a block number to its
   corresponding column, and columns that have the same number are imputed
   together, while columns that are not to be imputed have to carry the value
   \code{0}. All block numbers have to be integral, starting from \code{1},
   and the total number of imputation blocks has to be the maximum block
   number.
   
   \bold{Example:}\cr 
   Again we want to blockwise impute the column tuples \code{(sws, ps)} and
   \code{(mls, gt)} from \code{mammal_data}. To specify these blocks via the
   vector notation, we assign block number \code{1} to the columns of the
   first tuple and group number \code{2} to the columns of the second tuple,
   while all other columns have block number \code{0}. Hence, we would pass
   \cr
   \code{blocks = c(0,0,0,1,1,0,2,2,0,0,0)} \cr
   to \code{mice.post.matching()}.
   }

   \item{\strong{\code{weights}}}{
   \strong{1. List Format} \cr 
   To specify the imputation \code{weights} using the list format, the
   corresponding imputations blocks must have been specified in list format
   as well. In this case, the \code{weights} list has to be the same length
   as blocks, and each of its elements has a numeric vector of the same
   length as its corresponding column block, thereby assigning each of its
   columns a (strictly postitive) weight. If we do not want to apply weights
   to a single tuple, we can write a single \code{0}, \code{1} or \code{NULL}
   in the corresponding spot of the list.
   
   \bold{Example:}\cr 
   In our example, we want to assign the \code{sws} column a \code{1.5} times
   higher weight than \code{ps}, while not assigning any values to the second
   tuple at all. To achieve that, we specify \cr
   \code{weights = list(c(3,2), NULL)}.
   
   \strong{2. Vector Format} \cr
   When specifying the imputation \code{weights} via the vector format, once
   again a single vector with as many elements as there are columns in the
   data has to be used, in which each element assigns a weight to its
   corresponding column. Weights of columns that are not imputed on will have
   no effect, while in all blocks that weights should not be applied on, each
   column should carry the same value. In general, the value \code{1} should
   be used as the standard weight value for all columns that either are not
   imputed on or that weights are not applied on.
   
   \bold{Example:}\cr 
   We want to assign the same weights to the first tuple as in the previous
   example, while not assigning weights the second block again. In this case,
   we would use the vector \cr
   \code{weights = c(1,1,1,3,2,1,1,1,1,1,1)} \cr
   in which all the values of the unimputed and unweighted columns are
   \code{1}.
   }
 
 }
 
 Internally, \code{mice.post.matching()} converts both parameters into the
 list format as this is more natural to iterate on. Hence, the output
 \code{blocks} parameter also is in list format.
}

\section{Algorithmic Details}{

 The algorithm basically iterates over the \code{m} imputations of the
 input \code{mice::mids} object and each column tuple in \code{blocks}, each
 time following these two main steps:
 \describe{
   \item{\strong{Prediction}}{First, we iterate over the columns of the
   current block, collecting the complete column \eqn{y} on which we want to
   impute, along with the matrix \eqn{x} of its predictors. Among the values
   of eqn{y}, we identify the observed values \eqn{y_{obs}} and the missing
   values \code{y_{mis}} along with their corresponding predictors
   \eqn{x_{obs}} and \code{x_{mis}}, restricting ourselves to only those
   values whose predictors are complete, i.e. don't contain any missing
   values themselves. Note that if the sets of predictor variables strongly
   vary among all the columns in the current tuple, it may occur that there
   are no common predictors in which case the function breaks. \cr
   From the observed values and their predictors, we build a Bayesian linear
   regression model and, depending on the specified matching type, use the
   model's coefficients and drawn regression weights \eqn{\beta} to compute
   the predicted means \eqn{\hat{y}_{obs}} and \eqn{\hat{y}_{mis}} which we
   then save in a matrix \eqn{\hat{Y}}.}
   
   \item{\strong{Matching}}{After iterating over the columns in the current
   tuple and building the matrix \eqn{\hat{Y}}, we match the multivariate
   predictions of the missing values in \eqn{\hat{Y}} against the predictions
   of the observed values. More precisely, for each \eqn{\hat{y}_{mis}}, we
   perform a k-nearest-neighbor search among all values \eqn{\hat{y}_{obs}},
   where \eqn{k} is the number of donors specified by the user, and then
   randomly sample one element of those \eqn{k} nearest neighbors which is
   then used to impute the missing value tuple in \eqn{y}. The distance
   metric that is used to determine the nearest neighbors is specified by the
   user as well, and in case of Euclidian or Manhattan distance its
   computation is very straightforward. If the Mahalanobis distance or the
   residual distance (as proposed by Little) has been selected, we first
   compute the corresponding covariance matrix and its (pseudo) inverse via
   its eigen decomposition, and then use it to transform the values
   \eqn{\hat{y}_{mis}} and \eqn{\hat{y}_{obs}} that are fed into the nearest
   neighbor search. If weights have been specified as well, they are also
   applied before the kNN-search. \cr
   If an additional observed variable to match against has been specified via
   \code{match_vars}, the set of predictors and missing values are
   partitioned by the values of the external variable first, and then the
   matching is performed within each pair of corresponding subsets of that
   partition.}

   Note that the imputed values are only stored as a result and, other than
   in the \code{mice} algorithm, are not used to compute predictive means for
   any other missing value. We exclusively use the imputed values that are
   provided within the input \code{mids} object.
 }
}

\examples{


\dontrun{
#------------------------------------------------------------------------------
# Example on modified 'mammalsleep' data set from mice, that has identical
# missing data patterns on the column tuples ('ps','sws') and ('mls','gt')
#------------------------------------------------------------------------------

# run mice on data set 'mammal_data' and obtain a mids object to post-process
mids_mammal <- mice(mammal_data)


# run function, as blocks have not been specified, it will automatically detect
# the column tuples with identical missing data patterns and then impute on
# these
post_mammal <- mice.post.matching(mids_mammal)

# read out which column tuples have been imputed on
post_mammal$blocks

# look into imputations within resulting mice::mids object
post_mammal$midsobj$imp



#------------------------------------------------------------------------------
# Example on original 'mammalsleep' data set from mice, in which we
# want to post-process the imputations in column 'sws' by only imputing values
# from rows whose value in 'pi' matches the value of 'pi' in the row we impute
# on.
#------------------------------------------------------------------------------

# run mice on data set 'mammal_data' and obtain a mids object to post-process
mids_mammal <- mice(mammalsleep)


# run function, specify 'sws' as the column to impute on, and specify 'pi' as
# the observed variable to consider in the matching
post_mammal <- mice.post.matching(mids_mammal, blocks = "sws", match_vars = "pi")


# look into imputations within resulting mice::mids object
post_mammal$midsobj$imp



#------------------------------------------------------------------------------
# Examples illustrating the combined usage of blocks and weights, relating to
# the examples in the input format section above. As before, we want to impute  
# on the column tuples ('ps','sws') and ('mls','gt') from mammal_data, while  
# this time assigning weights to the first block, in which 'ps' gets a 1.5 times
# higher weight than 'sws'. The second tuple is not weighted. 
#------------------------------------------------------------------------------

# run mice() first
mids_mammal <- mice(mammal_data)

## Now there are five options to specify the blocks and weights:

# First option: specify blocks and weights in list format
post_mammal <- mice.post.matching(obj = mids_mammal,
                                 blocks = list(c("sws","ps"), c("mls","gt")),
                                 weights = list(c(3,2), NULL))

# or equivalently, using colums indices:                                  
post_mammal <- mice.post.matching(obj = mids_mammal,
                                 blocks = list(c(4,5), c(7,8)),
                                 weights = list(c(3,2), NULL))

# Second option: specify blocks and weights in vector format
post_mammal <- mice.post.matching(obj = mids_mammal,
                                 blocks = c(0,0,0,1,1,0,2,2,0,0,0),
                                 weights = c(1,1,1,3,2,1,1,1,1,1,1))

# Third option: specify blocks in list format and weights in vector format
post_mammal <- mice.post.matching(obj = mids_mammal,
                                 blocks = list(c("sws","ps"), c("mls","gt")),
                                 weights = c(1,1,1,3,2,1,1,1,1,1,1))

# Fourth option: specify blocks in vector format and weights in list format.
# Note that the block number determines which tuple in the weights list it
# corresponds to, and within each tuple in the list the weight correspondence is
# determinded by left to right order of the data columns
post_mammal <- mice.post.matching(obj = mids_mammal,
                                 blocks = c(0,0,0,1,1,0,2,2,0,0,0),
                                 weights = list(c(3,2), NULL))

# Fifth option: only specify weights in vector format. If the user knows
# beforehand that at least the column tuple he wants to impute and use weights
# on have the same missing value patterns, he can assign weights to these using
# the vector format, while letting mice.post.matching() find all other blocks
# with identical missing value patterns - possibly even more than just
# ('ps','sws') and ('mls','gt')
post_mammal <- mice.post.matching(obj = mids_mammal,
                                 weights = c(1,1,1,3,2,1,1,1,1,1,1))



#------------------------------------------------------------------------------
# Example that illustrates the combined functionalities of mice.binarize(),
# mice.factorize() and mice.post.matching() on the data set 'boys_data', which
# contains the column blocks ('hgt','bmi') and ('hc','gen','phb') that have
# identical missing value patterns, and out of which the columns 'gen' and
# 'phb' are factors. We are going to impute both tuples blockwise, while
# binarizing the factor columns first. Note that we never need to specify any
# blocks or columns to binarize, as these are all determined automatically 
#------------------------------------------------------------------------------

# By default, mice.binarize() expands all factor columns that contain NAs,
# so the columns 'gen' and 'phb' are automatically binarized
boys_bin <- mice.binarize(boys_data)

# Run mice on binarized data, note that we need to use boys_bin$data to grab
# the actual binarized data and that we use the output predictor matrix
# boys_bin$pred_matrix which is recommended for obtaining better imputation
# models
mids_boys <- mice(boys_bin$data, predictorMatrix = boys_bin$pred_matrix)

# It is very likely that mice imputed multiple ones among one set of dummy
# variables, so we need to post-process. As recommended, we also use the output
# weights from mice.binarize(), which yield a more balanced weighting on the
# column tuple ('hc','gen','phb') within the matching. As in previous examples,
# both tuples are automatically discovered and imputed on
post_boys <- mice.post.matching(mids_boys, weights = boys_bin$weights)

# Now we can safely retransform to the original data, with non-binarized
# imputations
res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)

# Analyze the distribution of imputed variables, e.g. of the column 'gen',
# using the mice version of with()
with(res_boys, table(gen))



#------------------------------------------------------------------------------
# Similar example to the previous, that also works on 'boys_data' and
# illustrates some more advanced funtionalities of all three functions in miceExt: 
# This time we only want to post-process the column block ('gen','phb'), while
# weighting the first of these tuples twice as much as the second. Within the
# matching, we want to avoid matrix computations by using the euclidian distance
# to determine the donor pool, and we want to draw from three donors only.
#------------------------------------------------------------------------------

# Binarize first, we specify blocks in list format with a single block, so we 
# can omit an enclosing list. Similarly, we also specify weights in list format.
# Both blocks and weights will be expanded and can be accessed from the output
# to use them in mice.post.matching() later
boys_bin <- mice.binarize(boys_data, 
                         blocks = c("gen", "phb"), 
                         weights = c(2,1))

# Run mice on binarized data, again use the output predictor matrix from
# mice.binarize()
mids_boys <- mice(boys_bin$data, predictorMatrix = boys_bin$pred_matrix)

# Post-process the binarized columns. We use blocks and weights from the previous
# output, and set 'distmetric' and 'donors' as announced in the example
# description
post_boys <- mice.post.matching(mids_boys,
                               blocks = boys_bin$blocks,
                               weights = boys_bin$weights,
                               distmetric = "euclidian",
                               donors = 3L)

# Finally, we can retransform to the original format
res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)
}


}
\references{
Little, R.J.A. (1988), Missing data adjustments in large surveys (with
 discussion), Journal of Business Economics and Statistics, 6, 287--301.

 Van Buuren, S. (2012). Flexible Imputation of Missing Data. CRC/Chapman \&
 Hall, Boca Raton, FL.

 Van Buuren, S., Groothuis-Oudshoorn, K. (2011). \code{mice}: Multivariate
 Imputation by Chained Equations in \code{R}. \emph{Journal of Statistical
 Software}, \bold{45}(3), 1-67. \url{http://www.jstatsoft.org/v45/i03/}
}
\seealso{
\code{\link[mice]{mice}}, \code{\link[mice]{mids-class}},
 \code{\link[miceExt]{mice.binarize}}, \code{\link[miceExt]{mice.factorize}}
}
\author{
Tobias Schumacher, Philipp Gaffert
}
