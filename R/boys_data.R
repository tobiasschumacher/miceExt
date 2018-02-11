#' Growth of Dutch boys
#'
#' Modified sample of the \code{boys} data set that is included in \code{mice}.
#'
#'The original data set \code{boys} that is included in \code{mice} is a random
#'sample of 10\% from the cross-sectional data used to construct the Dutch
#'growth references 1997, and contains information height, weight, head
#'circumference and puberty of 748 Dutch boys. From this data, a sample of 100
#'rows has been taken and modified such that in each row, the values in the
#'columns \code{hc}, \code{gen} and \code{phb} are either blockwise \code{NA} or
#'blockwise non-\code{NA}.
#'
#'
#'@name boys_data
#'@format A dataset with 100 observations of the following 9 variables:
#'\describe{
#' \item{age}{Decimal age (0-21 years)}
#' \item{hgt}{Height (cm)}
#' \item{wgt}{Weight (kg)}
#' \item{bmi}{Body mass index}
#' \item{hc}{Head circumference (cm)}
#' \item{gen}{Genital Tanner stage (G1-G5)}
#' \item{phb}{Pubic hair (Tanner P1-P6)}
#' \item{tv}{Testicular volume (ml)}
#' \item{reg}{Region (north, east, west, south, city)}
#'}
#'@source Fredriks, A.M,, van Buuren, S., Burgmeijer, R.J., Meulmeester JF,
#'Beuker, R.J., Brugman, E., Roede, M.J., Verloove-Vanhorick, S.P., Wit, J.M.
#'(2000) Continuing positive secular growth change in The Netherlands
#'1955-1997.  \emph{Pediatric Research}, \bold{47}, 316-323.
#'
#'Fredriks, A.M., van Buuren, S., Wit, J.M., Verloove-Vanhorick, S.P. (2000).
#'Body index measurements in 1996-7 compared with 1980.  \emph{Archives of
#'Disease in Childhood}, \bold{82}, 107-112.-734.
#''@keywords datasets
#'@seealso \code{\link[mice]{boys}}
NULL
