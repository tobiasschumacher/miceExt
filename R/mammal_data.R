#'Mammal sleep data
#'
#'Modified version of the mammal sleep data set that is included in \code{mice}.
#'
#'The original dataset was from a study by Allison and Cicchetti (1976) of 62
#'mammal species on the interrelationship between sleep, ecological, and
#'constitutional variables, and was adapted in the \code{mice}-package. It
#'contains missing values on five variables and has been modified such that for
#'each row, the entries in the column tuples (\code{sws},\code{ps}) and
#'(\code{mls},\code{gt}) are either pairwise\code{NA} or pairwise non-\code{NA}.
#'
#'@name mammal_data
#'@format A dataset with 62 observations of the following 11 variables:
#'\describe{
#'  \item{species}{Species of animal}
#'  \item{bw}{Body weight (kg)}
#'  \item{brw}{Brain weight (g)}
#'  \item{sws}{Slow wave ("nondreaming") sleep (hrs/day)}
#'  \item{ps}{Paradoxical ("dreaming") sleep (hrs/day)}
#'  \item{ts}{Total sleep (hrs/day) (sum of slow wave and paradoxical sleep)}
#'  \item{mls}{Maximum life span (years)}
#'  \item{gt}{Gestation time (days)}
#'  \item{pi}{Predation index (1-5), 1 = least likely to be preyed upon}
#'  \item{sei}{Sleep exposure index (1-5), 1 = least exposed (e.g. animal sleeps
#'  in a well-protected den), 5 = most exposed}
#'  \item{odi}{Overall danger index (1-5) based on the above two indices and
#'  other information, 1 = least danger (from other animals),
#'  5 = most danger (from other animals)}
#'}
#'@source Allison, T., Cicchetti, D.V. (1976). Sleep in Mammals: Ecological and
#'Constitutional Correlates. Science, 194(4266), 732-734.
#'@keywords datasets
#'@seealso \code{\link[mice]{mammalsleep}}
NULL
