\name{cereals}
\alias{cereals}
\docType{data}
\title{Cereals datset}
\description{
  Data with several variables of different brands of cereal
}
\usage{data(cereals)}
\format{
  A data frame with 77 observations on the following 15 variables.
  \describe{
    \item{\code{mfr}}{Manufacturer of cereal}
    \item{\code{type}}{type: cold or hot}
    \item{\code{calories}}{calories per serving}
    \item{\code{protein}}{grams of protein}
    \item{\code{fat}}{grams of fat}
    \item{\code{sodium}}{milligrams of sodium}
    \item{\code{fiber}}{grams of dietary fiber}
    \item{\code{carbo}}{grams of complex carbohydrates}
    \item{\code{sugars}}{grams of sugars}
    \item{\code{potass}}{milligrams of potassium}
    \item{\code{vitamins}}{vitamins and minerals - 0, 25, or 100, indicating the typical percentage of FDA recommended}
    \item{\code{shelf}}{display shelf (1, 2, or 3, counting from the floor)}
    \item{\code{weight}}{weight in ounces of one serving}
    \item{\code{cups}}{number of cups in one serving}
    \item{\code{rating}}{a rating of the cereals}
  }
}
\source{
  \url{http://lib.stat.cmu.edu/DASL/Datafiles/Cereals.html}
}
\references{
 \url{http://lib.stat.cmu.edu/DASL/Stories/HealthyBreakfast.html}
}
\examples{
# load data
data(cereals)

# take a peek
head(cereals)
}
\keyword{datasets}
