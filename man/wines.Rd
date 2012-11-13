\name{wines}
\alias{wines}
\docType{data}
\title{Wines dataset}
\description{
  These data are the results of a chemical analysis of wines grown in the same region in Italy but derived from three different cultivars. The analysis determined the quantities of 13 constituents found in each of the three types of wines. 
}
\usage{wines}
\format{
  A data frame with 178 observations and 14 variables. 

  \tabular{lll}{
    \emph{Num} \tab \emph{Variable}	\tab \emph{Description}		\cr
    1   \tab class 		\tab Type of wine  			\cr
    2   \tab alcohol 		\tab Alcohol  				\cr
    3   \tab malic.acid 	\tab Malic acid 			\cr
    4   \tab ash	 	\tab Ash				\cr
    5   \tab alcalinity 	\tab Alcalinity 			\cr
    6   \tab magnesium		\tab Magnesium  			\cr
    7   \tab phenols		\tab Total phenols 			\cr
    8   \tab flavanoids		\tab Flavanoids				\cr
    9   \tab nofla.phen		\tab Nonflavanoid phenols		\cr
    10  \tab proantho		\tab Proanthocyanins			\cr
    11  \tab col.intens		\tab Color intensity			\cr
    12  \tab hue		\tab Hue				\cr
    13  \tab diluted		\tab OD280/OD315 of diluted wines	\cr
    14  \tab proline		\tab Proline				\cr
  }
}
\source{
  Machine Learning Repository. 
  \url{http://archive.ics.uci.edu/ml/datasets/Wine}
}
\references{
  Forina, M. et al, PARVUS \emph{An Extendible Package for Data Exploration, Classification and Correlation.} Institute of Pharmaceutical and Food Analysis  and Technologies, Via Brigata Salerno, 16147 Genoa, Italy. 
}
\examples{
  data(wines)
  wines
}
\keyword{datasets}
