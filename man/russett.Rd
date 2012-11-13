\name{russett}
\alias{russett}
\docType{data}
\title{Russett dataset}
\description{
  Data set from Russett (1964) about agricultural inequality, industrial development and political instability.
}
\usage{russett}
\format{
  A data frame with 47 observations on the following 11 variables. The variables may be used to construct three latent concepts: 1) \code{AGRIN}=Agricultural Inequality, 2) \code{INDEV}=Industrial Development, 3) \code{POLINS}=Political Instability. 

  \tabular{llll}{
    \emph{Num} \tab \emph{Variable}  \tab \emph{Description}			\tab \emph{Concept}	\cr
    1   \tab gini 		\tab Inequality of land distribution 			\tab \code{AGRIN} 	\cr
    2   \tab farm 		\tab Percentage of farmers that own half of the land 	\tab \code{AGRIN} 	\cr
    3   \tab rent 		\tab Percentage of farmers that rent all their land 	\tab \code{AGRIN} 		\cr
    4   \tab gnpr 	  \tab Gross national product per capita 	\tab \code{INDEV} 		\cr
    5   \tab labo 	  \tab Percentage of labor force employed in agriculture 		\tab \code{INDEV} 		\cr    
    6   \tab inst 		\tab Instability of executive (1945-1961) 	\tab \code{POLINS} 		\cr
    7   \tab ecks 	  \tab Number of violent internal war incidents (1941-1961)		\tab \code{POLINS} 		\cr
    8   \tab death 		\tab Number of people killed as a result of civic group violence (1950-1962)	\tab \code{POLINS} 		\cr    
    9   \tab demostab \tab Political regime: stable democracy	\tab \code{POLINS} 		\cr
    10  \tab demoinst \tab Political regime: unstable democracy		\tab \code{POLINS} 		\cr
    11  \tab dictator \tab Political regime: dictatorship	\tab \code{POLINS} 		\cr
  }
}
\references{
  Russett B.M. (1964) Inequality and Instability: The Relation of Land Tenure to Politics. \emph{World Politics} \bold{16:3}, pp. 442-454.
}
\examples{
  data(russett)
  russett
}
\keyword{datasets}
