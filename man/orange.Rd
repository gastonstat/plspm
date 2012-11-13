\name{orange}
\alias{orange}
\docType{data}
\title{Orange Juice dataset}
\description{
  This data set contains the physico-chemical, sensory and hedonic measurements of 6 orange juices.
}
\usage{orange}
\format{
  A data frame with 6 observations and 112 variables. The variables refer to three latent concepts: 1) \code{PHYCHEM}=Physico-Chemical, 2) \code{SENSORY}=Sensory, and 3) \code{HEDONIC}=Hedonic.

  \tabular{llll}{
    \emph{Num} \tab \emph{Variable}	\tab \emph{Description}		\tab \emph{Concept}	\cr
    1   \tab glucose 		\tab Glucose (g/l) 		\tab physico-chemical 	\cr
    2   \tab fructose 		\tab Fructose (g/l) 		\tab physico-chemical 	\cr
    3   \tab saccharose 	\tab Saccharose (g/l) 		\tab physico-chemical 	\cr
    4   \tab sweet.power 	\tab Sweetening power (g/l) 	\tab physico-chemical 	\cr
    5   \tab ph1 		\tab pH before processing 	\tab physico-chemical 	\cr
    6   \tab ph2		\tab pH after centrifugation 	\tab physico-chemical 	\cr
    7   \tab titre		\tab Titre (meq/l) 		\tab physico-chemical 	\cr
    8   \tab citric.acid	\tab Citric acid (g/l)		\tab physico-chemical 	\cr
    9   \tab vitamin.c		\tab Vitamin C (mg/100g)	\tab physico-chemical 	\cr
    10  \tab smell.int		\tab Smell intensity 		\tab sensory		\cr
    11  \tab odor.typi		\tab Odor typicity 		\tab sensory		\cr
    12  \tab pulp		\tab Pulp 			\tab sensory		\cr
    13  \tab taste.int		\tab Taste intensity 		\tab sensory		\cr
    14  \tab acidity		\tab Acidity 			\tab sensory		\cr
    15  \tab bitter		\tab Bitterness			\tab sensory		\cr
    16  \tab sweet		\tab Sweetness			\tab sensory		\cr
    17  \tab judge1		\tab Ratings of judge 1 	\tab hedonic		\cr
    18  \tab judge2		\tab Ratings of judge 2 	\tab hedonic		\cr
    ... \tab ...		\tab ...			\tab ...		\cr
    112 \tab judge96		\tab Ratings of judge 96 	\tab hedonic		
  }
}
\source{
  Laboratoire de Mathematiques Appliques, Agrocampus, Rennes.
}
\references{
  Tenenhaus, M., Pages, J., Ambroisine, L., and Guinot, C. (2005) PLS methodology to study relationships between hedonic jedgements and product characteristics. \emph{Food Quality and Preference}, \bold{16}(4), pp. 315-325.

  Pages, J., and Tenenhaus, M. (2001) Multiple factor analysis combined with PLS path modelling. Application to the analysis of relationships between physicochemical, sensory profiles and hedonic judgements. \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{58}, pp. 261-273.

  Pages, J. (2004) Multiple Factor Analysis: Main Features and Application to Sensory Data. \emph{Revista Colombiana de Estadistica}, \bold{27}, pp. 1-26.
}
\examples{
  data(orange)
  orange
}
\keyword{datasets}
