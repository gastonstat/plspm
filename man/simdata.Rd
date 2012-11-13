\name{simdata}
\alias{simdata}
\docType{data}
\title{ Simulated data for REBUS with two groups}
\description{
  Simulated data with two latent classes showing different local models.  
}
\usage{data(simdata)}
\format{
  A data frame of simulated data with 400 observations on the following 14 variables.
  \describe{
    \item{\code{mv1}}{first manifest variable of the block \emph{Price Fairness} }
    \item{\code{mv2}}{second manifest variable of the block \emph{Price Fairness}}
    \item{\code{mv3}}{third manifest variable of the block \emph{Price Fairness}}
    \item{\code{mv4}}{fourth manifest variable of the block \emph{Price Fairness}}
    \item{\code{mv5}}{fifth manifest variable of the block \emph{Price Fairness}}
    \item{\code{mv6}}{first manifest variable of the block \emph{Quality} }
    \item{\code{mv7}}{second manifest variable of the block \emph{Quality}}
    \item{\code{mv8}}{third manifest variable of the block \emph{Quality}}
    \item{\code{mv9}}{fourth manifest variable of the block \emph{Quality}}
    \item{\code{mv10}}{fifth manifest variable of the block \emph{Quality}}
    \item{\code{mv11}}{first manifest variable of the block \emph{Customer Satisfaction}}
    \item{\code{mv12}}{second manifest variable of the block \emph{Customer Satisfaction}}
    \item{\code{mv13}}{third manifest variable of the block \emph{Customer Satisfaction}}
    \item{\code{group}}{a numeric vector}
  }
}
\details{
  The postulated model overlaps the one used by Jedidi \emph{et al.} (1997) and by 
  Esposito Vinzi \emph{et al.} (2007) for their numerical examples. 
  It is composed of one latent endogenous variable, \emph{Customer Satisfaction}, 
  and two latent exogenous variables, \emph{Price Fairness} and \emph{Quality}.
  Each latent exogenous variable (\emph{Price Fairness} and \emph{Quality}) has five
  manifest variables (reflective mode), and the latent endogenous variable 
  (\emph{Customer Satisfaction}) is measured by three indicators (reflective mode).
  
  Two latent classes showing different local models are supposed to exist.
  Each one is composed of 200 units. Thus, the data on the aggregate level 
  for each one of the numerical examples includes 400 units.
  
  The simulation scheme involves working with local models that are different 
  at both the measurement and the structural model levels. 
  In particular, the experimental sets of data consist of two latent classes with 
  the following characteristics:
  \cr
  (a) Class 1 - price fairness seeking customers - characterized by a strong relationship 
  between \emph{Price Fairness} and \emph{Customer Satisfaction} (close to 0.9) and a 
  weak relationship between \emph{Quality} and \emph{Customer Satisfaction} (close to 0.1),
  as well as by a weak correlation between the 3rd manifest variable 
  of the \emph{Price Fairness block} (mv3) and the corresponding latent variable;
  \cr
  (b) Class 2 - quality oriented customers - characterized by a strong relationship 
  between \emph{Quality} and \emph{Customer Satisfaction} (close to 0.1) and a weak 
  relationship between \emph{Price Fairness} and \emph{Customer Satisfaction} (close to 0.9), 
  as well as by a weak correlation between the 3rd manifest variable (mv8) of the \emph{Quality} 
  block and the corresponding latent variable.
}
\source{
  Simulated data from Trinchera (2007). See \bold{References} below.
}
\references{
Esposito Vinzi, V., Ringle, C., Squillacciotti, S. and Trinchera, L. (2007)
Capturing and treating unobserved heterogeneity by response based segmentation in
PLS path modeling. A comparison of alternative methods by computational experiments. 
\emph{Working paper}, ESSEC Business School.

Jedidi, K., Jagpal, S. and De Sarbo, W. (1997) 
STEMM: A general finite mixture structural equation model.
\emph{Journal of Classification} \bold{14}, pp. 23-50.

Trinchera, L. (2007)
Unobserved Heterogeneity in Structural Equation Models:  
a new approach to latent class detection in PLS Path Modeling. 
\emph{Ph.D. Thesis}, University of Naples "Federico II", Naples, Italy.
}
\examples{
data(simdata)
simdata
}
\keyword{datasets}
