\name{get_scores}
\alias{get_scores}
\title{Calculate Latent Variable Scores}
\usage{
  get_scores(X, W)
}
\arguments{
  \item{X}{matrix of scaled/quantified manifest variables}

  \item{W}{matrix with outer weights (ODM: Outer Design
  Matrix)}
}
\value{
  Scores matrix of latent variables scores
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_scores} is called by
  \code{plspm}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

