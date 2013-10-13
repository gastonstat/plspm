\name{get_inner_summary}
\alias{get_inner_summary}
\title{Inner summary assessment}
\usage{
  get_inner_summary(path_matrix, blocks, modes,
    communality, redundancy, R2)
}
\arguments{
  \item{path_matrix}{matrix of path connections}

  \item{blocks}{list indicating blocks of manifest
  variables}

  \item{modes}{vector of modes}

  \item{communality}{vector with communality values}

  \item{redundancy}{vector with redundancy values}

  \item{R2}{vector with R2 values}
}
\value{
  A data frame with the following columns:

  \item{Type}{Exogenous or Endogenous}

  \item{R2}{R2 coefficient}

  \item{Mean_Communality}{average communality}

  \item{Mean_Redundancy}{average redudancy}

  \item{AVE}{Average Variance Extracted}
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_inner_summary} is called by
  \code{plspm}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

