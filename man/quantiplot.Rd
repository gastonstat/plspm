\name{quantiplot}
\alias{quantiplot}
\title{Quantification Plot}
\usage{
  quantiplot(pls, lv = NULL, mv = NULL, pch = 16,
    col = "darkblue", lty = 2, ...)
}
\arguments{
  \item{pls}{a non-metric \code{"plspm"} object}

  \item{lv}{number or name of latent variable}

  \item{mv}{number or name of manifest variable}

  \item{pch}{Either an integer specifying a symbol or a
  single character to be used as the default in plotting
  points}

  \item{col}{color}

  \item{lty}{type of line}

  \item{\dots}{Further arguments passed on to
  \code{\link{plot}}.}
}
\description{
  Quantification Plots for Non-Metric PLS-PM
}
\details{
  If both \code{lv} and \code{mv} are specified, only the
  value of \code{lv} will be taken into account. \cr If the
  given \code{lv} have more than 15 variables, only the
  first 15 are plotted.
}

