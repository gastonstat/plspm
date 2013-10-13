\name{get_rank}
\alias{get_rank}
\title{Rank of a non-metric variable}
\usage{
  get_rank(X)
}
\arguments{
  \item{X}{an ordinal or nominal manifest variable to be
  ranked}
}
\value{
  vector with the corresponding rank for each category
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_rank} is called by
  \code{get_treated_data} and \code{get_weights}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

