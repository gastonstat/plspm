\name{get_ave}
\alias{get_ave}
\title{Calculate AVE (Average Variance Extracted)}
\usage{
  get_ave(communality)
}
\arguments{
  \item{communality}{list of communalities for each block}
}
\value{
  vector of average variance extracted for each block
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_ave} is called by
  \code{plspm}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

