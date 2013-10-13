\name{get_alpha}
\alias{get_alpha}
\title{Calculate Cronbach's alpha}
\usage{
  get_alpha(Block)
}
\arguments{
  \item{Block}{matrix (one block) of standardized manifest
  variables}
}
\value{
  Cronbach's alpha
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_alpha} is called by
  \code{get_unidim} and \code{alpha}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

