\name{get_numerics}
\alias{get_numerics}
\title{Transform factors in MV into numeric}
\usage{
  get_numerics(MV)
}
\arguments{
  \item{MV}{data frame with manifest variables as factors}
}
\value{
  list with transformed MV and categories of factors
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_numerics} is called by
  \code{get_treated_data}.
}
\note{
  This function is used when there are categorical
  variables defined as 'factors'. The idea is to convert
  those factors into numeric vectors while keeping the
  levels (ie categories) in a separate list
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

