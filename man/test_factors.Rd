\name{test_factors}
\alias{test_factors}
\title{Test presence of factors}
\usage{
  test_factors(DF)
}
\arguments{
  \item{DF}{data frame that might contain factors}
}
\value{
  whether DF contains factors or not
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{test_factors} is called by
  \code{get_treated_data}.
}
\note{
  Test presence of factors in a data.frame
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

