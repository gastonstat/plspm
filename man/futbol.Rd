\name{futbol}
\alias{futbol}
\docType{data}
\title{Futbol dataset from Spain-England-Italy}
\description{
  This data set contains the results of the teams in the Spanish, English, and Italian football leagues 2009-2010 season.
}
\usage{data(futbol)}
\format{
  A data frame with 60 observations on the following 12 variables. The variables may be used to construct three latent concepts: 1) \code{ATTACK}=Attack, 2) \code{DEFENSE}=Defense, 3) \code{SUCCESS}=Success.

  \tabular{llll}{
    \emph{Num} \tab \emph{Variable}	\tab \emph{Description}		\tab \emph{Concept}	\cr
    1   \tab GSH: Goals Scored at Home 		\tab total number of goals scored at home		\tab ATTACK 	\cr
    2   \tab GSA: Goals Scored Away		\tab total number of goals scored away 			\tab ATTACK 	\cr
    3   \tab SSH: Success to Score at Home	\tab percentage of matches with scores goals at home 	\tab ATTACK 	\cr
    4   \tab SSA: Success to Score Away		\tab percentage of matches with scores goals away 	\tab ATTACK 	\cr
    5   \tab NGCH: Goals Conceded at Home 	\tab total number (negative) of goals conceded at home	\tab DEFENSE 	\cr
    6   \tab NGCA: Goals Conceded Away		\tab total number (negative) of goals conceded away 	\tab DEFENSE	\cr
    7   \tab CSH: Clean Sheets at Home		\tab percentage of matches with no conceded goals at home 	\tab DEFENSE 	\cr
    8   \tab CSH: Clean Sheets Away		\tab percentage of matches with no conceded goals away 	\tab DEFENSE 	\cr
    9   \tab WMH: Won Matches at Home		\tab total number of matches won at home		\tab SUCCESS 	\cr
    10  \tab WMA: Won Matches Away		\tab total number of matches won away			\tab SUCCESS 	\cr
    11  \tab Country: Leangue Country		\tab country of the team's league			\tab none 	\cr
    12  \tab Rank: Rank Position		\tab final ranking position within its league		\tab none 	\cr
  }
}
\source{
  League Day. \url{http://www.leagueday.com} \cr
  Statto. \url{http://www.statto.com} \cr  
}
\references{
  Applying plspm Beginners Guide. \url{http://www.scribd.com/plspm}
}
\examples{
  data(futbol)
  futbol
}
\keyword{datasets}
