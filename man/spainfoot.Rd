\name{spainfoot}
\alias{spainfoot}
\docType{data}
\title{Spanish football dataset}
\description{
  This data set contains the results of the teams in the Spanish football league 2008-2009.
}
\usage{spainfoot}
\format{
  A data frame with 20 observations on 14 variables. The variables may be used to construct four latent concepts: 1) \code{ATTACK}=Attack, 2) \code{DEFENSE}=Defense, 3) \code{SUCCESS}=Success, 4) \code{INDIS}=Indiscipline. 
  
    \tabular{llll}{
    \emph{Num} \tab \emph{Variable}  \tab \emph{Description}			\tab \emph{Concept}	\cr
    1   \tab GSH 	 \tab Goals Scored Home: total number of goals scored at home	\tab \code{ATTACK} 	\cr
    2   \tab GSA 	 \tab Goals Scored Away: total number of goals scored away	\tab \code{ATTACK} 	\cr
    3   \tab SSH 	 \tab Success to Score Home: Percentage of matches with scores goals at home	\tab \code{ATTACK} 		\cr
    4   \tab SSA 	 \tab Success to Score Away: Percentage of matches with scores goals away	\tab \code{ATTACK} 		\cr
    5   \tab GCH 	 \tab Goals Conceded Home: total number of goals conceded at home	\tab \code{DEFENSE} 		\cr    
    6   \tab GCA 	 \tab Goals Conceded Away: total number of goals conceded away	\tab \code{DEFENSE} 		\cr
    7   \tab CSH 	 \tab Clean Sheets Home: percentage of matches with no conceded goals at home	\tab \code{DEFENSE} 		\cr
    8   \tab CSA 	 \tab Clean Sheets Away: percentage of matches with no conceded goals away	\tab \code{DEFENSE} 		\cr    
    9   \tab WMH   \tab Won Matches Home: total number of won matches at home	\tab \code{SUCCESS} 		\cr    
    10  \tab WMA   \tab Won matches Away: total number of won matches away	\tab \code{SUCCESS} 		\cr    
    11  \tab LWR   \tab Longest Winning Run: longest run of won matches	\tab \code{SUCCESS} 		\cr    
    12  \tab LRWR  \tab Longest Run Without Loss: longest run of matches without losing	\tab \code{SUCCESS} 		\cr    
    13  \tab YC    \tab Yellow Cards: total number of yellow cards	\tab \code{INDIS} 		\cr    
    14  \tab RC    \tab Red Cards: total number of red cards	\tab \code{INDIS} 		\cr    
  }
}
\source{
  League Day. \url{http://www.leagueday.com} \cr
  Statto. \url{http://www.statto.com} \cr
  BDFutbol. \url{http://www.statto.com} \cr
  Cero a cero. \url{http://www.ceroacero.es/} \cr
}
\references{
  Getting started with plspm. \url{http://www.gastonsanchez.com/plspm}
}
\examples{
  data(spainfoot)
  spainfoot
}
\keyword{datasets}
