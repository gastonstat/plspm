\name{offense}
\alias{offense}
\docType{data}
\title{Offense dataset}
\description{
  Dataset with offense statistics of American football teams from the NFL (2010-2011 season).
}
\usage{data(offense)}
\format{
  A data frame with 32 teams on the following 17 variables. The variables may be used to construct five suggested latent concepts: 1) \code{RUSH}=Rushing Quality, 2) \code{PASS}=Passing Quality, 3) \code{SPEC}=Special Teams and Other, 4) \code{SCORING}=Scoring Success, 5)\code{OFFENSE}=Offense Performance

  \tabular{llll}{
    \emph{Num} \tab \emph{Variable}  \tab \emph{Description}  		\tab \emph{Concept}	\cr
    1   \tab \code{YardsRushAtt}  \tab Yards per Rush Attempt 			\tab \code{RUSH} 	\cr
    2   \tab \code{RushYards} 		 \tab Rush Yards per game 	      \tab \code{RUSH} 	\cr
    3   \tab \code{RushFirstDown} \tab Rush First Downs per game 	\tab \code{RUSH} 		\cr
    4   \tab \code{YardsPassComp} \tab Yards Pass Completion 	    \tab \code{PASS} 		\cr
    5   \tab \code{PassYards} 	   \tab Passed Yards per game		    \tab \code{PASS} 		\cr    
    6   \tab \code{PassFirstDown} \tab Pass First Downs per game	  \tab \code{PASS} 		\cr
    7   \tab \code{FieldGoals} 	 \tab Field Goals per game		    \tab \code{SPEC} 		\cr
    8   \tab \code{OtherTDs} 		 \tab Other Touchdowns (non-offense) per game \tab \code{SPEC} 		\cr    
    9   \tab \code{PointsGame}    \tab Points per game	            \tab \code{SCORING} 		\cr
    10  \tab \code{OffensTD}      \tab Offense Touchdowns per game	\tab \code{SCORING} 		\cr
    11  \tab \code{TDGame}        \tab Touchdowns per game	        \tab \code{SCORING} 		\cr
    12  \tab \code{PassTDG}       \tab Passing Touchdowns per game          \tab \code{OFFENSE} 		\cr
    13  \tab \code{RushTDG}        \tab Rushing Touchdowns per game          \tab \code{OFFENSE} 		\cr
    14  \tab \code{PlaysGame}      \tab Plays per game          \tab \code{OFFENSE} 		\cr
    15  \tab \code{YardsPlay}      \tab Yards per Play          \tab \code{OFFENSE} 		\cr
    16  \tab \code{FirstDownPlay}  \tab First Downs per Play    \tab \code{OFFENSE} 		\cr
    17  \tab \code{OffTimePossPerc}\tab Offense Time Possession Percentage  \tab \code{OFFENSE} 		\cr
  }
}
\source{
 \url{http://www.teamrankings.com/nfl/stats/}
}
\examples{
# load data
data(offense)

# take a peek
head(offense)
}
\keyword{datasets}
