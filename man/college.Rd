\name{college}
\alias{college}
\docType{data}
\title{College datasets}
\description{
  Dataset with different scores (high school, undergrad basic, undergrad intermediate, and GPA) of graduated college student in life sciences majors
}
\usage{data(college)}
\format{
  A data frame with 352 students on the following 13 variables. The variables may be used to construct four suggested latent concepts: 1) \code{HighSchool}=High School related scores, 2) \code{Basic}=scores of basic courses, 3) \code{InterCourse}=Scores of intermediate courses, 4) \code{GPA}=Final GPA (Graduate Point Average)

  \tabular{llll}{
    \emph{Num} \tab \emph{Variable}  \tab \emph{Description}    	\tab \emph{Concept}	\cr
    1   \tab \code{HS_GPA}      \tab High School GPA 			\tab \code{HighSchool} 	\cr
    2   \tab \code{SAT_Verbal} 	\tab Verbal SAT score 	  \tab \code{HighSchool} 	\cr
    3   \tab \code{SAT_Math}    \tab Math SAT score 	    \tab \code{HighSchool} 		\cr
    4   \tab \code{Biology1}    \tab Ecology Evolution 	  \tab \code{BasicCourses} 		\cr
    5   \tab \code{Chemistry1} 	\tab General Chemistry    \tab \code{BasicCourses} 		\cr    
    6   \tab \code{Math1}       \tab Calculus1	          \tab \code{BasicCourses} 		\cr
    7   \tab \code{Physics1} 	  \tab Introductory Physics	\tab \code{BasicCourses} 		\cr
    8   \tab \code{Biology2} 		\tab Cell and molecular biology \tab \code{InterCourses} 		\cr    
    9   \tab \code{Chemistry2}  \tab Organic Chemistry	  \tab \code{InterCourses} 		\cr
    10  \tab \code{Math2}       \tab Calculus 2	          \tab \code{InterCourses} 		\cr
    11  \tab \code{Physics2}    \tab Intermediate Physics	\tab \code{InterCourses} 		\cr
    12  \tab \code{FinalGPA}    \tab Graduation GPA       \tab \code{FinalGPA} 		\cr
    13  \tab \code{Gender}      \tab Gender               \tab \code{none} 		\cr
  }
}
\references{
  \url{http://www.gastonsanchez.com/plspm}
}
\examples{
# load data
data(college)

# take a peek
head(college)
}
\keyword{datasets}
