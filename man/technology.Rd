\name{technology}
\alias{technology}
\docType{data}
\title{Technology data set}
\description{
  This data set contains the variables from a "user and acceptance of technology" model on 300 users.
}
\usage{data(technology)}
\format{
  A data frame with 300 observations and 21 variables. Variables can be grouped in six latent concepts: 1) \code{PERF_EXP}=Performance Expectancy, 2) \code{EFF_EXP}=Effort Expectancy, 3) \code{SUB_NORM}=Subjective Norm, 4) \code{FAC_COND}=Facilitating Conditions, 5) \code{BEH_INT}=Behavioral Intention, and 6) \code{USE_BEH}=Use Behavior. 
  
  \tabular{lll}{
    \emph{Num} \tab \emph{Variable}  \tab \emph{Description}		\cr
    1   \tab pe1 	\tab I find computers useful in my job  			\cr
    2   \tab pe2 	\tab Using computers in my job enables me to accomplish tasks more quickly  				\cr
    3   \tab pe3 	\tab Using computers in my job increases my productivity 			\cr
    4   \tab pe4	\tab Using computers enhances my effectiveness on the job				\cr
    5   \tab ee1 	\tab My interactions with computers are clear and understandable 			\cr
    6   \tab ee2	\tab It is easy for me to become skillful using computers  			\cr
    7   \tab ee3	\tab I find computers easy to use 			\cr
    8   \tab ee4	\tab Learning to use computers is easy for me \cr
    9   \tab sn1	\tab Most people who are important to me think I should use computers	\cr
    10  \tab sn2	\tab Most people who are important to me would want me to use computers			\cr
    11  \tab sn3	\tab People whose opinions I value would prefer me to use computers			\cr
    12  \tab fc1	\tab I have the resources and the knowledge and the ability to make use of the computer				\cr
    13  \tab fc1	\tab A central support was available to help with computer problems	\cr
    14  \tab fc1	\tab Management provided most of the necessary help and resources for computing				\cr
    15  \tab bi1  \tab I predict I will continue to use computers on a regular basis				\cr
    16  \tab bi2  \tab What are the chances in 100 that you will continue as a computer user?				\cr
    17  \tab bi3  \tab To do my work, I would use computers rather than any other means available				\cr
    18  \tab use1 \tab On an average working day, how much time do you spend using computers?				\cr
    19  \tab use2 \tab On average, how frequently do you use computers?				\cr
    20  \tab use3 \tab How many different computer applications have you worked with or used in your job?				\cr
    21  \tab use4 \tab According to your job requirements, indicate each task you use computers to perform?				\cr
  }    
}
\references{
  Venkatesh V., Morris M.G., Davis G.B., Davis F.D. (2003) User Acceptance of Information Technology: Toward a Unified View. \emph{MIS Quarterly}, Vol. 27 (3): 425-478.
}
\examples{
  data(technology)
  summary(technology)
}
\keyword{datasets}

