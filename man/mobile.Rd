\name{mobile}
\alias{mobile}
\docType{data}
\title{ ECSI Mobile Phone Provider dataset}
\description{
  This table contains data from the article by Tenenhaus \emph{et al.} (2005), see reference below.
}
\usage{data(mobile)}
\format{
  A data frame with 250 observations on 24 variables on a scale from 0 to 100. The variables refer to seven latent concepts: 1) \code{IMAG}=Image, 2) \code{EXPE}=Expectations, 3) \code{QUAL}=Quality, 4) \code{VAL}=Value, 5) \code{SAT}=Satisfaction, 6) \code{COM}=Complaints, and 7) \code{LOY}=Loyalty. 

  \code{IMAG}: Includes variables such as reputation, trustworthiness, seriousness, and caring about customer's needs. \cr
  \code{EXPE}: Includes variables such as products and services provided and expectations for the overall quality. \cr
  \code{QUAL}: Includes variables such as reliable products and services, range of products and services, and overall perceived quality. \cr
  \code{VAL}:  Includes variables such as quality relative to price, and price relative to quality. \cr
  \code{SAT}:  Includes variables such as overall rating of satisfaction, fulfillment of expectations, satisfaction relative to other phone providers.\cr
  \code{COM}:  Includes one variable defining how well or poorly custmer's complaints were handled. \cr
  \code{LOY}:  Includes variables such as propensity to choose the same phone provider again, propensity to switch to other phone provider, intention to recommend the phone provider to friends. \cr
  
  \describe{
    \item{\code{ima1}}{ First MV of the block Image}
    \item{\code{ima2}}{ Second MV of the block Image}
    \item{\code{ima3}}{ Third MV of the block Image}
    \item{\code{ima4}}{ Fourth MV of the block Image}
    \item{\code{ima5}}{ Fifth MV of the block Image}
    \item{\code{exp1}}{ First MV of the block Expectations}
    \item{\code{exp2}}{ Second MV of the block Expectations}
    \item{\code{exp3}}{ Third MV of the block Expectations}
    \item{\code{qua1}}{ First MV of the block Quality}
    \item{\code{qua2}}{ Second MV of the block Quality}
    \item{\code{qua3}}{ Third MV of the block Quality}
    \item{\code{qua4}}{ Fourth MV of the block Quality}
    \item{\code{qua5}}{ Fifth MV of the block Quality}
    \item{\code{qua6}}{ Sixth MV of the block Quality}
    \item{\code{qua7}}{ Seventh MV of the block Quality}
    \item{\code{val1}}{ First MV of the block Value}
    \item{\code{val2}}{ Second MV of the block Value}
    \item{\code{sat1}}{ First MV of the block Satisfaction}
    \item{\code{sat2}}{ Second MV of the block Satisfaction}
    \item{\code{sat3}}{ Third MV of the block Satisfaction}
    \item{\code{comp}}{ First MV of the block Complaints}
    \item{\code{loy1}}{ First MV of the block Loyalty}
    \item{\code{loy2}}{ Second MV of the block Loyalty}
    \item{\code{loy3}}{ Third MV of the block Loyalty}
  }
}
\references{
  Tenenhaus, M., Esposito Vinzi, V., Chatelin Y.M., and Lauro, C. (2005) PLS path modeling. \emph{Computational Statistics & Data Analysis}, \bold{48}, pp. 159-205. 
}
\examples{
data(mobile)
}
\keyword{datasets}
