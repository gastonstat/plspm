\name{satisfaction}
\alias{satisfaction}
\docType{data}
\title{Satisfaction dataset}
\description{
  This data set contains the variables from a customer satisfaction study of a Spanish credit institution on 250 customers.
}
\usage{satisfaction}
\format{
  A data frame with 250 observations and 28 variables. Variables from 1 to 27 refer to six latent concepts: 1) \code{IMAG}=Image, 2) \code{EXPE}=Expectations, 3) \code{QUAL}=Quality, 4) \code{VAL}=Value, 5) \code{SAT}=Satisfaction, and 6) \code{LOY}=Loyalty. The last variable is a categorical variable indicating the gender of the individual.

    \code{IMAG}: Includes variables such as reputation, trustworthiness, seriousness, solidness, and caring about customer's needs. \cr
    \code{EXPE}: Includes variables such as products and services provided, customer service, providing solutions, and expectations for the overall quality. \cr
    \code{QUAL}: Includes variables such as reliable products and services, range of products and services, personal advice, and overall perceived quality. \cr
    \code{VAL}: Includes variables such as beneficial services and products, valuable investments, quality relative to price, and price relative to quality. \cr
    \code{SAT}: Includes variables such as overall rating of satisfaction, fulfillment of expectations, satisfaction relative to other banks, and performance relative to customer's ideal bank. \cr
    \code{LOY}: Includes variables such as propensity to choose the same bank again, propensity to switch to other bank, intention to recommend the bank to friends, and sense of loyalty. \cr
}
\source{
  Laboratory of Information Analysis and Modeling (LIAM). Facultat d'Informatica de Barcelona, Universitat Politecnica de Catalunya.
}
\examples{
  data(satisfaction)
  satisfaction
}
\keyword{datasets}
