\name{Binomial Tests}
\alias{BinomialTestCodominant}
\alias{BinomialTestDominant}
\alias{BinomialTestContingency}
\title{Test for boundaries using the binomial test}
\description{These functions take the value of the CandidateBoundaries functions as input, and test to see
 if the number of candidate boundary elements at each point of the grid is significant with regards to the
  binomial distribution. The pvalue is chosen by the user. The resulting matrices of the boundaries and their directions are written in text files. And a map is also produced, showing the boundaries in grey.}
\usage{BinomialTestDominant(data, cbe, pvalue = 0.05, output_bounds = "bounds.txt", output_dir = "direction.txt")
BinomialTestCodominant(data, cbe, pvalue = 0.05, output_bounds = "bounds.txt", output_dir = "direction.txt")
BinomialTestContingency(data, cbe, pvalue = 0.05, output_bounds = "bounds.txt", output_dir = "direction.txt")
}
\arguments{
  \item{data}{the list of the six items provided by the Data and the Mirror functions.}
  \item{cbe}{the candidate boundary elements, results of the functions CandidateBoundaries.}
  \item{pvalue}{pvalue of the binomial test. Default value is 0.05.}
  \item{output_bounds, output_dir}{file names for writing the two result matrices. Default values are "bounds.txt" and "direction.txt" in the current directory.}
}
\value{none}
\seealso{\code{\link{Data}} and \code{\link{CandidateBoundaries}} for the format of the data and the cbe variables.}
\keyword{spatial}