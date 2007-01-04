\name{CandidateBoundaries}
\alias{CandidateBoundaries}
\alias{CandidateBoundariesCodominant}
\alias{CandidateBoundariesDominant}
\alias{CandidateBoundariesContingency}
\title{Number of alleles that contains a potential boundary at any point of the grid.}
\description{These functions count for each point of the grid, the number of candidate boundary elements at the given percentile pB at each point of the grid, and computes the overall direction.}
\usage{CandidateBoundariesCodominant(data, h, pB = 0.3)
CandidateBoundariesDominant(data, h, pB = 0.3)
CandidateBoundariesContingency(data, h, pB = 0.3)}
\arguments{
\item{data}{the list of the six items provided by the Data and the Mirror functions.}
\item{h}{bandwidth.}
\item{pB}{percentile of points to be in the candidate boundary for each allele.}
}
\value{a list of three components : the first is the matrix of the number of potential bounds at each point, the second component is the sum of the corresponding directions, and the third is the threshold.}
\seealso{\code{\link{Data}} for the format of the data variable.}
\keyword{spatial}
\keyword{regression}