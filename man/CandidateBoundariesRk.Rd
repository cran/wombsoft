\name{CandidateBoundariesRk}
\alias{CandidateBoundariesRk}
\alias{CandidateBoundariesCodominantRk}
\alias{CandidateBoundariesDominantRk}
\title{Number of alleles that contains a potential boundary at any point of the grid.}
\description{These functions count for each point of the grid, with a moving bandwidth, the number of candidate boundary elements at the given percentile pB at each point of the grid, and computes the overall direction.}
\usage{CandidateBoundariesCodominantRk(data, ranks, k, pB = 0.3)
CandidateBoundariesDominantRk(data, ranks, k, pB = 0.3)}
\arguments{
\item{data}{the list of the six items provided by the Data and the Mirror functions.}
\item{ranks}{the vector of the ranks of distances of individuals to the grid points, result of the function Ranks.}
\item{k}{number of individuals in the window.}
\item{pB}{percentile of points to be in the candidate boundary for each allele.}
}
\value{a list of three components : the first is the matrix of the number of potential bounds at each point, the second component is the sum of the corresponding directions, and the third is the threshold.}
\seealso{\code{\link{Data}} and \code{\link{Ranks}} for the format of the data and the ranks variables.}
\keyword{spatial}
\keyword{regression}
