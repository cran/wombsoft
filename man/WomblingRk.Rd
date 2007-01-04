\name{WomblingRk}
\alias{WomblingDominantRk}
\alias{WomblingCodominantRk}
\title{Matrices of the systemic function and the direction of the gradient.}
\description{These functions compute the wombling method (i.e. systemic function) with a moving bandwidth : the windows contain always the same number of individuals : k. The resulting matrix of the systemic function is written in the file systemic.txt, and the directions of the gradient in direction.txt. It also produces a coloured map showing the intensity of the gradient and its main directions.}
\usage{WomblingDominantRk(data,ranks,k,output_syst="syst.txt",output_dir="dir.txt")
WomblingCodominantRk(data,ranks,k,output_syst="syst.txt",output_dir="dir.txt")
}
\arguments{
\item{data}{the list of the six items provided by the Data and the Mirror functions.}
\item{ranks}{the vector of the ranks of distances of individuals to the grid points, result of the function Ranks.}
\item{k}{number of individuals in each window.}
\item{output_syst, output_dir}{file names of the two matrices. Default values are "syst.txt" and "dir.txt" in the current directory.}
}
\value{none}
\details{at each grid point, the weight assigned to an individual is \eqn{exp(-\frac{r_i}{k}}{exp(-Ri/k)}
where \eqn{r_i}{Ri} is the rank of the individual regarding the distances of individuals to the considered grid point, and k the number of individuals given by the user.}
\seealso{\code{\link{Data}} and \code{\link{Ranks}} for the format of the data and the ranks variables.}
\keyword{spatial}
\keyword{regression}
