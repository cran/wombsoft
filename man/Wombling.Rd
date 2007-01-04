\name{Wombling}
\alias{WomblingDominant}
\alias{WomblingCodominant}
\alias{WomblingContingency}
\title{Matrices of the systemic function and the direction of the gradient.}
\description{These functions compute the wombling method (i.e. systemic function). The resulting matrix of the systemic function is written in the file systemic.txt, and the directions of the gradient in direction.txt. It also produces a coloured map showing the intensity of the gradient and its main direction.}
\usage{WomblingDominant(data,h,output_syst="syst.txt",output_dir="dir.txt")
WomblingCodominant(data,h,output_syst="syst.txt",output_dir="dir.txt")
WomblingContingency(data,h,output_syst="syst.txt",output_dir="dir.txt")}
\arguments{
	\item{data}{the list of the six items provided by the Data and the Mirror functions.}
	\item{h}{bandwidth : corresponding to the scale we perform the analysis on.}
	\item{output_syst, output_dir}{file names of the two matrices. Default values are "syst.txt" and "dir.txt" in the current directory.}
}
\value{none}
\seealso{\code{\link{Data}} for the format of the data variable.}
\keyword{spatial}
\keyword{regression}