\name{Mirror}
\alias{MirrorDominant}
\alias{MirrorCodominant}
\alias{MirrorContingency}
\title{Provides mirror effect on the data}
\description{Duplicates individuals that are close to the border of the
  convex hull (at a distance smaller than the one given by the user), on the position of their symmetric about the border.}
\usage{MirrorDominant(data,m)
MirrorCodominant(data,m) 
MirrorContingency(data,m)}
\arguments{
\item{data}{the list of the six items provided by the Data functions.}
\item{m}{width of the domain where individuals will be duplicated.}
  }
\value{a new data list, with the last number corresponding to the number of distinct original individuals (the new individuals provided by the mirror effect are not counted).}
\seealso{\code{\link{Data}} for the computing of the data list.}
\keyword{manip}
\keyword{spatial}