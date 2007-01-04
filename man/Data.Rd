\name{Data}
\alias{Data}
\alias{DataCodominant}
\alias{DataDominant}
\alias{DataContingency}
\title{Data preparation}
\description{Read a text file with coordinates and markers in columns and individuals in rows.}
\usage{DataDominant(input_file,conversion,nb_x,nb_y,output_coords="coord_km.txt")
DataCodominant(input_file,conversion,nb_x,nb_y,output_coords="coord_km.txt")
DataContingency(input_file,conversion,nb_x,nb_y,output_coords="coord_km.txt")}
\arguments{
	\item{input_file}{Path of the input text file.
  For dominant or codominant data, each row contains the name of the individual, the two coordinates (either abscissa and ordinates, or longitude and
latitude), and the genetic data in succession.
For contingency table, each row corresponds to a sampled point, with the name of the point, its coordinates, and the number of individuals for each modality of each variable.}
	\item{conversion}{0 if the coordinates are cartesians, 1 if they are in degree and therefore need to be converted to cartesians.}
	\item{nb_x,nb_y}{number of pixels in width and length of the grid.}
	\item{output_coords}{the name of the file where the kilometer coordinates will be saved in. Default value is "coord\_indiv.txt".}
}
\value{a list of six items :
	\item{spatial coordinates of individuals}{a matrix with one line per individual, and two columns containing abscissa and ordinates of individuals, (x,y).}
	\item{genetic_encoded}{the genetic data, containing one column per locus. If data are dominant, it's the same table as the input file. 
       }
  \item{grid}{a list of the vector of x, and the vector of y.}
  \item{cvx_vertices}{the	vertices of the convex hull of sampling area (same format than individuals coordinates).}
  \item{cvx_matrix}{a matrix containing a 1 if the corresponding point of the grid is in the convex hull, and a 0 otherwise.}
  \item{nb_individual}{the number of individuals in the dataset.}  
      }
\keyword{manip}
