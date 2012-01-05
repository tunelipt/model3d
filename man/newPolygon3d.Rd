\name{newPolygon3d}
\alias{newPolygon3d}
\title{Create a new simple 3D polygon.}
\usage{
  newPolygon3d(x, y, z)
}
\arguments{
  \item{x}{Vector with x coordinates of the vertices}

  \item{y}{Vector with y coordinates of the vertices}

  \item{z}{Vector with z coordinates of the vertices}
}
\value{
  A matrix with class \code{polygon3d} where each column is
  a vertex of the polygon
}
\description{
  Creates a new simple 3D polygon from the coordinates. The
  polygon is stored as a matrix where the first row are the
  x coordinates of the vertices and the second row are the
  y coordinates anf the third row are the z coordinates.
}
\examples{
p <- newPolygon3d(c(0,1,2,1,2,0), c(0,0,1,2,3,3), c(0,0,0,0,0,0))
print(p)
}
\seealso{
  \code{\link{newPolygon}}, \code{\link{polygonArea}},
  \code{\link{polygonCentroid}}
}

