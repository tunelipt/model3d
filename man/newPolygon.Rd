\name{newPolygon}
\alias{newPolygon}
\title{Create a new simple 2D polygon.}
\usage{
  newPolygon(x, y)
}
\arguments{
  \item{x}{Vector with x coordinates of the vertices}

  \item{y}{Vector with y coordinates of the vertices}
}
\value{
  A matrix with class \code{polygon} where each column is a
  vertex of the polygon
}
\description{
  Creates a new simple 2D polygon from the coordinates. The
  polygon is stored as a matrix where the first row are the
  x coordinates of the vertices and the second row are the
  y coordinates.
}
\examples{
p <- newPolygon(c(0,1,2,1,2,0), c(0,0,1,2,3,3))
print(p)
}
\seealso{
  \code{\link{newPolygon3d}}, \code{\link{polygonArea}},
  \code{\link{polygonCentroid}},
  \code{\link{polygonIntersect}}
}

