\name{polygonArea}
\alias{polygonArea}
\title{Calculates the area of a polygon.}
\usage{
  polygonArea(p)
}
\arguments{
  \item{p}{Matrix containing the vertices of the polygon.}
}
\value{
  Area of the polygon.
}
\description{
  A generic S3 function that computes the area of polygons
  given the coordinates of its vertices arranged as a
  matrix where each column is a vertex numbered in
  counterclockwise manner.
}
\examples{
poly <- newPolygon(c(0, 1, 1, 0), c(0, 0, 1, 1))
print(polygonArea(poly))
}
\seealso{
  \code{\link{crossProduct}}, \code{\link{poly3dArea}}
}

