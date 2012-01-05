\name{polygonArea.polygon3d}
\alias{polygonArea.polygon3d}
\title{Calculates the area of a polygon.}
\usage{
  polygonArea.polygon3d(p)
}
\description{
  A method that computes the area of simple convex 3D
  polygons of class polygon3d.
}
\examples{
p <- newPolygon3d(c(0,1,2,1,2,0), c(0,0,1,2,3,3), c(0,0,0,0,0,0))
print(polygonArea(p))
}
\seealso{
  \code{\link{polygonArea}}
}

