\name{polygonArea.gpc.poly}
\alias{polygonArea.gpc.poly}
\title{Calculates the area of a polygon.}
\usage{
  polygonArea.gpc.poly(p)
}
\description{
  A method that computes the area of polygons of class
  gpc.poly.
}
\examples{
p <- newPolygon(c(0,1,2,1,2,0), c(0,0,1,2,3,3))
g <- poly2gpc(p)
print(polygonArea(g)
}
\seealso{
  \code{\link{polygonArea}},
  \code{\link{polygonArea.polygon}}
}

