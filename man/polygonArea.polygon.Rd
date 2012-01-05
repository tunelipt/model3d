\name{polygonArea.polygon}
\alias{polygonArea.polygon}
\title{Calculates the area of a polygon.}
\usage{
  polygonArea.polygon(p)
}
\description{
  A method that computes the area of simple polygons of
  class polygon.
}
\examples{
poly <- newPolygon(c(0, 1, 1, 0), c(0, 0, 1, 1))
print(polygonArea(poly))
}
\seealso{
  \code{\link{polygonArea}},
  \code{\link{polygonArea.gpc.poly}}
}

