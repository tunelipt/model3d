\name{poly2gpc}
\alias{poly2gpc}
\title{Converts a polygon to a gpclib polygon.}
\usage{
  poly2gpc(p)
}
\arguments{
  \item{p}{A polygon.}
}
\value{
  A gpclib polygon.
}
\description{
  The R package gpclib provides a few algorithms that are
  useful but it has an internal (very general)
  representation of polygons. This function converts a
  simple polygon of class polygon to a gpclib polygon.
}
\examples{
p <- newPolygon(c(0,1,2,1,2,0), c(0,0,1,2,3,3))
g <- poly2gpc(p)
plot(g)
}
\seealso{
  \code{\link{newPolygon}}, \code{\link{gpc2poly}}
}

