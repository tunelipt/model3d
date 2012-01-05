\name{gpc2poly}
\alias{gpc2poly}
\title{Converts a gpclib polygon to a simple polygon.}
\usage{
  gpc2poly(gpc)
}
\arguments{
  \item{p}{A gpclib polygon}
}
\value{
  A polygon of class polygon
}
\description{
  This function Assumes that the polygon is simple and has
  no holes.
}
\examples{
p <- newPolygon(c(0,1,2,1,2,0), c(0,0,1,2,3,3))
g <- poly2gpc(p)
print(gpc2poly(g))
}
\seealso{
  \code{\link{newPolygon}}, \code{\link{poly2gpc}}
}

