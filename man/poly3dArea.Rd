\name{poly3dArea}
\alias{poly3dArea}
\title{Calculates the area of a 3D polygon.}
\usage{
  poly3dArea(p, n = NULL)
}
\arguments{
  \item{p}{Matrix containing the vertices of the polygon.}

  \item{n}{Normal to the polygon (if it has been previously
  computed).}
}
\value{
  Area of the polygon.
}
\description{
  Computes the area of a simple 3D polygon given the
  coordinates of its vertices arranged as a matrix where
  each column is a vertex numbered in counterclockwise
  manner.
}
\examples{
poly <- matrix(c(0,0,0, 1,0,0, 1,1,0, 0,1,0), 3, 4)
print(poly3dArea(poly))
}
\seealso{
  \code{\link{crossProduct}}, \code{\link{poly3dNormal}}
}

