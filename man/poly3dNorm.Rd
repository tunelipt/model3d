\name{poly3dNorm}
\alias{poly3dNorm}
\title{Normal of a 3D polygon.}
\usage{
  poly3dNorm(p, eps = 1e-08)
}
\arguments{
  \item{p}{Polygon where each column are the coordinates of
  a vertex}

  \item{eps}{Relative error acceptable when determining if
  3 vertices are colinear}
}
\value{
  A vector containing the normal to the polygon when the
  vertices are assumed to be numbered in counterclockwise
  manner.
}
\description{
  Computes the normal of a simple polygon.
}
\examples{
poly <- matrix(c(0,0,0, 1,0,0, 1,1,0, 0,1,0), 3, 4)
print(poly3dNorm(poly))
}
\seealso{
  \code{\link{crossProduct}}, \code{\link{triNormal}}
}

