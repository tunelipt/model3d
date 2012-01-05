\name{triNormal}
\alias{triNormal}
\title{Compute the normal to a triangle.}
\usage{
  triNormal(x, invert = FALSE)
}
\arguments{
  \item{x}{Matrix 3x3 with each column as a vertex.}

  \item{invert}{A boolean that specifies whether the outer
  normal is wanted.}
}
\value{
  A vector containing the 3 components of the normal.
}
\description{
  Compute the normal to the plane formed by the edges of a
  triangle.
}
\examples{
tri <- matrix(c(0,0,0, 1,0,0, 1,1,0), 3, 3)
print(triNormal(tri))
}
\seealso{
  \code{\link{crossProduct}}
}

