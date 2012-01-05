\name{crossProduct}
\alias{crossProduct}
\title{Cross product of two vectors.}
\usage{
  crossProduct(u, v)
}
\arguments{
  \item{u}{First vector}

  \item{v}{Second vector}
}
\value{
  u X v
}
\description{
  Compute the cross product of two 3D vectors.
}
\examples{
ex <- c(1, 0, 0)
ey <- c(0, 1, 0)
ez <- crossProduct(ex, ey)
}
\seealso{
  \code{\link{dotProduct}}
}

