\name{dotProduct}
\alias{dotProduct}
\title{Dot product of two vectors.}
\usage{
  dotProduct(u, v)
}
\arguments{
  \item{u}{First vector}

  \item{v}{Second vector}
}
\value{
  \eqn{\sum_{i=1}^{n} u_i v_i}
}
\description{
  Compute the dot product of two n-dimensional vectors.
}
\examples{
print(dotProduct(c(1,0,0), c(0,1,0)))
}
\seealso{
  \code{\link{crossProduct}}, \code{\link{vnorm}}
}

