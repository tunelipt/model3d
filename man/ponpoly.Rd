\name{ponpoly}
\alias{ponpoly}
\title{Is the point inside or on the boundary of a polygon?}
\usage{
  ponpoly(vx, vy, x, y, eps0 = 1e-08)
}
\arguments{
  \item{vx}{x coordinates of the vertices of the polygon.}

  \item{vy}{y coordinates of the vertices of the polygon.}

  \item{x}{x coordinate of the point to be tested.}

  \item{y}{y coordinate of the point to be tested.}
}
\value{
  TRUE if the point is inside the polygon.
}
\description{
  Determines if a point is inside or on the boundary of a
  simple polygon. This function simply calls
  \code{\link{pnpoly}} and \code{\link{pnpolybnd}}.
}
\examples{
vx <- c(0, 2, 2, 0)
vy <- c(0, 0, 2, 2)
print(pnpoly(vx, vy, 1, 1))
print(pnpoly(vx, vy, 3, 1))
print(pnpoly(vx, vy, 2, 1))
}

