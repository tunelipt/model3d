\name{pnpolybnd}
\alias{pnpolybnd}
\title{Is the point in the boundary of a polygon?}
\usage{
  pnpolybnd(vx, vy, x, y, eps0 = 1e-08)
}
\arguments{
  \item{vx}{x coordinates of the vertices of the polygon.}

  \item{vy}{y coordinates of the vertices of the polygon.}

  \item{x}{x coordinate of the point to be tested.}

  \item{y}{y coordinate of the point to be tested.}

  \item{eps}{Admissible error.}
}
\value{
  TRUE if the point is on the segments.
}
\description{
  Determines if a point is on the boundary of a polygon.
}
\examples{
vx <- c(0, 2, 2, 0)
vy <- c(0, 0, 2, 2)
print(pnpolybnd(vx, vy, 1, 1))
print(pnpolybnd(vx, vy, 2, 1))
}

