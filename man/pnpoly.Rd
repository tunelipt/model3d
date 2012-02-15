\name{pnpoly}
\alias{pnpoly}
\title{Is the point inside the the polygon?}
\usage{
  pnpoly(vx, vy, x, y)
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
  Determines if a point is inside a simple polygon.
}
\examples{
vx <- c(0, 2, 2, 0)
vy <- c(0, 0, 2, 2)
print(pnpoly(vx, vy, 1, 1))
print(pnpoly(vx, vy, 3, 1))
}

