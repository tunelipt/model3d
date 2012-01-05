\name{voronoi2d}
\alias{voronoi2d}
\title{Computes the voronoi tesselation of a set of 2D points.}
\usage{
  voronoi2d(x, y, xb = NULL, yb = NULL)
}
\arguments{
  \item{x}{X coordinates of the points}

  \item{y}{Y coordinates of the points}

  \item{xb}{X limits of outer box of the tesselation}

  \item{yb}{Y limits of outer box of the tesselation}
}
\value{
  A list of polygons corresponding to the voronoi regions
  of each node.
}
\description{
  Uses the package \code{\link{deldir}} to compute the
  voronoi tesselation of a set of points.
}
\examples{
x <- runif(10)
y <- runif(10)
vor <- voronoi2d(x, y, c(0, 1), c(0, 1))
plot(x, y, ty='p', xlim=c(-.1, 1.1), ylim=c(-.1, 1.1))
for (p in vor) polygon(v[1,], v[2,])
}

