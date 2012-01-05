\name{getCoords}
\alias{getCoords}
\title{Returns the coordinates of a simple polygon.}
\usage{
  getCoords(p)
}
\arguments{
  \item{p}{Polygon}
}
\value{
  Matrix with coordinates of the vertices
}
\description{
  Returns a matrix where each row is a vertex of the
  polygon. Only works with simples polygons. It is actually
  a generic S3 function.
}

