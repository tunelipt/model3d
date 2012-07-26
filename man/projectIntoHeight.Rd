\name{projectIntoHeight}
\alias{projectIntoHeight}
\title{Chops a polygon at two different heights.}
\usage{
  projectIntoHeight(p, hmin, hmax, dnum = 3)
}
\arguments{
  \item{p}{Polygon, a matrix where the rows are x, y, z
  coordinates of the vertices.}

  \item{hmin}{Minimum height to chop the polygon.}

  \item{hmax}{Maximum height to chop the polygon.}

  \item{dnum}{Axis along which the heights are defined.}
}
\value{
  Chopped polygon or NULL if there is no intersection.
}
\description{
  Given a polygon and two levels, this function chops the
  polygon, returning a polygon that is the intersection
  between the original polygon and a strip from the lower
  height to the upper height. If there is not intersection
  this function returns \code{NULL}. The height can be any
  coordinate axis: x -> 1, y -> 2 and z -> 3 (default).
}

