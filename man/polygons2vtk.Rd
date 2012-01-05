\name{polygons2vtk}
\alias{polygons2vtk}
\title{Function that stores polygons in a vtk file.}
\usage{
  polygons2vtk(fname, mesh, ...)
}
\arguments{
  \item{fname}{Name of vtk file}

  \item{mesh}{Polygonal mesh}

  \item{...}{List of values associated with the points.}
}
\description{
  Create a vtk file (legacy format) from polygons. Values
  corresponding to the points can, optionally, be added as
  well.
}

