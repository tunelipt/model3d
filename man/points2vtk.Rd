\name{points2vtk}
\alias{points2vtk}
\title{Function that stores point data in a vtk file.}
\usage{
  points2vtk(fname, x, y, z, ...)
}
\arguments{
  \item{fname}{Name of vtk file}

  \item{x}{x coordinates of the points}

  \item{y}{y coordinates of the points}

  \item{z}{z coordinates of the points}

  \item{...}{List of values associated with the points.}
}
\description{
  Create a vtk file (legacy format) from points. Values
  corresponding to the points can, optionally, be added as
  well.
}

