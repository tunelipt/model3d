\name{points2tec}
\alias{points2tec}
\title{Function that stores point data in a tecplot file.}
\usage{
  points2tec(fname, x, y = NULL, z = NULL, ...,
    zname = NULL)
}
\arguments{
  \item{fname}{Name of vtk file}

  \item{x}{x coordinates of the points}

  \item{y}{y coordinates of the points}

  \item{z}{z coordinates of the points}

  \item{...}{List of values associated with the points.}

  \item{zname}{Zone name}
}
\description{
  Create a vtk file (legacy format) from points. Values
  corresponding to the points can, optionally, be added as
  well.
}

