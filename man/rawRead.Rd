\name{rawRead}
\alias{rawRead}
\title{Reads a raw file.}
\usage{
  rawRead(fname)
}
\arguments{
  \item{fname}{Name of the file to be read.}
}
\value{
  A list of faces that contain a 3D triangular mesh (mesh3d
  format).
}
\description{
  Raw files are a very simple file format to store 3D
  surface geometries. The surface geometries are decomposed
  into different surfaces called objects that store the
  surface as a triangulation of the surface.
}

