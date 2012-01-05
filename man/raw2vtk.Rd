\name{raw2vtk}
\alias{raw2vtk}
\title{Converts a raw mesh to vtk file format.}
\usage{
  raw2vtk(tri, fname)
}
\arguments{
  \item{A}{list of mesh3d objects.}

  \item{fname}{Name of vtk file to save the data.}
}
\description{
  The raw mesh is actually a list of mesh3d that was read
  with rawRead. This function saves this list as vtk file
  (legacy format).
}

