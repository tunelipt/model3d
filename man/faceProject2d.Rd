\name{faceProject2d}
\alias{faceProject2d}
\title{Project a 3D face into a plane.}
\usage{
  faceProject2d(f, base2d)
}
\arguments{
  \item{f}{Face}

  \item{base2d}{Orthonormal local basis}
}
\value{
  Projection of f into base2d
}
\description{
  Generic function that projects a 3D face (a polygon or a
  mesh) into a plane characterized by 2 orhtonormal
  vectors.
}

