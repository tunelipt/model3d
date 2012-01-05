\name{project3dTo2d}
\alias{project3dTo2d}
\title{Project a 3D face into a plane.}
\usage{
  project3dTo2d(p, base2d)
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
\seealso{
  \code{\link{faceProject2d}}
}

