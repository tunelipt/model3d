\name{project3dTo2d.default}
\alias{project3dTo2d.default}
\title{Method that projects a 3D points  into a plane.}
\usage{
  project3dTo2d.default(p, base2d)
}
\arguments{
  \item{f}{Matrix where each row is a 3D point}

  \item{base2d}{Orthonormal local basis}
}
\value{
  2D projection
}
\description{
  Method that projects a 3D points into a plane.
}
\seealso{
  \code{\link{faceProject2d}}, \code{\link{pointsProject}}
}

