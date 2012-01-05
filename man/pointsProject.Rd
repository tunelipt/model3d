\name{pointsProject}
\alias{pointsProject}
\title{Projects 3D points into a plane.}
\usage{
  pointsProject(p, base2d)
}
\arguments{
  \item{p}{Points}

  \item{base2d}{Local plane coordinate system}
}
\value{
  a matrix where each column is a projected 2D point.
}
\description{
  3D points arranged as a matrix where each column
  corresponds to a point is projected into a plane that has
  to direction unit vectors x and y given by the columns of
  \code{base2d}.
}

