\name{polygonProject}
\alias{polygonProject}
\title{Project a 3D polygon into a plane.}
\usage{
  polygonProject(p, base2d)
}
\arguments{
  \item{p}{3D polygon to be projected.}

  \item{base2d}{Matrix where each column corresponds to the
  local x and y dirs.}
}
\value{
  A 2D polygon.
}
\description{
  Projects all vertices of a 3D polygon into a plane
  characterized by 2 unit direction vectors.
}

