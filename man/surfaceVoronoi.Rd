\name{surfaceVoronoi}
\alias{surfaceVoronoi}
\title{Computes the voronoi diagram of points in a curved surface.}
\usage{
  surfaceVoronoi(m, pts, basefun = getBase)
}
\arguments{
  \item{m}{Polygonal 3D mesh.}

  \item{x}{X Coordinate of the points that will be used to
  compute the Voronoi diagram}

  \item{y}{Y Coordinate of the points that will be used to
  compute the Voronoi diagram}

  \item{z}{Z Coordinate of the points that will be used to
  compute the Voronoi diagram}

  \item{pnames}{Names of each point.}

  \item{basefun}{Function used to obtain the projection
  plane.}
}
\value{
  A list containing 3D polygonal meshes of the influence
  area for each node.
}
\description{
  Computing voronoi diagrams for 2D points is simple if
  libraries are available. If the points are contrained in
  a curved 3D surface, things are a little more
  complicated. This function projects the surface on a
  plane computes the voronoi diagram in this plane and then
  projects the resulting polygons back to 3D.
}

