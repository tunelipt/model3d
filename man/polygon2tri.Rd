\name{polygon2tri}
\alias{polygon2tri}
\title{Decomposes a polygon in triangles.}
\usage{
  polygon2tri(p)
}
\arguments{
  \item{p}{Polygon}
}
\value{
  Array containing the triangles.
}
\description{
  Given a a simple convex polygon, this function decomposes
  it into triangles. The result is a mesh object (mesh3d
  for 3d polygons) that is a 3D array where the first index
  is the coordinate (x, y or z), the second the vertex (1,
  2, 3) and the third the triangle number.
}
\examples{
p <- newPolygon(c(0, 1, 2, 1, 0), c(0, 0, 1, 2, 2))
tri <- polygon2tri(p)
plot(tri[1,,], tri[2,,])
for (i in 1:dim(tri)[3]) polygon(tri[1,,i], tri[2,,i])
}

