rawGetTri <- function(stri){
  nl <- length(stri)
  
  con <- textConnection(stri)
  on.exit(close(con))
  x <- scan(con, quiet=TRUE)
  dim(x) <- c(3, 3, nl)
  class(x) <- 'mesh3d'
  return(x)

}  


rawRead <- function(fname){

  s <- readLines(fname)

  header <- grep("Object[0-9]+", s)

  ntri <- length(header)

  lstart <- header+1
  lend <- c((header-1)[2:ntri], length(s))
  #return(list(s, lstart, lend))
  
  tri <- lapply(1:ntri, function(i) rawGetTri(s[lstart[i]:lend[i]]))
  n <- length(tri)
  names(tri) <- paste('F', 1:n, sep='')
  return(tri)
}
  

raw2vtk <- function(tri, fname){
  tri <- joinTri(tri)

  d <- dim(tri)
  n <- d[3]

  con <- file(fname, open='w')
  cat('# vtk DataFile Version 2.0\n Triangles\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS', n*3, 'float\n', file=con)
  for (i in 1:n)
    for (k in 1:3)
      cat(tri[,k,i], "\n", file=con)
  cat("\nCELLS", n, n*4, "\n", file=con)
  count <- 0
  for (i in 1:n){
    cat(3, count:(count+2), "\n", file=con)
    count <- count + 3
  }
  
  cat("\nCELL_TYPES", n, "\n", file=con)
  for (i in 1:n)
    cat(5, "\n", file=con)
  close(con)
  
  
}
