points2vtk <- function(fname, x, y, z, ...){


  
  n <- length(x)
  tab <- cbind(x, y, z)
  
  con <- file(fname, open='w')
  cat('# vtk DataFile Version 2.0\n Points\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS', n, 'float\n', file=con)
  write.table(tab, file=con, row.names=FALSE, col.names=FALSE, sep=' ')
  cat("\nCELLS", n, n*2, "\n", file=con)
  count <- 0
  tab <- cbind(1, (1:n)-1)
  write.table(tab, file=con, row.names=FALSE, col.names=FALSE, sep=' ')
  
  
  cat("\nCELL_TYPES", n, "\n", file=con)
  cat(rep(1,n), file=con, sep='\n')

  dados <- list(...)
  nd <- length(dados)
  if (nd > 0){
    dnames <- names(dados)
    if (is.list(dados[[1]])) dados <- dados[[1]]
    if (is.null(names(dados))) names(dados) <- paste('V', 1:nd, sep='')
    dnames <- names(dados)
    cat("\nCELL_DATA", length(dados[[1]]), file=con)
    for (dn in dnames){
      cat("\nSCALARS", dn, "float 1\nLOOKUP_TABLE default\n", file=con)
      
      cat(dados[[dn]], sep='\n', file=con)
    }
      
  }

  close(con)
  

  
}

points2tec <- function(fname, x, y=NULL, z=NULL, ..., zname=NULL){

  
  dados <- list(...)
  nd <- length(dados)
  if (nd > 0){
    if (is.list(dados[[1]])) dados <- dados[[1]]
    if (is.null(names(dados))) names(dados) <- paste('V', 1:nd, sep='')
    dnames <- names(dados)
  }
  if (!is.null(z) && !is.null(y))
    xyz <- list(x=x, y=y, z=z)
  else if (!is.null(y))
    xyz <- list(x=x, y=y)
  else
    xyz <- list(x=x)

  if (nd == 0)
    dados <- xyz
  else
    dados <- c(xyz, dados)
  vars <- names(dados)
  header <- paste('VARIABLES =', paste(vars, collapse=' '))

  d <- dim(x)
  if (is.null(d)) d <- length(x)
  ndim <- length(d)
  if (is.null(zname)) zname <- 'Points'
  zname <- paste('ZONE T="',zname, '",', sep='')
  
  bheader <-paste(zname, paste(paste(c('I', 'J', 'K'), '=', d, sep='')[1:ndim], collapse=', '))
  nd <- length(dados)
  nr <- length(x)
  tab <- matrix(0, nr=length(x), nc=length(dados))
  for (i in 1:nd)
    tab[,i] <- as.double(dados[[i]])

  con <- file(fname, 'w')
  cat(header, '\n', file=con)
  cat(bheader, '\n', file=con)
  write.table(tab, file=con, row.names=FALSE, col.names=FALSE)
  close(con)
  


}
polygons2vtk <- function(fname, mesh, ...){
  

  npoly <- meshSize(mesh)
  x <- double(0)
  y <- double(0)
  z <- double(0)
  nvert <- integer(npoly)
  i <- 0
  for (k in 1:npoly){
    p <- meshGet(mesh, k)
    i <- i+1
    x <- c(x, p[1,])
    y <- c(y, p[2,])
    z <- c(z, p[3,])
    nvert[i] <- dim(p)[2]
  }
  
  npts <- length(x)
  totalverts <- sum(nvert)
  
  con <- file(fname, open='w')
  cat('# vtk DataFile Version 2.0\n Points\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS',
      npts, 'float\n', file=con)
  write.table(cbind(x, y, z), file=con, row.names=FALSE, col.names=FALSE, sep=' ')
  cat('\nCELLS', npoly, totalverts+npoly, '\n', file=con)
  count <- 0
  for (i in 1:npoly){
    cat(nvert[i], count + (0:(nvert[i]-1)), '\n', file=con)
    count <- count + nvert[i]
  }
  
  cat('\nCELL_TYPES', npoly, '\n', file=con)
  for (i in 1:npoly)
    cat(7, '\n', file=con)

  dados <- list(...)
  nd <- length(dados)
  if (nd > 0){
    dnames <- names(dados)
    if (is.list(dados[[1]])) dados <- dados[[1]]
    if (is.null(names(dados))) names(dados) <- paste('V', 1:nd, sep='')
    dn <- names(dados)
    cat("\nCELL_DATA", length(dados[[1]]), file=con)
    for (dn in dnames){
      cat("\nSCALARS", dn, "float 1\nLOOKUP_TABLE default\n", file=con)
      
      cat(dados[[dn]], sep='\n', file=con)
    }
      
  }

  close(con)
}
    
  
cmap2vtk <- function(fname, ...){
  cmlst <- list(...)
  nm <- length(cmlst)

  con <- file(fname, open='w')
  
  cat('<doc>\n', file=con)
  for (i in 1:nm){
    cm <- cmlst[[i]]
    cmname <- names(cmlst)[i]
    rgb <- col2rgb(cm)
    nc <- length(cm)
    x <- seq(0, 1, len=nc)
    r <- round(rgb[1,] / 255, 6)
    g <- round(rgb[2,] / 255, 6)
    b <- round(rgb[3,] / 255, 6)
    
    cat('<ColorMap name="', cmname, '" space="RGB">\n', file=con, sep='')
    for (k in 1:nc){
      xx <- round(x[k],6)
      cat('<Point x="', xx, '" o="', xx, '" r="', r[k], '" g="', g[k], '" b="', b[k], '"/>\n',
          file=con, sep='')
    }

    cat("</ColorMap>\n", file=con)
  }
  cat('</doc>\n', file=con)
  close(con)
}


