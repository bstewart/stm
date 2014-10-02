writeLdac <- function(documents, file, zeroindex=TRUE) {
  if(missing(file)) stop("Please specify a filename with argument file.")
  lines <- lapply(documents, function(x) {
    Nd <-ncol(x)
    if(zeroindex) {
      paste(Nd, paste(x[1,]-1, x[2,], sep=":", collapse=" "))
    } else {
      paste(Nd, paste(x[1,], x[2,], sep=":", collapse=" "))
    }
  })
  lines <- unlist(lines)
  writeLines(lines, con=file)
}

