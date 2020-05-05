rowMergeMtx <- function(top,bottom) {
  
  trows = nrow(top)
  tcols = ncol(top)
  brows = nrow(bottom)
  bcols = ncol(bottom)
  
  stopifnot(trows == brows, tcols == bcols)
  mergemtx <- matrix(0.0, nrow=(trows+brows),ncol=tcols)
  
  for(j in 1:ncol(mergemtx)) {
    tidx = 1
    bidx = 1
    for(k in 1:nrow(mergemtx)) {
      if(k%%2==1) {
        mergemtx[k,j] = top[tidx,j]
        tidx = tidx + 1
      } else {
        mergemtx[k,j] = bottom[bidx,j]
        bidx = bidx + 1 
      }
    }
  }
  
  return(mergemtx)
}


rowSplitMtx <- function(m) {
  
  mrows = nrow(m)
  mcols = ncol(m)
  
  stopifnot(mrows%%2 == 0)
  
  top <- matrix(0.0, nrow=(mrows/2),ncol=mcols)
  bottom <- matrix(0.0, nrow=(mrows/2),ncol=mcols)
  
  for(j in 1:ncol(m)) {
    tidx = 1
    bidx = 1
    for(k in 1:nrow(m)) {
      if(k%%2==1) {
        top[tidx,j] = m[k,j]
        tidx = tidx + 1
      } else {
        bottom[bidx,j] = m[k,j]
        bidx = bidx + 1 
      }
    }
  }
  
  splitmtx <- list(top=top, bottom=bottom)
  return(splitmtx)
}