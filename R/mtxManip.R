rowMergeMtx <- function(a,b) {
  
  arows = nrow(a)
  acols = ncol(a)
  brows = nrow(b)
  bcols = ncol(b)
  
  if(arows != brows || acols != bcols) {
    return(FALSE)
  }
  
  mrgmtx <- matrix(0.0, nrow=(arows+brows),ncol=acols)
  
  for(j in 1:ncol(mrgmtx)) {
    aidx = 1
    bidx = 1
    for(k in 1:nrow(mrgmtx)) {
      if(k%%2==1) {
        mrgmtx[k,j] = a[aidx,j]
        aidx = aidx + 1
      } else {
        mrgmtx[k,j] = b[bidx,j]
        bidx = bidx + 1 
      }
    }
  }
  
  return(mrgmtx)
}