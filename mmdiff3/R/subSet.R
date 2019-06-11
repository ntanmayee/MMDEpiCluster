#' @return \code{setSet(x,idx)} returns a DBAmmd Object
#' with a subset of peak regions.
#' @rdname DBAmmd-Accessors
setMethod("subSet", "DBAmmd", function(x,Peak.IDs){

  idx <- match(Peak.IDs, names(Regions(x)))

  x@rowRanges = x@rowRanges[idx]
  if (length(x@Reads)>0){
    R <- x@Reads
    for (i in 1:length(R)){
      for (j in 1:length(R[[i]])){
        R[[i]][[j]] <- R[[i]][[j]][idx]
      }
    }
    x@Reads <- R
  }
  if (length(x@RawTotalCounts)>0){
    x@RawTotalCounts <- x@RawTotalCounts[idx,]
  }
  if (length(x@RawCounts.p)>0){
    x@RawCounts.p  <- x@RawCounts.p[idx,]
  }
  if (length(x@RawCounts.n)>0){
    x@RawCounts.n  <- x@RawCounts.n[idx,]
  }
  H <- x@Hists
  if (length(H)>0){
    for (i in 1:length(H)){
      H[[i]] <- H[[i]][idx]
    }
    x@Hists <- H
  }
  D <- x@DISTs
  if (length(D)>0){
    for (i in 1:length(D)){
      D[[i]] <- D[[i]][idx,]
    }
    x@DISTs <- D
  }

  if (length(x@mCounts)>0){
    x@mCounts <-  x@mCounts[idx,]
  }

  C <-  x@Contrasts
  if (length(C)>0){
    for (i in 1:length(C)){
      if (!is.null(C[[i]]$x.locfit)){
        C[[i]]$MMD.locfit$group1 <- C[[i]]$MMD.locfit$group1[idx]
        C[[i]]$MMD.locfit$group2 <- C[[i]]$MMD.locfit$group2[idx]
        C[[i]]$MMD.locfit$between.gr1 <- C[[i]]$MMD.locfit$between.gr1[idx]
        C[[i]]$MMD.locfit$between.gr2 <- C[[i]]$MMD.locfit$between.gr2[idx]
        C[[i]]$MMD.locfit$MMD.Table <- C[[i]]$MMD.locfit$MMD.Table[idx,]
        C[[i]]$MMD.locfit$outlier <- C[[i]]$MMD.locfit$outlier[idx]
      }
    }
  }
  x@Contrasts <- C
  return(x)
})

