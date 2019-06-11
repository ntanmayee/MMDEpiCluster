#' report results
#'
#' retrieve results of differential binding analysis
#'
#' @inheritParams DBAmmd-Accessors
#' @inheritParams compPvals
#' @inheritParams plotDists
#' @inheritParams plotPeak
#'
#' @param rm.oulier if TRUE, significant peaks with high
#' within-group distances are not reported. (DEFAULT: TRUE)
#' @examples
#'
#' data("MMD")
#' res <- reportResults(MMD)
#'
#' @export

reportResults <- function(MD, diff.method='MMD.locfit', th=0.1,Peak.id=NULL,
                          whichContrast=1,rm.oulier=TRUE,bUsePval=FALSE){
  
  if (is.null(Peak.id)){
    c <- Contrast(MD,whichContrast)
    c.m <- c[[diff.method]]
    if (is.null(c.m)){
      warning('No results found for this contrast')
      res <- Regions(MD)[NULL]
      return(res)
    }
    if (diff.method=='MMD.locfit'){
      if (!bUsePval){
        idx.sig <- which(c.m$mmd.Table[,'padj'] < th)
      } else{
        idx.sig <- which(c.m$mmd.Table[,'pval'] < th)
      }
      if (rm.oulier){
        i.o <- which(c.m$outlier==TRUE)
        idx.sig <- setdiff(idx.sig,i.o)
      }
      res <- Regions(MD)[idx.sig]
      tab <- c.m$mmd.Table[idx.sig,]
      ii <- sort(tab[,'padj'],index.return = TRUE)
      res <- res[ii$ix]
      tab <- tab[ii$ix,]
      res$pval <- tab[,'pval']
      res$padj <- tab[,'padj']
    }
  } else{
    Peak.idx <- match(Peak.id, names(Regions(MD)))
    res <- matrix(0,length(MD@Contrasts),6)
    rownames(res) = seq(1,length(MD@Contrasts))
    colnames(res) = c('ID','pvals','padj','log2BvsW','ranklog2BvsW','outlier')
    for (i in 1:length(MD@Contrasts)){
      C <- MD@Contrasts[[i]]
      name <- paste(C$name1,' vs ',C$name2,sep='')
      if (is.null(C[[diff.method]])){
        next
      }
      res[i,1:5] = C$MMD.locfit$mmd.Table[Peak.idx,]
      res[i,6] = C$MMD.locfit$outlier[Peak.idx]
    rownames(res)[i] = name
    }
    
  }
  
  return(res)
}