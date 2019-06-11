#' compute p-values
#'
#' This function determines peak-specific p-values based on distances between
#' sample histograms.
#'
#' @inheritParams compDists
#' @param diff.method method used to determine p-values and false discovery rates.
#' Currently only 'MMD.locfit' implemented.
#' (DEFAULT: 'MMD.locfit')
#'
#' @return DBAmmd object with updated \code{Contrasts} slot.
#'
#' @examples
#'
#' ## Example using a small data set provided with this package:
#'
#' data("MMD")
#' MMD.1 <- setContrast(MMD,contrast='byCondition')
#' MMD.1 <- compPvals(MMD.1)
#' reportResults(MMD.1)
#'
#' @seealso \code{\link{DBAmmd}},\code{\link{reportResults}},
#' \code{\link{plotDists}},\code{\link{compDists}}
#'
#'
#' @import Biobase locfit
#' @importFrom stats p.adjust residuals predict pchisq pt
#' @export
#
#
# subfunctions:
#         - compute.pvalues
#         - determineGroupComps
#         - findComps
#
# Gabriele Schweikert
# March 2016

compPvals <- function(MD, dist.method='MMD',diff.method='MMD.locfit',
                      plot.on=FALSE,verbose=1){

  if (missing(MD))
    stop("DBAmmd is missing")


  Meta <- metaData(MD)

  DISTs <- Dists(MD,dist.method)
  mCounts <- MD@mCounts
  contrasts <- MD@Contrasts
  if (is.null(contrasts)) {
    stop('contrasts missing')
  }


  for (c in 1:length(contrasts)){
    contrast <- contrasts[[c]]
    group1 <- names(which(contrasts[[c]]$group1))
    group2 <- names(which(contrasts[[c]]$group2))

    compNames <- colnames(DISTs)

    within1 <- determineGroupComps(group1,type='within')
    within2 <- determineGroupComps(group2,type='within')
    between <- determineGroupComps(group1,group2,type='between')

    if (is.null(within1)&&is.null(within2)){
      warning('no replicates given to estimate within group variance')
      next()
    }
    if (is.null(within1)){
      warning('no replicates given for group1')}
    if (is.null(within2)){
      warning('no replicates given for group2')}

    group1.ids <- findComps(within1,compNames,dist.method)
    group2.ids <- findComps(within2,compNames,dist.method)
    bw.group.ids <- findComps(between,compNames,dist.method)

    message("Computing p-values \n")
    if (verbose==2){
        message("group1.ids:")
        print(group1.ids)
        message("group2.ids:")
        print(group2.ids)
        message("bw.group.ids:")
        print(bw.group.ids )
    }
    if (diff.method=='MMD.locfit'){
      P <- compute.pvalues(mCounts, MMD=DISTs, group1.ids, group2.ids,
                           bw.group.ids,plot.on=plot.on)
      contrast$MMD.locfit <- list(fit.Group1=P$fit.Group1, fit.Group2=P$fit.Group2,
                                  fit.disp.Group1=P$fit.disp.Group1,
                                  fit.disp.Group2=P$fit.disp.Group2,
                                  group1=P$group1,
                                  group2=P$group2,
                                  between.gr1=P$between.gr1,
                                  between.gr2=P$between.gr2,
                                  mmd.Table=P$mmd.Table,
                                  outlier=P$outlier)
    }
    contrasts[[c]] <- contrast
  }
  MD@Contrasts=contrasts
  return(MD)
}



compute.pvalues <- function(Counts, MMD, group1.ids, group2.ids, bw.group.ids,
                            plot.on=FALSE ){


  Counts <- log10(Counts+1)

  mmd.Table <- matrix(1,dim(Counts),5)
  colnames(mmd.Table) <- c('id','pval','padj','log2BvsW','ranklog2BvsW')
  mmd.Table[,'id'] <- seq(1,dim(Counts)[1])

  W1.y <- as.matrix(MMD[,group1.ids])
  W2.y <- as.matrix(MMD[,group2.ids])
  W.y <- rowMeans(cbind(W1.y,W2.y))
  B.y <- rowMeans(as.matrix(MMD[,bw.group.ids]))
  mmd.Table[,'log2BvsW'] <- log2(B.y/W.y)
  mmd.Table[,'ranklog2BvsW'] <- rank(-log2(B.y/W.y))

  if (!is.null(group1.ids)){
    # Regression MMD on Counts

    # FILTER
    cc <- Counts[,group1.ids]
    cc <- cc[cc>0]
    quants <- quantile(cc,seq(0,1,0.02))

    mm <- MMD[,group1.ids]
    mm <- mm[!is.na(mm)]
    quants.m <- quantile(mm,seq(0,1,0.1))

    M <- rowMeans(as.matrix(MMD[,group1.ids]),na.rm=TRUE)
    idx = which(rowMin(as.matrix(Counts[,group1.ids])) > quants[2] &
                  rowMax(as.matrix(Counts[,group1.ids])) < quants[length(quants)-1] &
                  M < quants.m[length(quants.m)-1])

    C = as.matrix(Counts[idx,group1.ids])
    M = as.matrix(MMD[idx,group1.ids])


    data = data.frame(N = rowMeans(C),
                      D = log(rowMeans(M)))


    fit.Group1 <- locfit(D ~ N, data=data)
    #disp <- (exp(data$D) - exp(predict(fit.Group1,data$N)))^2
    #data <- cbind(data,disp = log(disp))
    disp <- abs(data$D - predict(fit.Group1,data$N))
    #fit.Group1 <- lm(D ~ 1/N, data=data)
    #disp <- abs(data$D - predict.lm(fit.Group1,data))

    data <- cbind(data,disp = disp)
    #data <- data[disp>0,]
    fit.disp.Group1 <- locfit(disp ~ N, data=data)
    if (plot.on){
      smoothScatter(data$N,data$D,main='group1')#,ylim=c(-10,10))
      NN <- seq(0,10^5,length.out = 10^6)
      yy <- predict(fit.Group1,log10(NN+1))
      yyd <- predict(fit.disp.Group1,log10(NN+1))
      points(log10(NN+1),yy,col='black',t='l',lwd=4)
      points(log10(NN+1),yy+yyd,col='grey',t='l',lwd=4)
      points(log10(NN+1),yy-yyd,col='grey',t='l',lwd=4)

      smoothScatter(data$N,exp(data$D),main='group1',ylim=c(0,0.5))#,ylim=c(-10,10))
      NN <- seq(0,10^5,length.out = 10^6)
      yy <- predict(fit.Group1,log10(NN+1))
      yyd <- predict(fit.disp.Group1,log10(NN+1))
      points(log10(NN+1),exp(yy),col='black',t='l',lwd=4)
      points(log10(NN+1),exp(yy+yyd),col='grey',t='l',lwd=4)
      points(log10(NN+1),exp(yy-yyd),col='grey',t='l',lwd=4)
    }
    #fit.disp.Group1 <- lm(disp ~ 1/(10^N), data=data)
    #yy=predict.lm(fit.disp.Group1,data)
    ## find outlier
    group1 = detpvals(MMD,Counts,fit.Group1,fit.disp.Group1,group1.ids)
    Group1 = group1$combined

    between.gr1 <- detpvals(MMD,Counts,fit.Group1,fit.disp.Group1,bw.group.ids)

  }  else {
    fit.Group1 <- NULL
    fit.disp.Group1 <- NULL
    between.gr1 <- NULL
    group1 = NULL
    Group1 = matrix(nrow=nrow(Counts),ncol=0)
  }

  if (!is.null(group2.ids)){

    # FILTER
    cc <- Counts[,group2.ids]
    cc <- cc[cc>0]
    quants <- quantile(cc,seq(0,1,0.02))

    mm <- MMD[,group2.ids]
    mm <- mm[!is.na(mm)]

    quants.m <- quantile(mm,seq(0,1,0.1))

    M <- rowMeans(as.matrix(MMD[,group2.ids]),na.rm=TRUE)
    idx = which(rowMin(as.matrix(Counts[,group2.ids])) > quants[2] &
                  rowMax(as.matrix(Counts[,group2.ids])) < quants[length(quants)-1] &
                  M < quants.m[length(quants.m)-1])

    C = as.matrix(Counts[idx,group2.ids])
    M = as.matrix(MMD[idx,group2.ids])


    data = data.frame(N = rowMeans(C),
                      D = log(rowMeans(M)))

    fit.Group2 <- locfit(D ~ N, data=data)
    # disp <- (exp(data$D) - exp(predict(fit.Group2,data$N)))^2
    # fit.Group2 <- lm(D ~ 1/N, data=data)
    disp <- abs(data$D - predict(fit.Group2,  data$N)) #### use predict instead of predict.lm
    data <- cbind(data,disp = disp)
    # data <- data[disp>0,]
    fit.disp.Group2 <- locfit(disp ~ N, data=data)
    if (plot.on){
      smoothScatter(data$N,data$D,main='group2')#,ylim=c(-10,10))
      NN <- seq(0,10^5,length.out = 10^6)
      yy <- predict(fit.Group2,log10(NN+1))
      yyd <- predict(fit.disp.Group2,log10(NN+1))
      points(log10(NN+1),yy,col='black',t='l',lwd=4)
      points(log10(NN+1),yy+yyd,col='grey',t='l',lwd=4)
      points(log10(NN+1),yy-yyd,col='grey',t='l',lwd=4)
    }
    group2 = detpvals(MMD,Counts,fit.Group2,fit.disp.Group2,group2.ids)
    Group2 = group2$combined
    between.gr2 <- detpvals(MMD,Counts,fit.Group2,fit.disp.Group2,bw.group.ids)

    if (!is.null(group1.ids)){
      comb <- cbind(between.gr1$combined,between.gr2$combined)
      comb[is.na(comb)] <- 10
      comb <- rowMax(comb)
    } else {
      comb <- between.gr2$combined
      comb[is.na(comb)] <- 10
    }

  } else {
    fit.Group2 <- NULL
    fit.disp.Group2 <- NULL
    group2 = NULL
    Group2 <- matrix(nrow=nrow(Counts),ncol=0)
    between.gr2 <- NULL
    comb <- between.gr1$combined
    comb[is.na(comb)] <- 10
  }
  ##   p-vals


  fdr <- comb
  fdr[fdr!=10] <- p.adjust(fdr[fdr!=10] ,method="BH")
  mmd.Table[,'pval'] <- comb
  mmd.Table[,'padj'] <- fdr

  withinGroup.pvals = cbind(p.adjust(Group1,method="BH"),
                            p.adjust(Group2,method="BH"))
  withinGroup.pvals[is.na(withinGroup.pvals)]=10
  w.p  = rowMin(withinGroup.pvals)
  pvals = mmd.Table[,'pval']
  outlier =  rep(FALSE,length(w.p))
  outlier[pvals>w.p & w.p<0.1]=TRUE


  #I=sort.int(mmd.Table[,'padj'] ,
  # na.last = NA, decreasing = FALSE, index.return = TRUE)
  # mmd.Table = mmd.Table[I$ix,]

  P=list(fit.Group1=fit.Group1, fit.Group2=fit.Group2,
         fit.disp.Group1=fit.disp.Group1, fit.disp.Group2=fit.disp.Group2,
         mmd.Table=mmd.Table,
         group1=group1,
         group2=group2,
         between.gr1=between.gr1,
         between.gr2=between.gr2,
         outlier=outlier)
  return(P)
}

detpvals <- function(M,C,fit,fit.disp,which.comps){

  p <- matrix(NA,dim(M)[1],length(which.comps))
  s <- summary(fit)
  #N <- length(s$residuals)
  #K <- s$df[2]
  N <- summary(fit)$n
  K <- fit$dp['df2']
  #RSDR <- quantile(abs(residuals(fit)),.6827)*N/(N-K)
  for (i in 1:length(which.comps)) {
    pred.y <- predict(fit,C[,which.comps[i]])
    #pred.y <- predict.lm(fit,data.frame(N=C[,which.comps[i]]))
    Dist <- log(M[,which.comps[i]])-pred.y
    DR <- predict(fit.disp,C[,which.comps[i]])
    #t1 <- Dist/RSDR
    idx <- which(Dist>0)
    t2 <- sign(Dist)*Dist^2/sqrt(DR)
    #p[,i] <- pchisq(t2,N-K,lower.tail = FALSE)
    p[,i] <- pt(t2,df=N-K,lower.tail = FALSE)
  }
  W <- -2*rowSums(log(p))
  pp <- 1-pchisq(W,2*length(which.comps))
  #Dist = rowMeans(Dist)
  #t <- Dist/RSDR
  #p <- pt(t,df=N-K,lower.tail = FALSE)
  P <- list(uncombined=p,combined=pp)
  return(P)
}
