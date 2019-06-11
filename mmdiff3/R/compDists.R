#' Compute distances between Peaks
#'
#' This function computes pairwise distances between histograms
#' according to the dist.method (MMD, KS). For large data sets it is a bit
#' time consuming.
#'
#' @inheritParams DBAmmd-Accessors
#' @inheritParams getPeakReads
#' @param sigma  sigma parameter of the RBF kernel,
#' determining the distance (along the genome) at which fragment counts
#' decorrelate. If set to NULL, 100 random Peaks are used to determine sigma
#' heuristically according to the method described in the MMDiff paper [1].
#' (DEFUALT: NULL)
#' @param CompIDs (2 x C) Matrix specifying for pairs of Sample IDs for which
#' distances should be computed. If NULL all pairwise distances will be computed.
#' (DEFAULT: NULL)
#' @return DBAmmd object with updated slot \code{Dists}
#' @seealso \code{\link{DBAmmd}}, \code{\link{plotDists}},
#' \code{\link{plotDISTS4Peak}}, \code{\link{compPvals}}
#'
#' @examples
#'
#' ## Example using a small data set provided with this package:
#'
#' data("MMD")
#' MMD.1 <- compDists(MMD)
#'
#' # To inspect the computed distances:
#' D <- Dists(MMD.1,dist.method='MMD')
#' head(D)
#'
#' # To analyse the result:
#' plotDists(MMD.1)
#'
#' @import BiocParallel
#' @importFrom stats quantile median ks.test
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#' @author Gabriele Schweikert \email{G.Schweikert@ed.ac.uk}
#' @references [1] Schweikert et al.  BMC Genomics 2013
#' ...
# subfunctions:
#         - compDist.R
#         - mmd.R
#
# check parameters
#
# Gabriele Schweikert
# March 2016




compDists <- function(MD, dist.method='MMD',sigma=NULL,CompIDs=NULL,
                      strand.specific=FALSE,
                      run.parallel = FALSE, BPPARAM=bpparam(),background_intensity=0.2){
  verbose=1
  message('checking parameters...')

  if (missing(MD))
    stop("DBAmmd object")

  Meta <- metaData(MD)

  PeakBoundary <- Meta$AnaData$PeakBoundary
  Peaks <- Regions(MD)
  
  print(length(Peaks))

  if (length(Peaks)<200){
    numPeaksForSigma <- length(Peaks)/2-1
  } else {
    numPeaksForSigma <- 100
  }
  
  print(numPeaksForSigma)
 
  if (!strand.specific){
    Pos.C <- Reads(MD,'Center')
    Counts <- Counts(MD)
  } else {
    STR <- strand(Peaks)
    Pos.C <- Reads(MD,'Center.p')
    Pos.Cn <- Reads(MD,'Center.n')
    idx <- which(as.vector(STR)=='-')
    for (i in 1:length(Pos.C)){
      Pos.C[[i]][idx] <- Pos.Cn[[i]][idx]
    }
    rm(Pos.Cn)
    
    Counts <- Counts(MD,'p')
    Counts.n <- Counts(MD,'n')
    Counts[idx,] <- Counts.n[idx,] 
    rm(Counts.n)
  }
  # flanking regions

  set.seed(1236436)
  Flanks <- vector(mode='list',length=ncol(Counts))
  names(Flanks) <- colnames(Counts)
  for (i in 1:ncol(Counts)){
    N <- median(Counts[,i]/width(Peaks)*PeakBoundary)
    Flanks[[i]] <- sample(PeakBoundary,N,replace=TRUE)
  }
  ###

  nSamples <- numSamples(MD)
  sampleIDs <- Samples(MD)$SampleID

  if (is.null(CompIDs)){

    ## ----------------
    ## establish all pairwise comparisons for given contrasts
    ## -----------------
    Contrasts <- MD@Contrasts
    if (length(Contrasts)>0){
      message('computing distances for all contrasts')
      CompIDxs <- compIDsFromContrasts(Contrasts)
      CompIDs <- rbind(sampleIDs[CompIDxs[1,]],sampleIDs[CompIDxs[2,]])
      nComps <- ncol(CompIDs)
    } else {

      ## ----------------
      ## establish all pairwise comparisons
      ## -----------------
      message('computing all distances')

      nComps <- (nSamples^2-nSamples)/2
      CompIDs <- matrix(0,2,nComps)
      CompIDxs <- matrix(0,2,nComps)
      st <- 0
      for (j in 1:(nSamples-1)){
        CompIDs[1,st+1:(nSamples-j)] <- sampleIDs[j]
        CompIDs[2,st+1:(nSamples-j)] <- sampleIDs[(j+1):nSamples]

        CompIDxs[1,st+1:(nSamples-j)] <- j
        CompIDxs[2,st+1:(nSamples-j)] <- (j+1):nSamples
        st <- st+nSamples-j
      }
    }
  } else{
    nComps <- ncol(CompIDs)
    CompIDxs <- matrix(0,2,nComps)
    for (j in 1:nComps){
      CompIDxs[1,j] = which(sampleIDs==CompIDs[1,j])
      CompIDxs[2,j] = which(sampleIDs==CompIDs[2,j])
    }
  }
  
  if(dist.method=='MMD2'){
    
    message('preparing negative contrast...')
    ### Mark Intensity
    NegativeContrast <- vector(mode='list',length=ncol(Counts))
    names(NegativeContrast) <- colnames(Counts)
    for (i in 1:ncol(Counts)){
      N <- median(Counts[,i]/width(Peaks))*width(Peaks)
      NegativeContrast[[i]] <- N*background_intensity
    }
    
    message('estimating sigma...')
    Sigma <- matrix(0,numPeaksForSigma,nSamples )
    colnames(Sigma) <- sampleIDs
    
    for (i in 1:nSamples){
      C <- Counts[,i]
      peak.ids <- which(C>quantile(C)[2]|C<=quantile(C)[4])
      peak.ids <- sample(peak.ids, numPeaksForSigma)
      R <- Pos.C[[i]][peak.ids]
      for(j in 1:numPeaksForSigma){
        Sigma[j,i] <- chooseSigma(R[[j]],R[[j]])
      }
    }
    
    summary(Sigma)
    sigma <- median(as.vector(Sigma), na.rm = TRUE)

    L <- max(width(Peaks))+2*PeakBoundary
    Meta$AnaData$NegativeContrast <- NegativeContrast
    Meta$AnaData$MMDKernelFinalSigma <- sigma
    Meta$AnaData$LUT<-buildMMDLUT(L, sigma)
    Meta$AnaData$maxval<-L
    MD@MetaData <- Meta
    
  }
  
  if (dist.method=='MMD'){
    KernelMatrix <- Meta$AnaData$KernelMatrix
    if (!is.null(KernelMatrix)){
      message('using previously computed KernelMatrix')
    } else {

      ## ----------------
      ## 1. estimate sigma
      ## ----------------
      message('estimating sigma...')

      Sigma <- matrix(0,numPeaksForSigma,nSamples )
      colnames(Sigma) <- sampleIDs

      for (i in 1:nSamples){
        C <- Counts[,i]
        peak.ids <- which(C>quantile(C)[2]|C<=quantile(C)[4])
        peak.ids <- sample(peak.ids, numPeaksForSigma)
        R <- Pos.C[[i]][peak.ids]
        for(j in 1:numPeaksForSigma){
          Sigma[j,i] <- chooseSigma(R[[j]],R[[j]])
        }
      }
      summary(Sigma)
      sigma <- median(as.vector(Sigma),na.rm = TRUE)
      print(sigma)

      ## ----------------
      ## 2. precompute Kernel matrix
      ## ----------------

      message('pre-computing Kernel matrix...')
      L <- max(width(Peaks))+2*PeakBoundary
      KernelMatrix <- compKernelMatrix(seq(1,L),sigma)

      Meta$AnaData$Flanks <- Flanks
      Meta$AnaData$MMDKernelSigma <- Sigma
      Meta$AnaData$MMDKernelFinalSigma <- sigma
      Meta$AnaData$KernelMatrix <- KernelMatrix

      MD@MetaData <- Meta
    }
  } #end if (dist.method=='MMD')



  ## ----------------
  ## 3. compute DISTs for all peaks and all sample pairs
  ## ----------------

  message(paste('computing', nComps ,'pair-wise distances...'))
  CompNames <- paste(CompIDs[1,],CompIDs[2,],sep=' vs ')
  #D <- matrix(NA,length(Peaks),nComps)


  Data <- lapply(seq_len(nComps),
                 function(i){list('i1' = CompIDs[1,i],
                                  'i2' = CompIDs[2,i])})
  MD.small <-  DBAmmd()
  MD.small@MetaData <- MD@MetaData
  MD.small@Reads <- list('Center'= Pos.C)
  MD.small@rowRanges <- MD@rowRanges

  if (run.parallel){
    message('starting parallel computing...')
    print(BPPARAM)

    D <- bplapply(Data, mmdWrapper,BPPARAM=BPPARAM,verbose=0,MD=MD.small,
                    dist.method=dist.method)


  } else {
    D <- vector(mode='list',length=nComps)
    for (i in 1:nComps){
      if (verbose>0){
        message(paste('computing distances for',
                      CompIDs[1,i],'vs',  CompIDs[2,i]))
      }



      D[[i]] <- mmdWrapper(Data[[i]],verbose=1,MD=MD.small,
                           dist.method=dist.method)
    } # end loop over comps
  }
  #D.L[is.na(D.L)] = D.R[is.na(D.L)]
  #D.R[is.na(D.R)] = D.L[is.na(D.R)]
  #DISTs <- (D.L+D.R)/2
  DISTs <- matrix(NA,length(Peaks),nComps)
  for (i in 1:nComps){
    DISTs[,i] <- D[[i]]
  }
  CompNames <- paste(CompIDs[1,],CompIDs[2,],sep=' vs ')
  colnames(DISTs) <- CompNames
  rownames(DISTs) <- names(Pos.C[[1]])


  MCounts <- matrix(NA,length(Peaks),nComps)
  for (i in 1:nComps){
    MCounts[,i] = rowMax(Counts[,c(CompIDxs[1,i],CompIDxs[2,i])])
  }
  colnames(MCounts) <- CompNames
  rownames(MCounts) <- names(Pos.C[[1]])
  if (length(MD@mCounts)==0){
    MD@mCounts <- MCounts
  } else {
    MD@mCounts <- cbind(MD@mCounts,MCounts)
  }

  ## ----------------
  ## 4. preparing Output
  ## ----------------

  message('preparing Output...')

  MD <- setDists(MD,DISTs,dist.method)

  return(MD)
}



chooseSigma <- function(x,y){
  x=as.matrix(x,nrow=(length(x)))
  y=as.matrix(y,nrow=(length(y)))
  Kxy=vectorized_pdist.sq(x,y)

  mdist <- median(Kxy[Kxy!=0])

  mdist <- sqrt(mdist/2)
  sigma <- mdist

  # X <- seq(1,max(c(x,y)))
  # D <- as.matrix(dist(X,diag=TRUE,upper=TRUE))
  # DISTs <- D[x,y]
  # Kxy <- DISTs^2
  if (is.na(sigma)){
    return(sigma)}

  if (sigma ==0 ){
    sigma <- 1}
  return(sigma)
}


compKernelMatrix <- function(x,sigma){
  ## should do the same as rbf_dot
  #DISTs <- as.matrix(dist(x),diag = TRUE, upper = TRUE)
  #K <-  exp(-1/2/sigma^2 * DISTs^2)
  x=as.matrix(x,nrow=(length(x)))
  DISTs.sq <- vectorized_pdist.sq(x,x)
  K <-  exp(-1/2/sigma^2 * DISTs.sq)
  return(K)
}

mmdWrapper <- function(Data,verbose=1,MD,dist.method) {

  i1 <- Data$i1
  i2 <- Data$i2

  Meta <- metaData(MD)
  PeakBoundary <- Meta$AnaData$PeakBoundary
  Peaks <- Regions(MD)
  Ls <- width(Peaks)+PeakBoundary

  KernelMatrix <- Meta$AnaData$KernelMatrix
  Flanks <- Meta$AnaData$Flanks
  NegativeContrast <- Meta$AnaData$NegativeContrast
  sigma <- Meta$AnaData$MMDKernelFinalSigma

  Pos.C <- Reads(MD,'Center')
  PosA <- Pos.C[[i1]]
  PosB <- Pos.C[[i2]]
  FlanksA <- Flanks[[i1]]
  FlanksB <- Flanks[[i2]]

  if (dist.method=='MMD'){
    Data <- lapply(seq_len(length(PosA)),
                 function(row) list('PosA'= c(PosA[[row]],FlanksA,Ls[row]+FlanksA),
                                    'PosB'= c(PosB[[row]],FlanksB,Ls[row]+FlanksB)))
  }
  
  if (dist.method=='MMD2'){
    Data <- lapply(seq_len(length(PosA)),
                   function(row) list('PosA'= c(PosA[[row]]),
                                      'PosB'= c(PosB[[row]])))
  }

  D <- rep(NA,length(PosA))
  for (j in seq_len(length(PosA))){
    if (dist.method=='MMD'){
      D[j] <- mmd(Data[[j]], KernelMatrix = KernelMatrix)
    } else if (dist.method=='KS'){
      stop('needs checking and fixing')
      KS <- ks.test(Data$posA,Data$posB)
      D[j] <- KS$statistic
    } else if (dist.method=='MMD2'){
      maxval <- Meta$AnaData$maxval
      lut <-Meta$AnaData$LUT
      if(NegativeContrast[[i1]][j] < 1 || NegativeContrast[[i2]][j] < 1) {
        message(paste("no contrast err",i1,"vs",i2))
      }
      
      bounds <- c(0, Ls[j]+PeakBoundary)
      D[j] <- computeDist(Data[[j]]$PosA,
                          Data[[j]]$PosB,
                          bounds, sigma,
                          0,
                          NegativeContrast[[i1]][j]/2,
                          NegativeContrast[[i2]][j]/2,
                          maxval,
                          lut
                          )
    }
  }
  return(D)
}




mmd <- function(Data,KernelMatrix) {

  posA <- Data$PosA
  posB <- Data$PosB
  m <- length(posA)
  n <- length(posB)

  if (m <2 | n< 2){
    return(NA)
  }
  Kxx <- KernelMatrix[posA,posA]
  Kyy <- KernelMatrix[posB,posB]
  Kxy <- KernelMatrix[posA,posB]


  ## Get test statistic (unbiased)
  # MMDf = ( Kxx + Kyy - Kxy - t(Kxy) )
  # MMDf = MMDf - diag(diag(MMDf))
  # testStat = 1/m/(m-1) * sum(sum( MMDf ));
  # testStat = testStat * m; #null distirbution on m*MMD

  ## MMD statistic. Here we use biased
  ## v-statistic. NOTE: this is m * MMD_b
  # testStat = 1/m * sum(sum(Kxx + Kyy - Kxy - t(Kxy)))
  biased <- sqrt(sum(Kxx)/(m*m) +  sum(Kyy)/(n * n) - 2/m/n * sum(Kxy))
  # diag(Kxx) <- 0
  # diag(Kxy) <- 0
  # unbiased <- sqrt(sum(Kxx)/(m*(m-1)) +  sum(Kyy)/(n * (n-1)) - 2/m/n * sum(Kxy))
  return(biased)
}

vectorized_pdist.sq <- function(A,B) {
  #from Alex Smola
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))

  m = nrow(A)
  n = nrow(B)

  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  d.sq = tmp - 2 * tcrossprod(A,B)
  return(d.sq)
}



