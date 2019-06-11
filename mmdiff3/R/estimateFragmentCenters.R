#' estimate center of fragments
#'
#' This function computes average shifts between
#' forward and reverse strand, which is equivalent to the median fragment lengths.
#' In case of single-end reads, it is applied it to estimate fragment centers.
#' For paired-end reads, this step can be omitted, but it can be used to determine
#' the median fragment length.
#'
#' @inheritParams getPeakReads
#' @param shift can be set if the offset between forward and reverse strand is
#' known (e.g. 1/2 median fragment size). In this case shift will not be
#' estimated  (DEFAULT: NULL)
#' @param draw.on plot scatterplots for counts on
#' forward vs reverse strand and histograms of determined shifts
#' (DEFAULT: FALSE)
#' @param bin.length determines resolution with which fragment centers can be
#' estimated. (DEFAULT: 10)
#' @param rm.Ends if TRUE only center of fragments will be stored,
#' Position of fragment ends will be removed. (DEAFULT: TRUE)
#' @return DBAmmd object with updated slots \code{Reads} and \code{MetaData}.
#'
#' @examples
#'
#' ## Example using a small data set provided with this package
#'
#' data("MMD")
#' MMD.1 <- estimateFragmentCenters(MMD)
#'
#' # To access centers of fragments:
#' Reads.C <- Reads(MMD.1,'Center')
#'
#' # To access the determined shifts for each sample:
#' meta <- metaData(MMD.1)
#' meta$AnaData$Shifts
#'
#'
#' @seealso \code{\link{DBAmmd}}, \code{\link{getPeakReads}},
#' \code{\link{compHists}}
#'
#' @import graphics
#' @importFrom stats ccf median
#' @export

#
#
# subfunctions:
#         - getPeakIds4shift
#         - getShift
#         - shiftReads
#
#
# Gabriele Schweikert
# January 2014
# needs paired end!! checking cluster usage, speed-ups + memory improvement


estimateFragmentCenters <- function(MD, shift=NULL, draw.on=FALSE, bin.length=10,
                                    rm.Ends=TRUE, strand.specific=FALSE){

  message('checking parameters...')

  if (missing(MD))
    stop("MD object")

  Meta <- metaData(MD)
  PeakBoundary <- Meta$AnaData$PeakBoundary
  Peaks <- Regions(MD)

  Left <- Reads(MD,'Left')
  Right <- Reads(MD, 'Right')
  Center.p <- Reads(MD,'Center.p')
  Center.n <- Reads(MD, 'Center.n')
  Center <- Reads(MD, 'Center')

  numSamples <-  length(Left)
  SampleIDs <- names(Left)


  Shifts <- rep(0, numSamples)
  names(Shifts) <- SampleIDs
  for (i in 1:numSamples){

    Idx <- getPeakIds4shift(Counts(MD,'p')[,i],Counts(MD,'n')[,i], draw.on,
                            name = SampleIDs[i],strand.specific=strand.specific)
    if (is.null(shift)){
      if (length(Idx)==0){
        message('can\'t determine shift between forward and
                reverse strand. Not enough peaks!')
      } else{
        PeakLengths <- width(Peaks)[Idx] +PeakBoundary*2
        Shifts[i] <- getShift(Left[[i]][Idx],Right[[i]][Idx],
                              PeakLengths, bin.length=bin.length,
                              draw.on = draw.on,
                              name=SampleIDs[i])
      }
    } else  Shifts[i] <- shift



    if (!Meta$AnaData$pairedEnd[i]){
      R <- list(Left.p=Left[[i]],Right.n = Right[[i]])
      R <- shiftReads(R, Shifts[i])

      Center.p[[i]] <- R$Center.p
      Center.n[[i]] <- R$Center.n
      Center[[i]] <- R$Center
    }
    message("")
  } # end loop over files


  ## prepare output

  if (rm.Ends & !strand.specific){
    R <- list(Center = Center)
  } else if (rm.Ends & strand.specific){
    R <- list(Center = Center,
              Center.p = Center.p,
              Center.n = Center.n)
  } else if (!rm.Ends & !strand.specific){
    R <- list(Left = Reads(MD,'Left'),
              Right = Reads(MD,'Right'),
              Center = Center)
  } else  {
    R <- list(Left.p = Reads(MD,'Left.p'),
              Left.n = Reads(MD,'Left.n'),
              Left = Reads(MD,'Left'),
              Right.p = Reads(MD,'Right.p'),
              Right.n = Reads(MD,'Right.n'),
              Right = Reads(MD,'Right'),
              Center.p = Center.p,
              Center.n = Center.n,
              Center = Center)
  }

  Meta$AnaData$Shifts <- Shifts

  MD@MetaData <- Meta
  MD@Reads <- R
  return(MD)
}

######################
# getPeakIds4shifts: Determines which Peaks to use to for strand shift
# correction. Peaks are selected wich have total counts on positve and
# negative strand in the 9th 10-quantiles of all total counts.
#
# INPUT   - Counts.p  (nb of reads mapping to peak on + strand)
#         - Counts.n  (nb of reads mapping to peak on + strand)
#         - draw.on=TRUE
#         - name (for figure caption)
#
# OUTPUT  - Idx: index of peaks
#
#
# Gabriele Schweikert
# January 2014

getPeakIds4shift <- function(Counts.p,Counts.n,draw.on=TRUE,name=NULL,strand.specific){

  summary(Counts.p)
  summary(Counts.n)

  x <- quantile(Counts.p, probs=seq(0,1,0.1))
  ii1 <- which(Counts.p>x[9] & Counts.p<x[10])
  x <- quantile(Counts.n, probs=seq(0,1,0.1))
  ii2 <- which(Counts.n>x[9] & Counts.n<x[10])
  if (!strand.specific){
    Idx <- intersect(ii2, ii1)
  } else{
    Idx <- unique(c(ii2, ii1))
  }
  
  message('Using ',length(Idx), ' peaks to determine median fragment length')
  if (draw.on){
    if (is.null(name)){
      main = ''
    } else{
      main = name
    }
    dev.new()
    ## fn <- paste(pics.dir,'forwardVSReverse_',main,'pdf',sep='')
    ## pdf(fn)
    smoothScatter(Counts.p, Counts.n, xlab='total counts on + strand',
                  ylab='total counts on - strand', main=main)
    if (length(Idx)>0){
      points(Counts.p[Idx], Counts.n[Idx], col='red')
    }
    legend("topleft", pch=c(21),col='red', "peaks to determine strand shift")
    # dev.off()
  }
  return(Idx)
}



######################
# getshift: Determines strand shift between reads mapping to positive
# and negative strand
#
# INPUT   - Reads: Reads$Left.p
#                  Reads$Right.n
#         - Idx: index of peaks to use
#         - bin.length = 10
#         - draw.on = TRUE
#
# OUTPUT  - shift
#
#
# Gabriele Schweikert
# January 2014

getShift <- function(p,n, PeakLengths, bin.length=10, draw.on, name ){
  message('using bin.length = ', bin.length, ' to determine fragment lengths.')

  H.p <- compHist(p, PeakLengths, bin.length=bin.length)
  H.n <- compHist(n, PeakLengths, bin.length=bin.length)
  # H.p <- compHist(p, bin.length=bin.length)
  # H.n <- compHist(n, bin.length=bin.length)

  shift <- rep(0, length(H.p$Counts))
  for (i in 1:length(H.p$Counts)){
    a <- ccf(H.n$Counts[[i]], H.p$Counts[[i]], lag.max=40, plot=FALSE);
    if (length(which.max(a$acf))>0){
      shift[i] <- a$lag[which.max(a$acf)]
    } else {
      shift[i] <- 0
    }
  }

  shift <- shift*bin.length
  message('lengths statistics:')
  print(summary(shift))
  if (draw.on){
    if (is.null(name)){
      main = ''
    } else{
      main = name
    }
    dev.new()
    ## fn <- paste(picss.dir,'/strandShifts_',main,'pdf',sep='')
    ## pdf(fn)
    hist(shift, 100, xlab='shift between + and - strand [bp]', main=main)
    lines(x=c(median(shift), median(shift)), y=c(0,length(shift)),
          col='red', lwd=2)
    # dev.off()

  }
  shift <- median(shift)
  message('determined fragment lengths
          (shift between positive and reverse strand) is: ', shift )
  return(shift)
}

######################
# shiftReads: shifts Reads$p by +shift/2 and Reads$n by
# -shift/2
#
# INPUT   - Reads: Reads$p
#                  Reads$n
#         - shift
#
# OUTPUT  - Reads
#
#
# Gabriele Schweikert
# August 2012

shiftReads <- function(Reads, shift){

  p <- Reads$Left.p
  #n <- Reads$Left.n

  p <- lapply(p, function(P){P <- P+shift/2})
  #n <- lapply(n, function(P){P <- P-shift/2})

  Reads$Center.p <- p
  #Reads$Left.n <- n

  #p <- Reads$Right.p
  n <- Reads$Right.n

  #p <- lapply(p, function(P){P <- P+shift/2})
  n <- lapply(n, function(P){P <- P-shift/2})

  #Reads$Right.p <- p
  Reads$Center.n <- n

  Reads$Center <- mapply(c, p,n, SIMPLIFY=FALSE)

  return(Reads)
}


