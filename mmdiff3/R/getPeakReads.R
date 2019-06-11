#' Get reads from indexed bam files for defined regions
#'
#' This function collects all short reads from bam files that map
#' to pre-defined regions of interest. Note, that it fetches the exact start and
#' end positions of mapped fragments, not the coverage. In the case of
#' single-end reads, the left most postions of fragments mapping to the positive
#' strands and the right most positions of fragments
#' mapping to the negative strands are stored. To find centers of fragments use
#' \code{estimateFragmentCenters()}. Positions are given relative to the start
#' of the peak. Also computed are TotalCounts, i.e. number of fragments mapping
#' to a peak region, as well as number of fragments mapping to forward and
#' reverse strand.
#'
#' @param MD DBAmmd Object. This Object can be created using \code{DBAmmd()}.
#' @param PeakBoundary Defines flanking regions of peaks.
#' The estimated fragment length is a good guess
#' (DEFAULT: 200).
#' @param strand.specific if strand specific reads should be kept apart.
#' (DEFAULT:FALSE)
#' @param filter.duplicates wether duplicate fragments should be removed.
#' (DEFAULT: FALSE)
#' @param run.parallel whether to run in parallel
#' (currently no parallelization implemented)
#' (DEFAULT: FALSE)
#' @param BPPARAM an optional parameter object passed internally
#' to \code{\link{bplapply}} when \code{parallel=TRUE}.
#' If not specified, the parameters last registered with
#' \code{\link{register}} will be used.
#'
#' @return DBAmmd object with updated slots
#'
#' @seealso \code{\link{DBAmmd}}, \code{\link{estimateFragmentCenters}}
#'
#' @examples
#'
#' ## Example using a small data set provided in the MMDiffBamSubset package
#'
#' # setting the Experiment meta data:
#' ExpData <- list(dataDir=system.file("extdata", package="MMDiffBamSubset"),
#'            sampleSheet="Cfp1.csv")
#'
#' MetaData <- list('ExpData' = ExpData)
#'
#' # Creating a DBAmmd data set:
#' MMD <- DBAmmd(MetaData)
#'
#' # defining a small Region for which to get reads:
#' Regions <- GRanges(seqnames=c('chr1'),
#'            IRanges(start = c(4560912,4677889), end = c(4562680,4679681)))
#' MMD <- setRegions(MMD,Regions)
#' MMD <- getPeakReads(MMD)
#'
#' # To access Left ends of fragments mapping to positive strand:
#' Reads.L <- Reads(MMD,'Left.p')
#'
#' # To access Right ends of fragments mapping to negative strand:
#' Reads.R <- Reads(MMD,'Right.n')
#'
#' # To access Matrix of TotalCounts:
#' C.t <- Counts(MMD,whichCounts='T')
#'
#' # Counts on positive strand:
#' C.p <- Counts(MMD,whichCounts='p')
#'
#' # Counts on negative strand:
#' C.n <- Counts(MMD,whichCounts='n')
#'
#' @import GenomicRanges Rsamtools grDevices BiocParallel
#' @export
#
# subfunctions:
#         - getReadsWrapper
#         - getReads
#         - sortReads
#
#
# Gabriele Schweikert
# March 2016
# needs checking for paired-end and cluster usage

getPeakReads <- function(MD, PeakBoundary = 200, strand.specific = FALSE,
                         filter.duplicates = FALSE,
                         run.parallel = FALSE, BPPARAM=bpparam()){

  message('checking parameters...')
  if (missing(MD))
    stop("DBAmmd object required")

  message('using Peak boundary:',PeakBoundary,"bp \n")
  Meta <- metaData(MD)
  wd <- setwd(metaData(MD)$ExpData$dataDir)
  Peaks <- Regions(MD)
  numPeaks <- length(Peaks)
  samples <- Samples(MD)

  #-----------
  # Get all bam.files including input control
  #-----------
  numChips <- nrow(samples)
  chips <- as.matrix(samples$bamReads,ncol=1)
  rownames(chips) = samples$SampleID
  chips <- unique(chips)
  Todo <- chips
  if (!is.null(samples$bamControl)){
    inputs <- as.matrix(samples$bamControl,ncol=1)
    rownames(inputs) <- samples$ControlID
    if (is.null(samples$ControlID)){
      rownames(inputs) = paste("Input",1:length(inputs),sep='')
    }
    inputs <- unique(inputs)
    if (length(inputs)>0){
      Todo <- rbind(chips, inputs)
    }
  }
  todo <- as.vector(Todo)
  names(todo) <- rownames(Todo)
  #check sample and control ids
  if (any(is.null(names(todo)))) {
    print(todo)
    stop('missing SampleIDs or ControlIDs')}
  if (anyDuplicated(names(todo))>0){
    print(todo)
    stop('SampleIDs, ControlIDs not unique')
  }

  #-----------
  #Check if all files can be found
  #-----------

  if (!all(file.exists(todo))){
    idx <- which(file.exists(todo)==FALSE)
    message('Missing files:')
    message(todo[[idx]])
    stop('Check sample sheet')
  }

  #-----------
  # CHECK IF FILES ARE PAIRED-END
  #-----------
  if (is.element("PairedEnd",colnames(samples))){
    pairedEnd <- samples[,"PairedEnd"]
    pairedEnd  <- ifelse( pairedEnd  %in% c(0, "FALSE"), FALSE, TRUE)
    if  (length(pairedEnd)<length(todo)){
      pairedEnd <- c(pairedEnd,rep(NA,length(todo)-length(pairedEnd)))
    }
  } else {
    pairedEnd <- rep(FALSE,length(todo))
  }
  for (i in 1:length(todo)){
    if (!is.element("PairedEnd",colnames(samples)) |
        is.na(pairedEnd[i])){
      suppressMessages(
      pairedEnd[i] <- testPairedEndBam(todo[i])
      )
    }
  }

  #-----------
  # Initialization
  #-----------
  RawTotalCounts <- matrix(0,nrow=numPeaks,ncol=length(todo))
  colnames(RawTotalCounts) <- names(todo)
  class(RawTotalCounts) <- "integer"
  RawCounts.p <- RawTotalCounts
  RawCounts.n <- RawTotalCounts
  reads <- vector("list", length(todo))


  #-----------
  # Prepare Data for parallel use
  #-----------
  Data <- lapply(seq_len(length(todo)),
                 function(fn.id) list('pairedEnd'= pairedEnd[fn.id],
                                      'bam.file'= todo[[fn.id]]))

  if (run.parallel){
    message('starting parallel computing...')
    print(BPPARAM)

    Reads <- bplapply(Data, getReadsWrapper,BPPARAM=BPPARAM,
                      Peaks,PeakBoundary,filter.duplicates,verbose=0)

  } else{
    Reads <- vector("list", length(todo))
    for (i in 1:length(todo)){
      message('loading reads from file: ', todo[i])
      if (pairedEnd[i]){
        message('file contains paired-end reads.')
      } else{
        message('file contains single-end reads.')
      }

      Reads[[i]] <- getReadsWrapper(Data[[i]],Peaks,PeakBoundary,
                                    filter.duplicates,
                                    verbose=1)}
  }


  #---------
  # Resorting output and combining strands
  #----------

  Names <- paste(seqnames(Peaks), ':', start(Peaks), '-', end(Peaks), sep='')
  IDX <- match(Names,names(Reads[[1]]$Left.p))

  Left.p <- vector("list", length(todo))
  names(Left.p) <- names(todo)
  Left.n <- Left.p
  Left <- Left.p
  Right.n <- Left.p
  Right.p <- Left.p
  Right <- Left.p
  Center.p <- Left.p
  Center.n <- Left.p
  Center <- Left.p

  for (i in 1:length(todo)){
    Left.p[[i]] <- Reads[[i]]$Left.p[IDX]
    Right.n[[i]] <- Reads[[i]]$Right.n[IDX]

    if (pairedEnd[i]){
      Left.n[[i]] <- Reads[[i]]$Left.n[IDX]
      Right.p[[i]] <- Reads[[i]]$Right.p[IDX]
      Center.p[[i]] <- Reads[[i]]$Center.p[IDX]
      Center.n[[i]] <- Reads[[i]]$Center.n[IDX]

      Right[[i]] <- mapply(c, Right.p[[i]], Right.n[[i]], SIMPLIFY=FALSE)
      Left[[i]] <- mapply(c, Left.p[[i]], Left.n[[i]], SIMPLIFY=FALSE)
      Center[[i]] <- mapply(c, Center.p[[i]], Center.n[[i]], SIMPLIFY=FALSE)
    } else{
      Left[[i]] <- Left.p[[i]]
      Right[[i]] <- Right.n[[i]]
    }

    RawTotalCounts[,i] <- Reads[[i]]$RawTotalCounts[IDX]
    RawCounts.p[,i] <- Reads[[i]]$RawCounts.p[IDX]
    RawCounts.n[,i] <- Reads[[i]]$RawCounts.n[IDX]
    names <- names(Left.p[[i]])
    rNames1 <- names(Reads[[i]]$RawTotalCounts)[IDX]
    rNames2 <- names(Reads[[i]]$RawCounts.p)[IDX]
    rNames3 <- names(Reads[[i]]$RawCounts.n)[IDX]
    if (!identical(Names ,names) | !identical(Names ,rNames1) |  !identical(rNames1 ,rNames2) | !identical(rNames1 ,rNames3)){
      stop('something wrong')
    }
    if (is.null(rownames(RawTotalCounts) )){
      rownames(RawTotalCounts) <- rNames1
      rownames(RawCounts.p) <- rNames1
      rownames(RawCounts.n) <- rNames1
    } else {
      if (!identical(rNames1 ,rownames(RawTotalCounts)) ){
        stop('something wrong')
      }
    }
    class(RawTotalCounts) <- "integer"
    class(RawCounts.p) <- "integer"
    class(RawCounts.n) <- "integer"
  }


  if (strand.specific){
    Reads <- list(Left.p = Left.p,
                  Left.n = Left.n,
                  Left = Left,
                  Right.p = Right.p,
                  Right.n = Right.n,
                  Right = Right,
                  Center.p = Center.p,
                  Center.n = Center.n,
                  Center = Center)
  } else{
    Reads <- list(Left = Left,
                  Right = Right,
                  Center = Center)
  }

  Meta <- metaData(MD)
  Meta$AnaData$pairedEnd <- pairedEnd
  Meta$AnaData$PeakBoundary <- PeakBoundary
  Meta$AnaData$filter.duplicates <- filter.duplicates
  MD@MetaData <- Meta
  MD@Reads <- Reads
  MD@RawTotalCounts = RawTotalCounts
  MD@RawCounts.p=RawCounts.p
  MD@RawCounts.n=RawCounts.n

  setwd(wd)
  return(MD)
}



######################
# getReads: This function uses Rsamtools to collect all short reads on
# chromosome chrom in the bam.file that match to regions defined by
# Peaks. 5' Positions of reads are returned in Reads$p and
# Reads$n for reads mapping to positive or negative strand,
# respectively. Also computes total number of reads mapping to a given peak.
#
#
# INPUT   - Peaks: GRanges
#         - chrom: chromosome id
#         - bam.file
#         - pairedEnd
#
# OUTPUT  - Reads$Reads  for single-end reads contains 'Left.p' (Left end of fragments mapping to positive strand)
#                        and 'Right.n' (Right end of fragments mapping to negative strand)
#
#                        for paired-end reads additionally 'Right.p','Left.n', 'Center.p','Center.n',
#
#         - Reads$RawTotalCounts
#         - Reads$RawCounts.p
#         - Reads$RawCounts.n
#
# Gabriele Schweikert
# January 2014
getReads <- function(Data, bam.file, pairedEnd, PeakBoundary,filter.duplicates){
  Peaks <- Data$Peaks
  chrom <- Data$chrom
  Names <- paste(seqnames(Peaks), ':', start(Peaks), '-', end(Peaks), sep='')
  header <- scanBamHeader(bam.file)
  if (filter.duplicates){
    isDuplicate = FALSE
  } else{
    isDuplicate = NA
  }

  if (length(which(names(header[[1]]$targets)==chrom))==0) {
    chrchrom <- paste('chr', chrom, sep='')
    rmchrchrom <- unlist(strsplit(chrom, 'chr'))[2]
    if (length(which(names(header[[1]]$targets)==chrchrom))>0) {
      message('changing ', chrom, ' to ', chrchrom)
      chrom <- chrchrom
    } else if (length(which(names(header[[1]]$targets)==rmchrchrom))>0) {
      message('changing ', chrom, ' to ', rmchrchrom)
      chrom <- rmchrchrom
    } else {
      message('no reads found mapping to ', chrom)
      Left.p <- vector('list', length(Peaks))
      names(Left.p) <- Names
      Left.n <- Left.p
      Right.p <- Left.p
      Right.n <- Left.p
      Center.p <- Left.p
      Center.n <- Left.p
      Reads <- list(Left.p, Left.n,Right.p, Right.n,Center.p,Center.n)
      names(Reads) <- c('Left.p', 'Left.n','Right.p', 'Right.n','Center.p','Center.n')

      RawCounts.p <- rep(0,length(Peaks))
      RawCounts.n <- rep(0,length(Peaks))
      RawTotalCounts <- rep(0,length(Peaks))
      Reads <-  list(Reads=Reads,
                     RawTotalCounts=RawTotalCounts,
                     RawCounts.p=RawCounts.p,
                     RawCounts.n=RawCounts.n)
      return(Reads)
    }
  }

  which <- IRangesList(IRanges(start(Peaks), end(Peaks)), compress=FALSE)
  names(which) <- c(chrom)
  if (pairedEnd==TRUE){
    what <- c("qname","strand", "pos","isize")
    param <- ScanBamParam(which=which, what=what,
                          scanBamFlag(isProperPair = TRUE,
                                      isSecondMateRead = FALSE,
                                      isMinusStrand = FALSE,
                                      isDuplicate = isDuplicate))
    bam <- scanBam(bam.file, param=param)
    if (!all(bam$strand=='+')){
      warning("something wrong with bam filter")
    }

    Left.p <- lapply(bam, function(bam){bam$pos})
    Right.p <- lapply(bam, function(bam){bam$pos + bam$isize})
    Center.p <- lapply(bam, function(bam){bam$pos + .5*bam$isize})

    ## Negative Strand
    param <- ScanBamParam(which=which, what=what,
                          scanBamFlag(isProperPair = TRUE,
                                      isSecondMateRead = TRUE,
                                      isMinusStrand = FALSE,
                                      isDuplicate = isDuplicate))
    bam <- scanBam(bam.file, param=param)
    if (!all(bam$strand=='+')){
      warning("something wrong with bam filter")
    }
    Left.n <- lapply(bam, function(bam){bam$pos})
    Right.n <- lapply(bam, function(bam){bam$pos + bam$isize })
    Center.n <- lapply(bam, function(bam){bam$pos + .5*bam$isize})
    # get total coverage of Peaks
    Counts.Left.n <-  mapply(Left.n, FUN=function(pos){length(pos)})
    Counts.Right.p <-  mapply(Right.p, FUN=function(pos){length(pos)})

  } else {
    what <- c("rname", "strand", "pos","seq")
    param <- ScanBamParam(which=which, what=what,
                          scanBamFlag(isDuplicate = isDuplicate))
    bam <- scanBam(bam.file, param=param)
    Left.p <- lapply(bam, function(bam){bam$pos[bam$strand=='+']})
    read.Length <- unique(width(bam[[1]]$seq))
    if (length(read.Length)>1){
      warning('reads have different read lengths:')
      print(read.Length)
      read.Length = mean(read.Length)
    }
    Right.n <- lapply(bam, function(bam){bam$pos[bam$strand=='-']+read.Length})
    Left.n = NULL
    Right.p = NULL
    Center.p = NULL
    Center.n = NULL
  } #single-end

  ## Read positions are shifted to peak coordinates
  ## i.e. for each peak reads mapping exactly to the 5'end of the peak region
  ## are at pos = PeakBoundary+1
  starts <- start(Peaks)
  Left.p <- lapply(seq_len(length(Peaks)),
                   function(row) Left.p[[row]] - starts[row]+PeakBoundary+1)
  Right.n <- lapply(seq_len(length(Peaks)),
                    function(row) Right.n[[row]] - starts[row]+PeakBoundary+1)
  if ( pairedEnd==TRUE){
    Left.n <- lapply(seq_len(length(Peaks)),
                     function(row) Left.n[[row]] - starts[row]+PeakBoundary+1)
    Right.p <- lapply(seq_len(length(Peaks)),
                      function(row) Right.p[[row]] - starts[row]+PeakBoundary+1)
    Center.n <- lapply(seq_len(length(Peaks)),
                       function(row) Center.n[[row]] - starts[row]+PeakBoundary+1)
    Center.p <- lapply(seq_len(length(Peaks)),
                       function(row) Center.p[[row]] - starts[row]+PeakBoundary+1)
    names(Left.n) <- Names
    names(Right.p) <- Names
    names(Center.p) <- Names
    names(Center.n) <- Names

  }
  ## get total coverage of Peaks
  RawCounts.p <-  mapply(Left.p, FUN=function(pos){length(pos)})
  RawCounts.n <-  mapply(Right.n, FUN=function(pos){length(pos)})
  RawTotalCounts <- RawCounts.p+ RawCounts.n
  if ( pairedEnd==TRUE &&
       ( sum(RawCounts.n-Counts.Left.n )!=0  |
         sum(RawCounts.n-Counts.Left.n)!=0 ))
    stop('')
  names(Left.p) <- Names
  names(Right.n) <- Names
  names(RawTotalCounts ) <- Names
  names(RawCounts.p) <- Names
  names(RawCounts.n) <- Names

  Reads <-  list(Left.p = Left.p,
                 Left.n = Left.n,
                 Right.p = Right.p,
                 Right.n = Right.n,
                 Center.p = Center.p,
                 Center.n = Center.n,
                 RawTotalCounts=RawTotalCounts,
                 RawCounts.p=RawCounts.p,
                 RawCounts.n=RawCounts.n)


  return(Reads)
}

getReadsWrapper <- function(Data,Peaks,PeakBoundary,filter.duplicates,verbose=1){

  bam.file <- Data$bam.file
  pairedEnd <- Data$pairedEnd
  chroms <- unique(as.character(seqnames(Peaks)))

  Left.p <- vector("list", length(Peaks))
  Left.n <- Left.p
  Right.p <- Left.p
  Right.n <- Left.p
  Center.p <- Left.p
  Center.n <- Left.p
  RawTotalCounts <- rep(0,length(Peaks))
  RawCounts.p <- RawTotalCounts
  RawCounts.n <- RawTotalCounts
  peak.id <- 0
  if (length(chroms)==0){
    return()
  }

  for (j in 1:length(chroms)){
    if (is.element(chroms[j],seqnames(seqinfo(Peaks)))){
      P <- Peaks[seqnames(Peaks)==chroms[j]]
      Data = list('Peaks' = P,'chrom' = chroms[j])
      if (verbose>0){
        message('-chromosome: ', chroms[j], ' containing ', length(P), ' peaks')
      }
      reads <- getReads(Data, bam.file = bam.file,
                        pairedEnd = pairedEnd,
                        PeakBoundary = PeakBoundary,filter.duplicates)

      ids <- peak.id+seq(1, length(P))

      if (length(ids) != length(reads$Left.p)){
        message('trouble on chromosome', chroms[j]) }

      Left.p[ids] <- reads$Left.p
      Right.n[ids] <- reads$Right.n
      names(Left.p)[ids] <- names(reads$Left.p)
      names(Right.n)[ids] <- names(reads$Right.n)
      RawTotalCounts[ids] <- reads$RawTotalCounts
      RawCounts.p[ids] <- reads$RawCounts.p
      RawCounts.n[ids] <- reads$RawCounts.n

      names(RawTotalCounts)[ids] <- names(reads$RawTotalCounts)
      names(RawCounts.p)[ids] <- names(reads$RawCounts.p)
      names(RawCounts.n)[ids] <- names(reads$RawCounts.n)

      if (pairedEnd){
        Left.n[ids] <- reads$Left.n
        Right.p[ids] <- reads$Right.p
        Center.n[ids] <- reads$Center.n
        Center.p[ids] <- reads$Center.p

        names(Left.n)[ids] <- names(reads$Left.n)
        names(Right.p)[ids] <- names(reads$Right.p)
        names(Center.p)[ids] <- names(reads$Center.p)
        names(Center.n)[ids] <- names(reads$Center.n)
      }

      peak.id <- ids[length(ids)]
    }

  }
  if (pairedEnd){
    Reads <-  list(Left.p = Left.p,
                   Left.n = Left.n,
                   Right.p = Right.p,
                   Right.n = Right.n,
                   Center.p = Center.p,
                   Center.n = Center.n,
                   RawTotalCounts = RawTotalCounts,
                   RawCounts.p = RawCounts.p,
                   RawCounts.n = RawCounts.n)
  }else{
    Reads <-  list(Left.p = Left.p,
                   Right.n = Right.n,
                   RawTotalCounts = RawTotalCounts,
                   RawCounts.p = RawCounts.p,
                   RawCounts.n = RawCounts.n)
  }

  return(Reads)
}


sortReadsWrapper <- function(Reads,Peaks,todo,pairedEnd,BPPARAM){
  Left.p <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                     Peaks=Peaks,
                     which='Left.p')
  names(Left.p) <- names(todo)
  Right.n <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                      Peaks=Peaks,
                      which='Right.n')
  names(Right.n) <- names(todo)

  out <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                  Peaks=Peaks,
                  which='RawTotalCounts')
  RawTotalCounts <- matrix(unlist(out),
                           ncol = length(out), byrow = FALSE)
  class(RawTotalCounts) <- "integer"
  rownames(RawTotalCounts) <- names(out[[1]])
  colnames(RawTotalCounts) <- names(todo)

  out <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                  Peaks=Peaks,
                  which='RawCounts.p')
  RawCounts.p <- matrix(unlist(out), ncol = length(out), byrow = FALSE)
  class(RawCounts.p) <- "integer"
  rownames(RawCounts.p) <- names(out[[1]])
  colnames(RawCounts.p) <- names(todo)

  out <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                  Peaks=Peaks,
                  which='RawCounts.n')
  RawCounts.n <- matrix(unlist(out), ncol = length(out), byrow = FALSE)
  class(RawCounts.n) <- "integer"
  rownames(RawCounts.n) <- names(out[[1]])
  colnames(RawCounts.n) <- names(todo)

  if (any(pairedEnd)){
    Left.n <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                       Peaks=Peaks,
                       which='Left.n')

    Right.p <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                        Peaks=Peaks,
                        which='Right.p')


    Center.n <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                         Peaks=Peaks,
                         which='Center.n')

    Center.p <- bplapply(Reads, sortReads,BPPARAM=BPPARAM,
                         Peaks=Peaks,
                         which='Center.p')


  } else {
    Left.n <- vector(mode="list",length=length(todo))
    Right.p <- Left.n
    Center.n <- Left.n
    Center.p <- Left.n
  }
  names(Left.n) <- names(todo)
  names(Right.p) <- names(todo)
  names(Center.n) <- names(todo)
  names(Center.p) <- names(todo)
  Reads <- list(Left.p = Left.p,
                Left.n = Left.n,
                Right.p = Right.p,
                Right.n = Right.n,
                Center.p = Center.p,
                Center.n = Center.n)
  OUT=list(Reads=Reads,
           RawCounts.n=RawCounts.n,
           RawCounts.p=RawCounts.p,
           RawTotalCounts=RawTotalCounts)
  return(OUT)
}


sortReads <- function(r, Peaks, which){
  numPeaks <- length(Peaks)
  pos <- vector("list", numPeaks)

  peak.id <- 0
  for (chr in 1:length(r)){
    if (length(r[[chr]])==0){
      next
    }

    reads <- r[[chr]][[which]]

    if (is.null(reads)){
      pos=NULL
      return(pos)
    }
    ids <- peak.id+seq(1, length(reads))
    pos[ids] <- reads
    names(pos)[ids] <- names(reads)
    peak.id <- ids[length(ids)]
  }


  if (length(pos)!=length(Peaks)){
    stop('some peaks missing')
  }

  ## compare names with peak coordinates
  tmp <- lapply(names(pos),function(n){
    n=unlist(strsplit(n,':'));pos=n[2]})
  st1 <- sapply(tmp,function(pos){
    n=unlist(strsplit(pos,'-'));pos=as.integer(n[1])})
  en1 <- sapply(tmp,function(pos){
    n=unlist(strsplit(pos,'-'));pos=as.integer(n[2])})

  st2 <- start(Peaks)
  en2 <- end(Peaks)
  if (st1 !=st2 || en1!=en2){
    message("check if some peaks start < 0")
    stop("something wrong with peak sorting")

  }


  return(pos)
}




