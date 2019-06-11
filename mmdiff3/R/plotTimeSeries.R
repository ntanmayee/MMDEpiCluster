#' plot Time Series
#'
#' This function plots histograms of fragment positions over a pre defined
#' regions of interests / peaks. Can also show occurences of Sequence motifs and
#' annotated objects (e.g. genes).
#'
#' @inheritParams getPeakReads
#' @inheritParams DBAmmd-Accessors
#' @param Peak.id Peak id to specify which Peak to plot.
#' (coresponding to names of Regions(MD))
#' @param NormMethod whether to apply normailzation factors.
#' currently no normalization method implemented (DEFAULT: None)
#' @param plot.input whether to plot input controls (DEFAULT: TRUE)
#' @param Motifs TF binding sites (DEFAULT: NULL)
#' @param Motifcutoff (Default: "80\%")
#' @param anno either a GRanges objects containing annotated objects, e.g. genes,
#' or a list of GRanges Objects. (Default: NULL)
#' @param xaxt (Default: NULL)
#' @param xlim (Default: NULL)
#' @param ylim (Default: NULL)
#'
#' @examples
#' dev.off()
#' data("MMD")
#' plotPeak(MMD,Peak.id='6',plot.input=FALSE)
#'
#' # add annotation (Overlapping genes)
#' data("mm9-Genes")
#' GR <- list(UCSCKnownGenes = GR)
#' plotPeak(MMD, Peak.id='6', plot.input = FALSE, anno=GR)
#'
#' # add TF binding sites
#' library('MotifDb')
#' motifs <- query(query(MotifDb, 'Mmusculus'), 'E2F')
#' plotPeak(MMD, Peak.id='6', plot.input = FALSE,
#'        Motifs=motifs,Motifcutoff="80%")
#'
#' # split peaks by contrast
#' plotPeak(MMD, Peak.id='6', plot.input = FALSE, whichContrast=1,
#'        Motifs=motifs,Motifcutoff="80%",anno=GR)
#'
#'
#' @import BSgenome GenomicRanges Biostrings  S4Vectors
#' @export
#
# TO DO check NormFactors!! from DBA
# check Input is working
# Potentially add seqlogos
# check genome in MD if motif
# make an option to plot strands independently

plotTimeSeries <- function(MD, Peak.id, whichCondition=NULL,whichTime=NULL,
                           NormMethod=NULL,
                           plot.input=FALSE,strand.specific=FALSE,whichPos="Center",
                           whichContrast=NULL,
                           Motifs=NULL,
                           Motifcutoff="80%",
                           anno=NULL,
                           xaxt=NULL,xlim=NULL,ylim=NULL){
  
  if (missing(MD))
    stop("MD is missing")
  
  if (missing(Peak.id))
    stop("No Peak specified")
  
  Peak.idx <- match(Peak.id, names(Regions(MD)))
  
  Meta <- metaData(MD)
  PeakBoundary <- Meta$AnaData$PeakBoundary
  if (strand.specific){
    R <- Regions(MD)[Peak.idx]
    if (as.vector(strand(R))=='+'){
      whichPos <- paste(whichPos,'.p',sep='')
    } else if (as.vector(strand(R))=='-') {
      whichPos <- paste(whichPos,'.n',sep='')
    } else if (as.vector(strand(R))=='*'){
      message('strand of region not specified')
    }
    
  }
  
  hists <- Hists(MD,whichPos)
  
  samples <- Samples(MD)
  
  if (is.null(hists)){
    message('whichPos = ', whichPos)
    stop("Peak histograms missing. run 'compHists', first. ")
  }
  
  Peak <- hists[Peak.idx]
  PeakName <- names(Peak)
  Peak <- Peak[[1]]
  
  if (is.null(whichCondition)){
    message("No Condition specified, plotting all Conditions")
    whichCondition <- unique(Meta$ExpData$samples$Condition)
  }
  if (is.null(whichTime)){
    message("No Time points specified, plotting all Time Points")
    whichTime <- sort(as.double(unique(Meta$ExpData$samples$Time)))
  }
  
  
  
  if (!is.null(Motifs)){
    if(!class(Motifs)=='MotifList')
      stop("Motifs has to be of class 'MotifList'")
    if (is.null(Genome(MD)))
      stop("genome needed to find motifs")
  }
  if (!is.null(anno)){
    if(!class(anno)=='list' & !class(anno)=='GRanges')
      stop("anno has to be a 'GRanges' object or a list of GRanges objects")
    if (class(anno)=='GRanges')
      anno <- list(anno=anno)
  }
  
  ################
  ## Define layouts
  ################
  numplots.row <- length(whichTime)
  widths <- rep(c(.1,rep(.3,numplots.col),numplots.row)) 
  numplots.col <- length(whichCondition)+1
  numplots <- numplots.row *numplots.col
  heigts <- rep(.2,numplots)
  
  par(oma = c(4,1,0,0) + 0.1,
      mar = c(0,0,1,1) + 0.1)
  leg.size <- .6
  
  
  if (!is.null(Motifs)){
    numplots.row <- numplots.row + 1
    widths <- c(widths,width(1:numplots.col))
    heigts <- c(heigts,heights(1:numplots.col))
  }
  if (!is.null(anno)){
    numplots.row <- numplots.row + 1
    widths <- c(widths,width(1:numplots.col))
    heigts <- c(heigts,heights(1:numplots.col))
  }
  numplots <- numplots.row *numplots.col
  
  layout(matrix(seq(1,numplots),  nrow = numplots.row, 
                ncol = numplots.col, byrow = TRUE),
           widths=widths, heights=heigts)
  leg.size <- 0.8
 
  
  if (is.null(xaxt)){
    xlab=PeakName
  } else{
    xlab=''
  }
  
  #first plot: condition 1 and time 1
  for (t.idx in 1:length(whichTime)){
    for (cond.idx in 1:length(whichCondition)){
      Sample.idx <- which(samples$Condition==whichCondition[cond.idx]&samples$Time==whichTime[t.idx])
      Sample.ids <- samples$SampleID[Sample.idx]
      SamplePeaks <- Peak[Sample.idx,,drop=FALSE]
      nSamples <- length(Sample.idx)

      ##Normalization
      if (is.null(NormMethod)){
        message('No normalization factors applied')
        Factors <- rep(1,length(Sample.idx))
        ylab <- 'counts'
      } else {
        message('normalization currently not implemented. Set NormMethod=NULL')
        NormFactors <- MD$NormFactors[[NormMethod]][[1]]
        ## check if factor is available for all relevant samples
        if (!all( is.element(rownames(SamplePeaks),names(NormFactors)))){
          tmp <- which(!is.element(rownames(SamplePeaks),names(NormFactors)))
          message('no normalization factor found for sample: ',
                  rownames(SamplePeaks)[tmp])
          stop('check normalization factors')
        }
        Factors <- NormFactors[sapply(rownames(SamplePeaks),function(name){
          idx=which(name==names(NormFactors))})]
        ylab='normalized counts'
      }
    }
  }
  
  ##Normalization
  if (is.null(NormMethod)){
    message('No normalization factors applied')
    Factors <- rep(1,length(Sample.idx))
    if (add.contplot==1){
      Factors2 <- rep(1,length(Sample.idx2))
    }
    ylab <- 'counts'
  } else {
    message('normalization currently not implemented. Set NormMethod=NULL')
    NormFactors <- MD$NormFactors[[NormMethod]][[1]]
    ## check if factor is available for all relevant samples
    if (!all( is.element(rownames(SamplePeaks),names(NormFactors)))){
      
      tmp <- which(!is.element(rownames(SamplePeaks),names(NormFactors)))
      message('no normalization factor found for sample: ',
              rownames(SamplePeaks)[tmp])
      stop('check normalization factors')
    }
    Factors <- NormFactors[sapply(rownames(SamplePeaks),function(name){
      idx=which(name==names(NormFactors))})]
    ylab='normalized counts'
  }
  
  ##colors, line type, legend text
  lty <- c(rep(1,nSamples))
  cols <- rainbow(nSamples,start=0,end=0.2)
  leg <- Sample.ids
  title <- ''
  if (add.contplot==1){
    lty2 <- c(rep(1,nSamples2))
    cols2 <- rainbow(nSamples2,start=0.5,end=0.7)
    leg2 <- Sample.ids2
    title <- contrast$name1
    title2 <- contrast$name2
  }
  
  
  ##INPUT
  if (plot.input){
    inputs <- unique(samples$ControlID[Sample.idx])
    input.idx <- sapply(inputs,function(id){ix=which(rownames(Peak)==id)})
    if (length(input.idx)==0){
      stop('no inputs given. Check column ControlID
           in your SampleSheet.csv file')
    }
    InputPeaks <- Peak[input.idx,, drop=FALSE]
    leg <- c(leg,inputs)
    cols <- c(cols,rep('black',length(input.idx)))
    lty <- c(lty,rep(2,length(input.idx)))
    
    
    if (add.contplot==1){
      inputs2 <- unique(samples$ControlID[Sample.idx2])
      input.idx2 <- sapply(inputs2,function(id){ix=which(rownames(Peak)==id)})
      InputPeaks2 <- Peak[input.idx2,, drop=FALSE]
      leg2 <- c(leg2,inputs2)
      cols2 <- c(cols2,rep('black',length(input.idx)))
      lty2 <- c(lty2,rep(2,length(input.idx)))
    }
    }
  
  ## get genome coordinates for x-values.
  
  start <- as.integer(strsplit(strsplit(PeakName,':')[[1]][2], '-')[[1]][1])
  end <- as.integer(strsplit(strsplit(PeakName,':')[[1]][2], '-')[[1]][2])
  chr <- strsplit(PeakName,':')[[1]][1]
  
  x.coords <- as.integer(colnames(SamplePeaks)) + start - PeakBoundary
  
  ################
  ### FIND Overlapping MOTIFS
  ################
  if (!is.null(Motifs)) {
    package <- Genome(MD)
    message('checking availability of package ',package)
    available <- suppressMessages(suppressWarnings(
      sapply(package, require, quietly = TRUE, character.only = TRUE,
             warn.conflicts = FALSE)))
    if (available){
      Genome <- getBSgenome(Genome(MD))
      chrom_seq <- Genome[[chr]]
      Seq <- Views(chrom_seq,start- PeakBoundary,end+ PeakBoundary)
      TFBS <- vector("list",length(Motifs))
      names(TFBS) <- names(Motifs)
      TFBS2 <- TFBS
      for (j in 1:length(Motifs)){
        TFBS[[j]] <- matchPWM(Motifs[[j]],Seq,Motifcutoff)
        TFBS2[[j]] <- matchPWM(reverseComplement(Motifs[[j]]),Seq,Motifcutoff)
      }
    } else {
      TFBS = vector("list",0)
      TFBS2 = vector("list",0)
      warning('Genome package is not available, needs to be installed
              first in order to find TF Binding sites')
    }
    
    }
  
  ################
  ### FIND Overlapping Annotation
  ################
  
  if (!is.null(anno)) {
    IRanges <- NULL
    rm(IRanges)
    
    Reg <- GRanges(seqnames=chr,IRanges(start- PeakBoundary,end+ PeakBoundary))
    GR <- vector("list",length(anno))
    names(GR) <- names(anno)
    for (j in 1:length(anno)){
      ov <- findOverlaps(anno[[j]],Reg)
      GR[[j]] <- anno[[j]][unique(queryHits(ov))]
      
    }
  }
  
  if (is.null(ylim)){
    # ylim <- c(min(SamplePeaks)-.4*(max(SamplePeaks)-min(SamplePeaks)),
    # 1.2*max(SamplePeaks))
    if ((add.contplot==0)){
      ylim <- c(min(SamplePeaks),1.2*max(SamplePeaks))
    } else {
      ylim <- c(min(c(SamplePeaks,SamplePeaks2)),
                1.2*max(c(SamplePeaks,SamplePeaks2)))
    }
  }
  
  if (is.null(xlim)){
    xlim <- c(start- PeakBoundary,end+ PeakBoundary)
  }
  
  ## first plot legend
  plot.new()
  
  legend("topleft",legend=leg,col=cols,lty=lty,lwd=1.4,inset=.01,cex=leg.size)
  
  ### plot Peaks (group1 if contrast given)
  plot(x.coords,SamplePeaks[1,]/Factors[1],
       t='l',col=cols[1],ylim=ylim,xlim=xlim,
       ylab=ylab,main=title,bty='n',
       cex.main=1.5,cex=1,cex.lab=1,cex.axis=1,lwd=1.4,xaxt=xaxt.1)
  if (nSamples>1){
    for (i in 2:nSamples){
      points(x.coords,SamplePeaks[i,]/Factors[i],t='l',col=cols[i],lwd=1.4)
    }
  }
  if (plot.input){
    for (i in 1:length(input.idx)){
      points(x.coords,InputPeaks[i,],lty=2,t='l',col='black',lwd=1.4)
    }
  }
  
  text(xlim[2]-0.1*(xlim[2]-xlim[1]),
       ylim[2]-0.1*(ylim[2]-ylim[1]),
       paste('Peak id:', Peak.id))
  text(xlim[2]-0.1*(xlim[2]-xlim[1]),
       ylim[2]-0.2*(ylim[2]-ylim[1]),
       paste('Fragment', whichPos))
  
  grid()
  
  # plot group2
  
  if (add.contplot==1){
    ## first plot legend
    plot.new()
    legend("topleft",legend=leg2,col=cols2,lty=lty2,lwd=1.4,inset=.01,
           cex=leg.size)
    
    ### plot Peaks (group2 if contrast given)
    plot(x.coords,SamplePeaks2[1,]/Factors2[1],
         t='l',col=cols2[1],ylim=ylim,xlim=xlim,
         ylab=ylab,main=title2,bty='n',
         cex.main=1.5,cex=1,cex.lab=1,cex.axis=1,lwd=1.4,xaxt=xaxt.2)
    if (nSamples2>1){
      for (i in 2:nSamples2){
        points(x.coords,SamplePeaks2[i,]/Factors2[i],t='l',col=cols2[i],lwd=1.4)
      }
    }
    if (plot.input){
      for (i in 1:length(input.idx2)){
        points(x.coords,InputPeaks2[i,],lty=2,t='l',col='black',lwd=1.4)
      }
    }
    grid()
    
  }
  
  ## ---------
  ## PLOT MOTIFS
  ## ---------
  if (!is.null(Motifs)&& (length(TFBS)>0 | length(TFBS2)>0)){
    #seqLogo(Motifs[[1]])
    ylim=c(0,length(TFBS)+1)
    plot(1, t='l',col='white',ylim=ylim,xlim=c(0,100),
         bty='n',
         cex=1.5,cex.lab=1.5,cex.axis=1.5,yaxt='n',xaxt='n')
    for (j in 1:length(TFBS)){
      text(0,j,names(TFBS[j]),adj=c(0,0),cex=.8)
    }
    plot(1, t='l',col='black',ylim=ylim,xlim=xlim,
         bty='n',
         cex=1.5,cex.lab=1.5,cex.axis=1.5,yaxt='n',xaxt=xaxt.3)
    
    grid()
    for (j in 1:length(TFBS)){
      tfbs = TFBS[[j]]
      tfbs2 = TFBS2[[j]]
      if (length(tfbs)==0 & length(tfbs2)==0  ) next
      
      if (!is.null(tfbs) & length(tfbs)>0){
        for (k in 1:length(tfbs)){
          lines(c(start(tfbs)[k],end(tfbs)[k]),c(j,j),lw=5)
        }}
      
      if (!is.null(tfbs2) & length(tfbs2)>0){
        for (k in 1:length(tfbs2)){
          lines(c(start(tfbs2)[k],end(tfbs2)[k]),c(j,j),lw=5,col='red')
        }
      }
    }
  }
  
  ####GENE Annotation
  if (!is.null(anno)){
    ylim=c(0,length(GR)+1)
    plot(1, t='l',col='white',ylim=ylim,xlim=c(0,100),
         bty='n',
         cex=1.5,cex.lab=1.5,cex.axis=1.5,yaxt='n',xaxt='n')
    for (j in 1:length(GR)){
      text(0,j,names(GR)[j],adj=c(0,0),cex=1)
    }
    plot(1, t='l',col='black',ylim=ylim,xlim=xlim,
         bty='n',
         cex=1.5,cex.lab=1.5,cex.axis=1.5,yaxt='n',xaxt=xaxt.4)
    grid()
    for (j in 1:length(GR)){
      gr = GR[[j]]
      
      if (length(gr)==0) next
      
      for (k in 1:length(gr)){
        lines(c(start(gr)[k],end(gr)[k]),c(j,j),lw=5,col='blue')
        if (as.vector(strand(gr))[k]=='+'){
          points(c(start(gr)[k],end(gr)[k]),c(j,j),cex=3,col='blue',pch='>')
        }
        if (as.vector(strand(gr))[k]=='-'){
          points(c(start(gr)[k],end(gr)[k]),c(j,j),cex=3,col='blue',pch='<')
        }
        if (k==1){
          ov <- intersect(IRanges(start(gr[k]),end(gr[k])),
                          IRanges(xlim[1],xlim[2]))
          text(c(start(ov)+.5*width(ov)),j+0.5,names(gr)[k])
        }
        
      }
    }
  }
  mtext(text=xlab,side=1,line=3)
}

