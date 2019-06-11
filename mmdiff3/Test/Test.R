
setwd('~/Edinburgh Biocourse/Trieste_Data/ChIP-Seq/')

library(devtools)
install("~/Documents/GitHub/mmdiff3/")

library(DiffBind)
#sampleSheet = "SampleSheetPromoters.csv"

sampleSheet = "SampleSheet.csv"
DBA <- dba(sampleSheet=sampleSheet)
DBA <- dba.count(DBA)

DBA <- dba.contrast(DBA,categories=DBA_TISSUE,minMembers=2)
DBA <- dba.analyze(DBA)

#------------
Peaks = dba.peakset(DBA,bRetrieve=TRUE)

#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("BSgenome.Hsapiens.UCSC.hg38")

library('MMDiff3')
ExperimentData <- list(genome='BSgenome.Hsapiens.UCSC.hg38',
                       dataDir=".",
                       sampleSheet = "SampleSheet.csv")

MetaData <- list('ExpData' = ExperimentData)
MMD <- DBAmmd(MetaData)

MMD <- setRegions(MMD, Peaks)
MMD <- getPeakReads(MMD)
MMD <- estimateFragmentCenters(MMD)
MMD <- compHists(MMD)
#MMD<- compDists(MMD, dist.method = "MMD2")
MMD<- compDists(MMD, dist.method = "MMD2")
MMD<- compDists(MMD, dist.method = "MMD")
MMD <- setContrast(MMD, contrast='byTissue')
MMD <- compPvals(MMD, dist.method = "MMD2")
res <- reportResults(MMD)

plotPeak(MMD, Peak.id='241', plot.input = FALSE, whichPos="Center")
