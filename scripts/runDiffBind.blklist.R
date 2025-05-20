#RY: DiffBind3 has gone through some major updates. This version only takes ~2min to run.

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("USAGE: runDiffBind.R diffBind_bins_BINSIZE_ BASENAME", call.=FALSE)
}

library("DiffBind")
library("BiocParallel")
library("foreach")
library("genomation")
#register(DoparParam())
#registered()
#bpparam("SerialParam")

prefix <- args[1]
bn <- args[2]

fileNames <- list.files(pattern = paste(prefix, "[1-9].txt", sep=""))
blacklist <-  readBed("/home/ruofany/diffbind/hg38-blacklist.v2.bed")
n = 1
for (file in fileNames){
    SON <- dba(sampleSheet=file, minOverlap=1)
    #RY: Set design to false to keep consistency with previous diffbind version.
    SON$config$design <-FALSE
    SON$config$doGreylist <-FALSE
    SON$config$factor <- "normalize"
    SON <- dba.blacklist(SON, blacklist = blacklist, greylist = FALSE)
    peakdata <- dba.show(SON)$Intervals
    #RY: Added parameters to remove RPKM filter etc. (diffbind3 filters for target with rpkm>1) and speedup calculation.
    SON <- dba.count(SON,  summits = FALSE, bUseSummarizeOverlaps=FALSE,filter=0,bParallel=TRUE)
    SON <- dba.normalize(SON,normalize=DBA_NORM_RLE, library=DBA_LIBSIZE_PEAKREADS, background=FALSE)
    SON <- dba.contrast(SON,categories = DBA_CONDITION, minMembers = 2)
    SON <- dba.analyze(SON)
    SON.DB <- dba.report(SON, th=1)
    write.table(as.data.frame(SON.DB), sep="\t", file=paste("DiffBindResults_", bn, n, ".txt", sep=""),  quote=FALSE)
    counts <- dba.peakset(SON, bRetrieve = TRUE)
    write.table(as.data.frame(counts), sep="\t", file=paste("counts_", bn, n, ".txt", sep=""), quote=FALSE)
    dba.save(SON, paste(bn, n, ".txt", sep=""))
    n = n + 1
}
