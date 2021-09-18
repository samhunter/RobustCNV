#Roxygen docs:
#' Plot GC bias plot
#'
#' @param sampleName sample name
#' @param pGC vector of percent GC for each interval
#' @param coverage vector of per-interval coverage
#' @return A gc-score value which is the sum of abs(predicted values).
#' @examples
#' \dontrun{
#' plotGCBias(sampleName, pGC, coverage)
#' }
plotGCBias = function(sampleName, pGC, coverage){
  #cat(paste(s, '\n'))
  i = seq(0, 1, by=0.001)
  lfit = loess(coverage ~ pGC, span=0.03)
  #smooth the model a bit:
  lfit2 = loess(predict(lfit, i) ~ i, span = 0.3)
  #calculte the gc-score:
  gcs = signif(sum(abs(predict(lfit2, seq(.3, .7, by=.001)))), 2)
  #png(file=paste(plot_path, "/GCbias_", gcs, "_", s, ".png", sep=''), width=1000, height=500)
  plot(y=coverage, pGC, pch=20,
       col=rgb(20,20,20,60, maxColorValue=255), ylab='Coverage', xlab='% GC',
       main=paste('GC bias:', sampleName, '\nGC bias score:', gcs))
  lines(predict(lfit2, i), x=i, col='red', lwd=2)
  #dev.off()
  return(gcs)
}

#Roxygen docs:
#' Make sample-cluster plots.
#'
#' @param rawNormalCoverage Data.frame containing raw coverage for normal samples
#' @param rawTumorCoverage Data.frame containing raw coverage for tumor samples
#' @return nothing.
#' @examples
#' \dontrun{
#' plotSampleClusters(rawNormalCoverage, rawTumorCoverage)
#' }
plotSampleClusters = function(rawNormalCoverage, rawTumorCoverage){
  #require(WGCNA)
  sampleTypeColors = c(rep('red', dim(rawNormalCoverage)[2]), rep('blue', dim(rawTumorCoverage)[2]))
  allSamples = c(colnames(rawNormalCoverage), colnames(rawTumorCoverage))
  tmpid = sapply(strsplit(allSamples, '_'), '[', 1)

  hc = hclust(as.dist(1-cor(cbind(rawNormalCoverage, rawTumorCoverage), method='spearman')))
  colors = cbind(sampleType=sampleTypeColors)
  layout(matrix(c(1,2), nrow = 2), heights = c(80,20))
  op = par(mar=c(1,6,4,.5))
  plot(x=hc, main="Heirarchical clustering, normal and tumor samples\nSpearman correlation of raw mapping depth")
  par(mar=c(1,6,0,.5))
  WGCNA::plotColorUnderTree(hc, colors=colors, rowLabels=colnames(colors), addTextGuide=T, cex.rowLabels=.6)
  par(op)
}

#Roxygen docs:
#' Make genome plots.
#'
#' @param perIntervalCalls data.frame containing per-Interval calls and other data.
#' @param main title of plot.
#' @return nothing.
#' @examples
#' \dontrun{
#' plotGenomeIntervalCalls(perIntervalCalls)
#' }
plotGenomeIntervalCalls = function(perIntervalCalls, main=''){
  # Plot a sample: #blue is loss, red is gain:
  simucolor = c("0"=rgb(20,20,20,120, maxColorValue=255),
                "-"=rgb(0,0,150,200, maxColorValue=255),
                "+"=rgb(150,0,0,200, maxColorValue=255))

  limited_coverage = perIntervalCalls$log2ratio
  pch = rep(20, length(limited_coverage))
  pch[limited_coverage > 2] = 24
  pch[limited_coverage < -2] = 25
  limited_coverage[limited_coverage > 2] = 2
  limited_coverage[limited_coverage < -2] = -2
  plot(x=perIntervalCalls$Index, y=limited_coverage,
       pch=pch, col=simucolor[perIntervalCalls$Segment_Call], bg=simucolor[perIntervalCalls$Segment_Call],
       main=main,
       cex=.9, ylab='Log2 Copy Ratio', xlab='Chromosome', ylim=c(-2, 2), xaxt='n', font.lab=2, font.axis=2)
  abline(h=0)
  chx.idx = which(!duplicated(perIntervalCalls$Chrom))
  abline(v=chx.idx, col="grey50", lwd=2)
  axis(side=1, at=chx.idx, labels=perIntervalCalls$Chrom[chx.idx], cex=.9, font=2)
  points(perIntervalCalls$Segment_Mean, x=perIntervalCalls$Index, col='orange', pch='.', cex=2)
  points(perIntervalCalls$Cutoff, x=perIntervalCalls$Index, col='green', pch='-', cex=1)
  points(-perIntervalCalls$Cutoff, x=perIntervalCalls$Index, col='green', pch='-', cex=1)
  legend("topright", col=c(simucolor, 'orange', 'green'), pch=c(20,20,20,20), legend=c(names(simucolor), 'SegMean', 'Cutoff'), pt.cex=1.2)
}

