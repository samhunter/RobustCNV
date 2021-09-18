#Core functions for RobustCNV normalization, segmentation, and calling

#Roxygen docs:
#' Normalize coverage data.
#'
#' Do a normalization with interval length, gcContent and panel of normals using a robust regression strategy
#' followed by a GC correction using loess.
#'
#' @param y a vector of coverage values for the sample of interest.
#' @param x a dataframe or matrix with coverage values for normal samples (columns are samples).
#' @param EVs a vector or dataframe of explanatory variables (e.g. interval lengths).
#' @param pGC, the percent GC content for each interval
#' @param idx a vector used to subset the input during model generation only (all intervals are reported)
#' @return A new vector normalized against the panel of normals, interval lengths, and for GC content
#' @examples
#' \dontrun{
#' robustNorm(y, x, ilengths, pGC)
#' }
robustNorm = function(y, x, pGC, EVs=NA, idx=NA){
  #take care of situations where no idx is passed in
  if(length(idx) == 1){ if(is.na(idx)){ idx <- 1:nrow(x) } }

  #Add EVs to the data.frame
  if(!is.na(EVs)[1]){
    x2 = data.frame(x, EVs)
  } else { x2 = x }

  #make colnames of x2 generic and compatible with rlm:
  colnames(x2) <- paste('s', 1:dim(x2)[2], sep='')

  #fit robust linear model from MASS:
  rlm1 <- MASS::rlm(y[idx] ~ 0 + ., data=x2[idx,], maxit=30)

  #generate fitted values:
  p = predict(rlm1, x2)
  #handle strange situations where predicted is negative or 0
  if(sum(p <= 0) > 0){
    p[p <= 0] = median(y)  # set values to make normalization robust
  }
  rNorm = log2(y / p)
  # Set all non-numeric values to 0:
  rNorm[is.infinite(rNorm)] = 0
  rNorm[is.nan(rNorm)] = 0
  rNorm[is.na(rNorm)] = 0

  #### Now remove any remaining GC bias ######
  lfit = loess(rNorm ~ pGC, span=0.03)
  i = seq(0, 1, by=0.001)
  #smooth the fit:
  lfit2 = loess(predict(lfit, i) ~ i, span=0.3)
  lfitPred = predict(lfit2, pGC)
  #fill in the NA values with local average:
  for(n in which(is.na(lfitPred))){
    covs = lfitPred[(pGC >(pGC[n] - .1)) & (pGC < (pGC[n] + .1))]
    lfitPred[n] = mean(covs, na.rm=T)
  }
  rNorm2 = rNorm - lfitPred
  return(rNorm2)
}

#Roxygen docs:
#' Segment normalized log2 coverage values for one or more samples.
#'
#' @param normCoverage a data.frame or matrix with normalized coverage. Columns should be samples.
#' @param orderedChrom a vector of ordered factors indicating the chromosome, must be nrow(norm.coverage) long.
#' @param maploc a vector containing chromosome locus for intervals (end of each interval)
#' @param idx segments can be called using a subset of targeted intervals, idx provides the subset.
#' @return A new dataframe with columns: Sample, Chromosome, Start, End, NumProbes, SegmentMean, SegmentCall
#' @examples
#' \dontrun{
#' segment(normalizedTumorCoverage, orderedChrom, maploc)
#' }
segment = function(normCoverage, orderedChrom, maploc, idx=NA){
  #Make segment work with vectors:
  if(is.null(dim(normCoverage))){
    normCoverage <- data.frame(normCoverage)
    colnames(normCoverage) <- c("s1")
  }
  #Set up an all-intervals index if no index is supplied
  if(length(idx) == 1){ if(is.na(idx)){ idx <- 1:dim(normCoverage)[1] } }

  #Create subsets:
  normCoverage <- normCoverage[idx,, drop=F]  # two commas is NOT a typo
  normCoverage2 <- 2^normCoverage  # the mean of a set of log2 values is not well defined
  orderedChrom <- orderedChrom[idx]
  maploc <- maploc[idx]

  #CNA screws up sample names, make a backup:
  samples <- colnames(normCoverage)
  names(samples) <- gsub('-', '.', samples)
  #TODO: re-write per-sample calling, figure out a way to tune segmentation by some normalization QC criterion.
  #make segment calls:
  CNAobject <- DNAcopy::CNA(normCoverage2, chrom=orderedChrom, maploc=maploc, data.type='logratio',
                   presorted=T, sampleid=colnames(normCoverage2))

  #Make results reproducible
  set.seed(68467046)
  #TODO evaluate the effect of smoothing
  smoothedCNA <- DNAcopy::smooth.CNA(CNAobject)
  #TODO DNAcopy supports many options for segmentation, the defaults may not be optimal.
  segmentSmoothedCNA <- DNAcopy::segment(smoothedCNA)

  returnValue <- data.frame(
    Sample=samples[segmentSmoothedCNA$output$ID],
    Chrom=segmentSmoothedCNA$out$chrom,
    Start=segmentSmoothedCNA$output$loc.start,
    End=segmentSmoothedCNA$output$loc.end,
    NumProbes=segmentSmoothedCNA$output$num.mark,
    SegmentMean=log2(segmentSmoothedCNA$output$seg.mean)  # put SegmentMean back in log2 space
  )

  #R does funny things with division by zero:
  returnValue$SegmentMean[returnValue$SegmentMean > 0 & is.infinite(returnValue$SegmentMean)] <- 0
  returnValue$SegmentMean[returnValue$SegmentMean < 0 & is.infinite(returnValue$SegmentMean)] <- -2
  return(returnValue)
}

#Roxygen docs:
#' Make calls on segments based on a cutoff calculated from segment SD and a tuning factor.
#'
#' @param gsegments a dataframe with columns: Sample,Chromosome,Start,End,NumProbes,Segment,SegmentMean
#' @param normCoverage a vector or dataframe containing normalized coverage
#' @param tuning a parameter for adjusting sensitivity, higher values are less sensitive, more specific.
#' @return A new dataframe with a "SegmentCall" column consisting of (-, 0, +).
#' @examples
#' \dontrun{
#' callSegments(gsegments, normCoverage)
#' callSegments(gsegments, normCoverage, tuning=.2)
#' }
callSegments = function(gsegments, normCoverage, tuning=.8){
  gsegments$SegmentCall = '0'
  #Make call_segments work with vectors:
  if(is.null(dim(normCoverage))){
    normCoverage=data.frame(normCoverage)
    colnames(normCoverage) = unique(gsegments$Sample)
  }

  for(s in unique(gsegments$Sample)){
    sidx = gsegments$Sample == s
    segmentidx <- c(1, cumsum(gsegments[sidx,]$NumProbes))
    tau <- median(sapply(1:(length(segmentidx) - 1), function(x){sd(normCoverage[segmentidx[x]: segmentidx[x+1], s])})) * tuning
    cat(paste(s, " Tau: ", tau, '\n'))
    calls <- gsegments$SegmentCall[sidx]
    calls[gsegments$SegmentMean[sidx] > tau] <- '+'
    calls[gsegments$SegmentMean[sidx] < -tau] <- '-'
    gsegments$SegmentCall[sidx] <- calls
    gsegments$Cutoff[sidx] <- tau
  }
  return(gsegments)
}

#Roxygen docs:
#' Given a set of segments with calls, intersect the capture intervals with these segments and assign
#' calls to them.
#'
#' @param intervals a dataframe with a definition of the intervals used in the capture experiment
#' @param gsegments dataframe returned from the callSegments() function.
#' @param normalized_tumor_coverage a vector returned from robustNorm() function.
#' @return a dataframe combining intervals and calls
assignIntervalCalls = function(intervals, gsegments, normalized_tumor_coverage){
  segmentGRs = GenomicRanges::GRanges(seqnames=gsegments$Chrom,
                                      ranges=IRanges::IRanges(start=gsegments$Start,
                                      end=gsegments$End), strand="*")
  intervalGRs = GenomicRanges::GRanges(seqnames=intervals$Chrom, ranges=IRanges::IRanges(start=intervals$Start,
                                      end=intervals$End), strand="*")
  ov = as.data.frame(GenomicRanges::findOverlaps(segmentGRs, intervalGRs))
  colnames(ov) = c("queryHits", "subjectHits")  #old versions of GRanges have other column names
  ov = ov[!duplicated(ov$subjectHits), ]
  interval.calls = intervals
  interval.calls$Segment_Call = '0'
  interval.calls$Segment_Mean = 0
  interval.calls$Cutoff = gsegments$Cutoff[1]
  interval.calls$Segment_Call[ov$subjectHits] = gsegments$SegmentCall[ov$queryHits]
  interval.calls$Segment_Mean[ov$subjectHits] = gsegments$SegmentMean[ov$queryHits]
  interval.calls$Segment_Index = 0
  interval.calls$Segment_Index[ov$subjectHits] = ov$queryHits
  interval.calls$log2ratio = normalized_tumor_coverage
  return(interval.calls)
}


#Roxygen docs:
#' Given a set of interval calls with gene information, use a set of rules to summarize them to
#' gene level calls.
#'
#' @param intervals a dataframe with a definition of the intervals used in the capture experiment
#' @param gsegments dataframe returned from the callSegments() function.
#' @param normalized_tumor_coverage a vector returned from robustNorm() function.
#' @return a dataframe combining intervals and calls
summarizeToGeneLevel = function(interval.calls){
  x[which(x == 'NA')] = '0'
  b = rle(unlist(x))
  icalls = paste(paste(b$lengths, b$values, sep=':'), collapse=', ')
  gcall = 'NormalCopy'
  tbl = table(unlist(x))
  if (('-' %in% names(tbl) & tbl['-'] > 2) & ('+' %in% names(tbl) & tbl['+']/sum(tbl) > .5)){
    gcall='gain+loss'
  } else if ('-' %in% names(tbl) & ((tbl['-'] > 2) | (tbl['-']/sum(tbl) == 1))){
    gcall='loss'
  } else if ('+' %in% names(tbl) & tbl['+']/sum(tbl) > .5) {
    gcall='gain'
  } else if('+' %in% names(tbl) & '-' %in% names(tbl)) {
    gcall='mixed'
  } else if('+' %in% names(tbl)){
    gcall='NormalCopy+'
  } else if('-' %in% names(tbl)){
    gcall='NormalCopy-'
  }
  #cat(paste(paste(x, collapse=''), " ", icalls, " ", gcall, '\n'))
  return(c(paste(x, collapse=''),icalls, gcall))
}


#Roxygen docs:
#' Some samples which are highly aneuploid are not well centered after normalization. This function
#' attempts to center them using allele fraction information to identify putatively diploid segments.
#'
#' @param normCoverage normalized log2 ratios from robustNorm()
#' @param intervals dataframe defining the intervals
#' @param sampleAlleles dataframe with allele fractions and index
#' @return a vector which has been re-centered.
allelicRecentering = function(normCoverage, intervals, sampleAlleles){
  AF = rep(NA, nrow(intervals))
  AF[match(sampleAlleles$Index, intervals$Index)] = sampleAlleles$af
  #set up weights, these values seem to work well for WES and targeted capture:
  w = rep(.001, nrow(intervals))
  w[AF < .6 & AF > .4] = 1
  w[AF > .45 & AF < .55] = 10
  w[AF > .49 & AF < .51] = 200
  adjustment = sum(normCoverage*w)/sum(w)
  cat(paste("adjustment:", adjustment, '\n'))
  adjusted = normCoverage - adjustment
  return(adjusted)
}

