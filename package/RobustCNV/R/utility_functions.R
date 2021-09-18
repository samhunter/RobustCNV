#
# A set of utility & helper functions not core to RobustCNV.

#Roxygen docs:
#' Create an output formatted for Nexus compatibility.
#'
#' @param perIntervalCalls a data.frame containing per-interval calls
#' @return A data.frame formatted for Nexus compatibility
#' @examples
#' \dontrun{
#' formatForNexus(perIntervalCalls)
#' }
formatForNexus = function(perIntervalCalls){
  nexus = data.frame(Chromosome=perIntervalCalls$chrom,
                     Start=perIntervalCalls$start,
                     End=perIntervalCalls$stop,
                     Value=perIntervalCalls$Normalized_Coverage,
                     Sample=perIntervalCalls$Sample,
                     Gene=perIntervalCalls$name,
                     GCpct=perIntervalCalls$pGC,
                     Index=perIntervalCalls$index,
                     Call=perIntervalCalls$Call,
                     SegmentMean=perIntervalCalls$SegmentMean )
  return(nexus)
}

#Roxygen docs:
#' Summarize interval level calls to gene level call. This is an interval function that
#' is used by summarizeToGeneCalls.
#'
#' @param x a vector of calls for a single gene
#' @return A data.frame containing QC metrics
#' @examples
#' \dontrun{
#' CNVqc(data)
#' }
summarizeIntervalCalls = function(x){
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
#' Summarize per-IntervalCalls to gene level calls.
#'
#' @param perIntervalCalls a vector of calls for a single gene
#' @return A data.frame containing calls summarized to gene level
#' @examples
#' \dontrun{
#' CNVqc(data)
#' }
summarizeToGeneCalls = function(perIntervalCalls){
  pergene.out = aggregate(perIntervalCalls$Call,
                          by=list(perIntervalCalls$Sample, perIntervalCalls$Gene), F=summarize_calls)
  pergene = cbind(pergene.out[,1:2], pergene.out[,3][,1], pergene.out[,3][,2], pergene.out[,3][,3])
  colnames(pergene) = c("Sample","Gene","IntervalPattern","CompressedCalls","Gene_call")

  med_seg = aggregate(as.numeric(perIntervalCalls$SegmentMean),
                      by=list(perIntervalCalls$Sample, perIntervalCalls$Gene), F=median, na.rm=T)
  colnames(med_seg) = c("Sample","Gene","GeneSegmentMedian")

  pergene$GeneSegmentMedian = signif(med_seg$GeneSegmentMedian, 3)

  #Sort by sample and gene:
  geneOrder = 1:length(unique(perIntervalCalls$Gene))
  #this extracts the index of only the first instance of each gene:
  idx = !duplicated(perIntervalCalls$Gene)

  #get chromosome and start position for each gene
  gChx = gcstats$chrom[idx]
  gStart = gcstats$start[idx]
  geneOrder = data.frame(geneOrder, gChx, gStart)
  geneOrder$Gene = perIntervalCalls$Gene[idx]

  #order the calls:
  pergene = pergene[order(pergene$Sample, geneOrder$geneOrder[match(pergene$Gene, geneOrder$Gene)] ), ]
  pergene = cbind(geneOrder[match(pergene$Gene, geneOrder$Gene), ], pergene)
  pergene2 = data.frame(Sample=pergene$Sample, Gene=pergene$Gene,
                        Chromosome=pergene$gChx, geneStart=pergene$gStart,
                        CompressedIntervalCalls=pergene$CompressedCalls,
                        GeneCall=pergene$Gene_call, GeneSegmentMedian=pergene$GeneSegmentMedian,
                        IntervalPattern=pergene$IntervalPattern)

  pergene2 = pergene2[grep("NA_", pergene2$Gene, invert=T),]
  return(pergene2)
}
