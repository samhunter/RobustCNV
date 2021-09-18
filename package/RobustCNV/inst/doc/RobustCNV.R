## ------------------------------------------------------------------------
require(RobustCNV)
sinf = function(x, sd=0, CNV=F, begin=0, end=0, magnitude=1, gc){
  #add systematic and random noise
  m= 3*(sin(x/100)) + rnorm(n=1, mean=0, sd=sd)
  m = m + 20
  #Simulate gc bias:
  m = m * (1 + gc^3 - gc^4)
  #add CNV
  if(CNV & x >= begin & x <= end){
    m = m + abs(m) * magnitude
  }
  
  return(m)
}

## ------------------------------------------------------------------------

intervals = data.frame(Chrom=rep(1:23, each=250), 
                      Gene=rep(paste('gene', 1:1150, sep=''), each=5),
                      pGC=sample(x=seq(0,1, .001), size=5750, replace = T),
                      Start=rep(1:250, 23) * 120,
                      End=rep(1:250, 23) * 120 + 70,
                      Index=1:(250*23))

normals = data.frame(n1 = sapply(seq(1,5750, 1), function(x){sinf(x, sd=1, gc=intervals$pGC[x])}), 
                     n2 = sapply(seq(1,5750, 1), function(x){sinf(x, sd=1, gc=intervals$pGC[x])}), 
                     n3 = sapply(seq(1,5750, 1), function(x){sinf(x, sd=1, gc=intervals$pGC[x])}),
                     n4 = sapply(seq(1,5750, 1), function(x){sinf(x, sd=1, gc=intervals$pGC[x])}),
                     n5 = sapply(seq(1,5750, 1), function(x){sinf(x, sd=1, gc=intervals$pGC[x])}))

high_gain_tumor = unlist(sapply(seq(1,5750, 1), function(x)
                                  {sinf(x, sd=1, CNV=T, begin=400, end=500, magnitude=1, gc=intervals$pGC[x])}))
low_gain_tumor = unlist(sapply(seq(1, 5750, 1), function(x)
                                 {sinf(x, sd=1, CNV=T, begin=400, end=500, magnitude=.5,gc=intervals$pGC[x])}))
low_loss_tumor = unlist(sapply(seq(1, 5750, 1), function(x)
                                 {sinf(x, sd=1, CNV=T, begin=50, end=150, magnitude=-.5, gc=intervals$pGC[x])}))
aneuploid_tumor = c(unlist(sapply(seq(1,1000, 1), function(x)
                                  {sinf(x, sd=1, CNV=T, begin=0, end=500, magnitude=.5, gc=intervals$pGC[x])})),
                    unlist(sapply(seq(1001,2000, 1), function(x)
                                  {sinf(x, sd=1, CNV=T, begin=1001, end=1500, magnitude=.5, gc=intervals$pGC[x])})),
                    unlist(sapply(seq(2001,3000, 1), function(x)
                                  {sinf(x, sd=1, CNV=T, begin=2001, end=2500, magnitude=-.5, gc=intervals$pGC[x])})),
                    unlist(sapply(seq(3001,4000, 1), function(x)
                                  {sinf(x, sd=1, CNV=T, begin=3001, end=3500, magnitude=.5, gc=intervals$pGC[x])})),
                    unlist(sapply(seq(4001,5000, 1), function(x)
                                  {sinf(x, sd=1, CNV=T, begin=4001, end=4500, magnitude=.5, gc=intervals$pGC[x])})),
                    unlist(sapply(seq(5001,5750, 1), function(x)
                                  {sinf(x, sd=1, CNV=T, begin=5001, end=5500, magnitude=.5, gc=intervals$pGC[x])}))
                    )

aneuploid_tumor_alleles = c(sample(c(2/3, 1/3), replace = T, size = 500), rep(.5, 500),
                            sample(c(2/3, 1/3), replace = T, size = 500), rep(.5, 500),
                            rep(1, 500), rep(.5, 500),
                            sample(c(2/3, 1/3), replace = T, size = 500), rep(.5, 500),
                            sample(c(2/3, 1/3), replace = T, size = 500), rep(.5, 500),
                            sample(c(2/3, 1/3), replace = T, size = 500), rep(.5, 250)
                            )

aneuploid_tumor_alleles = data.frame(Index=intervals$Index, af=aneuploid_tumor_alleles)

## ---- fig.show='hold'----------------------------------------------------
plotGCBias("GC bias plot", intervals$pGC, normals[,1])
plotGCBias("GC bias plot", intervals$pGC, high_gain_tumor)

## ---- fig.width=6, fig.height=5------------------------------------------
plotSampleClusters(rawNormalCoverage=normals, 
                   rawTumorCoverage=cbind(high_gain_tumor, low_gain_tumor, low_loss_tumor, aneuploid_tumor))

## ------------------------------------------------------------------------
high_gain_tumor.normalized <- robustNorm(high_gain_tumor, normals, pGC=intervals$pGC, EVs=rep(100, 5750))
low_gain_tumor.normalized <- robustNorm(low_gain_tumor, normals, pGC=intervals$pGC, EVs=rep(100, 5750))
low_loss_tumor.normalized <- robustNorm(low_loss_tumor, normals, pGC=intervals$pGC, EVs=rep(100, 5750))
aneuploid_tumor.normalized <- robustNorm(aneuploid_tumor, normals, pGC=intervals$pGC, EVs=rep(100, 5750))

## ---- fig.show='hold'----------------------------------------------------
plot(high_gain_tumor, main="High gain raw data", pch='.', xlab='index', ylab='')
plot(high_gain_tumor.normalized, main="High gain normalized", pch='.', xlab='index', ylab='')
abline(h=0)
plotGCBias("GC bias plot", intervals$pGC, high_gain_tumor.normalized)

## ---- fig.show='hold'----------------------------------------------------
plot(low_gain_tumor, main="Low gain raw data", pch='.', xlab='index', ylab='')
plot(low_gain_tumor.normalized, main="Low gain normalized", pch='.', xlab='index', ylab='')
abline(h=0)
plotGCBias("GC bias plot", intervals$pGC, low_gain_tumor.normalized)

## ---- fig.show='hold'----------------------------------------------------
plot(low_loss_tumor, main="Low loss raw data", pch='.', xlab='index', ylab='')
plot(low_loss_tumor.normalized, main="Low loss normalized", pch='.', xlab='index', ylab='')
abline(h=0)
plotGCBias("GC bias plot", intervals$pGC, low_gain_tumor.normalized)

## ---- fig.show='hold'----------------------------------------------------
plot(aneuploid_tumor, main="Aneuploid raw data", pch='.', xlab='index', ylab='')
plot(aneuploid_tumor.normalized, main="Aneuploid normalized", pch='.', xlab='index', ylab='')
abline(h=0)
plotGCBias("GC bias plot", intervals$pGC, low_gain_tumor.normalized)

## ---- fig.show='hold'----------------------------------------------------
aneuploid_tumor.recentered = allelicRecentering(aneuploid_tumor.normalized, intervals, aneuploid_tumor_alleles)

plot(aneuploid_tumor, main="Aneuploid raw data", pch='.', xlab='index', ylab='')
plot(aneuploid_tumor.recentered, main="Aneuploid normalized and recentered", pch='.', xlab='index', ylab='')
abline(h=0)
plotGCBias("GC bias plot", intervals$pGC, low_gain_tumor.normalized)

## ------------------------------------------------------------------------

high_gain_tumor.segmented <- segment(high_gain_tumor.normalized, 
                                     orderedChrom=as.ordered(intervals$Chrom), maploc=intervals$End)
low_gain_tumor.segmented <- segment(low_gain_tumor.normalized, 
                                    orderedChrom=as.ordered(intervals$Chrom), maploc=intervals$End)
low_loss_tumor.segmented <- segment(low_loss_tumor.normalized, 
                                    orderedChrom=as.ordered(intervals$Chrom), maploc=intervals$End)
aneuploid_tumor.segmented <- segment(aneuploid_tumor.recentered, 
                                    orderedChrom=as.ordered(intervals$Chrom), maploc=intervals$End)


## ------------------------------------------------------------------------
high_gain_tumor.called <- callSegments(high_gain_tumor.segmented, 
                                       high_gain_tumor.normalized, tuning=0.8)
low_gain_tumor.called <- callSegments(low_gain_tumor.segmented, 
                                      low_gain_tumor.normalized, tuning=0.8)
low_loss_tumor.called <- callSegments(low_loss_tumor.segmented, 
                                      low_loss_tumor.normalized, tuning=0.8)
aneuploid_tumor.called <- callSegments(aneuploid_tumor.segmented, 
                                      aneuploid_tumor.recentered, tuning=0.8)


## ---- fig.width=7, fig.height=5------------------------------------------
high_gain_tumor.called.intervals <- assignIntervalCalls(intervals, high_gain_tumor.called, high_gain_tumor.normalized)
low_gain_tumor.called.intervals <- assignIntervalCalls(intervals, low_gain_tumor.called, low_gain_tumor.normalized)
low_loss_tumor.called.intervals <- assignIntervalCalls(intervals, low_loss_tumor.called, low_loss_tumor.normalized)
aneuploid_tumor.called.intervals <- assignIntervalCalls(intervals, aneuploid_tumor.called, aneuploid_tumor.recentered)

plotGenomeIntervalCalls(perIntervalCalls=high_gain_tumor.called.intervals, main='High gain tumor')
plotGenomeIntervalCalls(perIntervalCalls=low_gain_tumor.called.intervals, main='Low gain tumor')
plotGenomeIntervalCalls(perIntervalCalls=low_loss_tumor.called.intervals, main='Low loss tumor')
plotGenomeIntervalCalls(perIntervalCalls=aneuploid_tumor.called.intervals, main='Aneuploid tumor')


