% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutation.R
\name{runCOCOAPerm}
\alias{runCOCOAPerm}
\title{Run COCOA permutations to get p-values}
\usage{
runCOCOAPerm(
  genomicSignal,
  signalCoord,
  GRList,
  rsScores,
  targetVar,
  signalCol = c("PC1", "PC2"),
  scoringMetric = "default",
  absVal = TRUE,
  olList = NULL,
  centerGenomicSignal = TRUE,
  centerTargetVar = TRUE,
  variationMetric = "cor",
  nPerm = 300,
  useSimpleCache = TRUE,
  cacheDir = getwd(),
  dataID = "",
  testType = "greater",
  gammaFitMethod = "mme",
  realScoreInDist = TRUE,
  force = FALSE,
  verbose = TRUE,
  returnCovInfo = FALSE,
  ...
)
}
\arguments{
\item{genomicSignal}{Matrix/data.frame. 
The genomic signal (e.g. DNA methylation levels)
Columns of genomicSignal should be samples/patients. 
Rows should be individual signal/features
(each row corresponds to one genomic coordinate/range)}

\item{signalCoord}{A GRanges object or data frame with coordinates 
for the genomic signal/original epigenetic data. 
Coordinates should be in the 
same order as the original data and the feature contribution scores 
(each item/row in signalCoord
corresponds to a row in signal). If a data.frame, 
must have chr and start columns (optionally can have end column, 
depending on the epigenetic data type).}

\item{GRList}{GRangesList object. Each list item is 
a distinct region set to test (region set: regions that correspond to 
the same biological annotation). The region set database
must be from the same reference genome
as the coordinates for the actual data/samples (signalCoord).}

\item{rsScores}{data.frame. A data.frame with region set
scores. The output of the 'aggregateSignalGRList' function.
Each row is a region set. One column for each sample
variable of interest (e.g. PC or sample phenotype).
Also can have columns with info on the overlap between the 
region set and the epigenetic data. 
Rows should be in the same order as the region sets in GRList
(the list of region sets used to create rsScores.)}

\item{targetVar}{Matrix or data.frame. Rows should be samples. 
Columns should be the target variables 
(whatever variable you want to test for association with
the epigenetic signal: e.g. PC scores),}

\item{signalCol}{A character vector with the names of the sample variables
of interest/target variables (e.g. PCs or sample phenotypes). 

The columns in `sampleLabels` for which to calculate
the variation related to the epigenetic data
(e.g. correlation) and then to run COCOA on.}

\item{scoringMetric}{A character object with the scoring metric.
There are different methods available for 
signalCoordType="singleBase" vs  signalCoordType="multiBase".
For "singleBase", the available methods are "regionMean", 
"regionMedian", "simpleMean", and "simpleMedian". 
The default method is "regionMean".
For "multiBase", the methods are "proportionWeightedMean", 
"simpleMean", and "simpleMedian". The default is "proportionWeightedMean".
"regionMean" is a weighted
average of the signal, weighted by region (absolute value of signal 
if absVal=TRUE). First the signal is
averaged within each regionSet region, 
then all the regions are averaged. With
"regionMean" method, be cautious in interpretation for
region sets with low number of regions that overlap signalCoord. The
"regionMedian" method is the same as "regionMean" but the median is taken
at each step instead of the mean.
The "simpleMean"
method is just the unweighted average of all (absolute) signal values that
overlap the given region set. For multiBase data, this includes
signal regions that overlap a regionSet region at all (1 base
overlap or more) and the signal for each overlapping region is
given the same weight for the average regardless of how much it overlaps.
The "simpleMedian" method is the same as "simpleMean" but takes the median 
instead of the mean. 
"proportionWeightedMean" is a weighted average of all signalCoord 
regions that overlap with regionSet regions. For each signalCoord region
that overlaps with a regionSet region, we calculate what proportion
of the regionSet region is covered. Then this proportion is used to
weight the signal value when calculating the mean. 
The denominator of the mean
is the sum of all the proportion overlaps.}

\item{absVal}{Logical. If TRUE, take the absolute value of values in
signal. Choose TRUE if you think there may be some 
genomic loci in a region set that will increase and others
will decrease (if there may be anticorrelation between
regions in a region set). Choose FALSE if you expect regions in a 
given region set to all change in the same direction (all be positively
correlated with each other).}

\item{olList}{list. Each list item should be a "SortedByQueryHits" object 
(output of findOverlaps function). Each hits object should have the overlap
information between signalCoord and one item of GRList (one unique region set).
The region sets from GRList must be the "subject" in findOverlaps 
and signalCoord must be the "query". E.g. findOverlaps(subject=regionSet,
query=signalCoord).
Providing this information can greatly improve permutation speed since the 
overlaps will not have to be calculated for each permutation. 
The "runCOCOAPerm" function calculates this information only once, internally,
so this does not have to be provided when using that function. When using 
this parameter, signalCoord, 
genomicSignal, and each region set must be in the same order as they were
when olList was created. Otherwise, the wrong genomic loci will be referenced
(e.g. if epigenetic features were filtered out of genomicSignal after olList
was created.)}

\item{centerGenomicSignal}{Logical. Should rows in genomicSignal
be centered based on
their means? (subtracting row mean from each row)}

\item{centerTargetVar}{Logical. Should columns in targetVar be 
centered based
on their means? (subtract column mean from each column)}

\item{variationMetric}{Character. The metric to use to quantify the
association between each feature in genomicSignal and each target
variable in sampleLabels.
Either "cor" (Pearson correlation), 
"cov" (covariation), or "spearmanCor" (Spearman correlation).}

\item{nPerm}{Numeric. The number of permutations to do.}

\item{useSimpleCache}{Logical. Whether to use save caches. Caches
will be created for each permutation so that if the function is disrupted
it can restart where it left off. The final results are also saved 
as a cache. See simpleCache package for more details.}

\item{cacheDir}{Character. The path for the directory in which the
caches should be saved.}

\item{dataID}{Character. A unique identifier for this dataset 
(for saving results with simpleCache)}

\item{testType}{Character. Parameter for `getPermStat`. Whether to
create p values based on one a two sided test or a lesser/greater one
sided test. Options are: "greater", "lesser", "two-sided"}

\item{gammaFitMethod}{Character. method to use for fitting the gamma
distribution to null distribution. Options are 
"mme" (moment matching estimation), "mle" (maximum likelihood estimation), 
"qme" (quantile matching estimation), and "mge" (maximum goodness-of-fit 
estimation). See ?COCOA::getGammaPVal and 
?fitdistrplus::fitdist() for more info.}

\item{realScoreInDist}{Logical. Should the actual score (from 
test with no permutations) be included in the null distribution 
when fitting the gamma distribution. realScoreInDist=TRUE is 
recommended.}

\item{force}{Logical. If force=TRUE, when fitting the gamma distribution
returns an error (as may happen when a method other than "mme"
is used) then allow the error. If force=FALSE, when fitting the 
gamma distribution returns an error then don't return an error but 
instead use the "mme" method
for fitting that specific gamma distribution.}

\item{verbose}{A "logical" object. Whether progress 
of the function should be shown. One
bar indicates the region set is completed.}

\item{returnCovInfo}{logical. If TRUE, the following coverage and 
region set info will
be calculated and included in function output: regionSetCoverage, 
signalCoverage, totalRegionNumber, and meanRegionSize. For the
proportionWeightedMean scoring method, 
sumProportionOverlap will also be calculated.}

\item{...}{Character. Optional additional arguments for simpleCache.}
}
\value{
Returns a list with the following 4 items: 1. a list of length nPerm
where each item is a data.frame of the COCOA scores from a single 
permutation. Each data.frame is the output of `runCOCOA()` 
2. a data.table/data.frame of empirical p-values (the
output of `getPermStat`) 3. a 
data.table/data.frame of z-scores (the output of `getPermStat`. 
4. a data.frame of p-values based on
the gamma approximation (the output of getGammaPVal().
}
\description{
This is a convenience function that runs multiple steps of the 
permutation process together: it runs COCOA permutations, converts these
to null distributions, gets the empirical p value (which is limited by the
number of permutations), gets z scores, and fits a gamma distribution 
to each null distribution to estimate p values (not limited by the 
number of permutations),
Requires that the user has previously calculated the real COCOA scores. 
See these individual functions for more info on each step: runCOCOA, 
convertToFromNullDist, getPermStat, and getGammaPVal.
}
\details{
For reproducibility, set seed with 'set.seed()' function before running.
}
\examples{
data("esr1_chr1")
data("nrf1_chr1")
data("brcaMethylData1")
data("brcaMCoord1")
pcScores <- prcomp(t(brcaMethylData1))$x
targetVarCols <- c("PC1", "PC2")
targetVar <- pcScores[, targetVarCols]

# give the actual order of samples to `runCOCOA` to get the real scores
correctSampleOrder=1:nrow(targetVar)
realRSScores <- runCOCOA(genomicSignal=brcaMethylData1,
                        signalCoord=brcaMCoord1,
                        GRList=GRangesList(esr1_chr1, nrf1_chr1),
                        signalCol=c("PC1", "PC2"),
                        targetVar=targetVar,
                        sampleOrder=correctSampleOrder,
                        variationMetric="cor")
        
# give random order of samples to get random COCOA scores 
# so you start building a null distribution for each region set 
# (see vignette for example of building a null distribution with `runCOCOA`)
randomOrder <- sample(1:nrow(targetVar), 
                      size=nrow(targetVar),
                      replace=FALSE)
randomRSScores <- runCOCOA(genomicSignal=brcaMethylData1,
                          signalCoord=brcaMCoord1,
                          GRList=GRangesList(esr1_chr1, nrf1_chr1),
                          signalCol=c("PC1", "PC2"),
                          targetVar=targetVar,
                          sampleOrder=randomOrder,
                          variationMetric="cor")

# runCOCOAPerm
permResults <- runCOCOAPerm(genomicSignal=brcaMethylData1,
                           signalCoord=brcaMCoord1,
                           GRList=GRangesList(esr1_chr1, nrf1_chr1),
                           rsScores=realRSScores,
                           targetVar=targetVar,
                           signalCol=c("PC1", "PC2"),
                           variationMetric="cor",
                           nPerm = 10,
                           useSimpleCache=FALSE)
permResults
  

}
