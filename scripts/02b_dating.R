library(ape)
library(phytools)
library(glue)
library(dplyr)
library(treedater) # install.github
library(lubridate)
library(ggtree) # from bioconductor
library(ggplot2)
library(DiagnoDating,quietly=T)

TR_DIR <- "results/2026_07_09_5pm/"
tr <-read.tree(glue("{TR_DIR}/bdbv.treefile"))
# 18940 - 10
# 18940 - 50 - 36
aln_len <- 18900 - 65 - 432 

# 2026_05_27_6pm, 2026_06_01_2pm
DATA_DIR <- "data/2026_07_09_5pm/"
# ebola-bdbv_metadata_2026-05-27T1700.tsv, ebola-bdbv_metadata_2026-06-01T1246.tsv
F_NAME <- "ebola-bdbv_metadata_2026-07-09T1550.tsv"
md <- read.csv(glue("{DATA_DIR}/{F_NAME}"), sep="\t")
colnames(md)
table(md$geoLocAdmin1, md$geoLocAdmin2)

length(intersect(tr$tip.label, md$accessionVersion)) #148
md_match <- md %>% filter(accessionVersion %in% tr$tip.label)
md_match <- md_match %>% dplyr::slice(match(tr$tip.label, accessionVersion))
all(tr$tip.label == md_match$accessionVersion)
md_match$sampleCollectionDate <- as.Date(md_match$sampleCollectionDate)
md_match$num_date <- decimal_date(md_match$sampleCollectionDate)

sts <- md_match$num_date
names(sts) <- md_match$accessionVersion

.rootanddrop <- function( tr0 ){
	tr_rooted <- root(tr0, outgroup = "NC_014373", resolve.root = T)
	drop.tip(tr_rooted, tip="NC_014373")
}

tr <- .rootanddrop( tr )

# load and root bootstrap trees 
nboot = 200 
btres <- read.tree( glue("{TR_DIR}/bdbv.ufboot") )
btres <- btres[ sample.int(length(btres), size=nboot, replace=F )]
btres <- lapply( btres, .rootanddrop)


# initial fit 
td0<- dater(tr, sts, s=aln_len,  clock="strict", meanRateLimits = c(0.0005, 0.002), numStartConditions = 20, quiet=F, ncpu=4)

# iterate to remove outliers 
.otips <- function(td)   {
	outlierTips( td ) -> ot 
	ot[ ot$q < .05 , 'taxon']
}
ot <- .otips( td0 ) 
td1 <- td0; tr1 <- tr 
while( length( ot ) > 0)  {
	tr1 <- drop.tip( tr1, ot )
	td1<- dater(tr1, sts, s=aln_len,  clock="strict", meanRateLimits = c(0.0005, 0.002), numStartConditions = 20, quiet=F, ncpu=4)
	ot <- .otips( td1 ) 
}
# any outliers found? 
print( setdiff( td0$tip.label, td1$tip.label))
# character(0)


# clock test 
rct <- relaxedClockTest(tr1, sts, aln_len, meanRateLimits = c(0.0009, 0.0019), quiet=F, ncpu=8, overrideTempConstraint = FALSE)

td1u <- dater(tr1, sts, s=aln_len,  clock="uncorrelated", meanRateLimits = c(0.0005, 0.002), numStartConditions = 20, quiet=F, ncpu=4)
td1a <- dater(tr1, sts, s=aln_len,  clock="additive", meanRateLimits = c(0.0005, 0.002), numStartConditions = 20, quiet=F, ncpu=4)

c( td1u$loglik, td1a$loglik )
# [1] -249.7674 -215.9803


plot( td1a ); axisPhylo(root.time = td1a$timeOf, backward = F)
td1a
# 
# NOTE: The estimated coefficient of variation of clock rates is high (>1). Sometimes this indicates a poor fit and/or a problem with the data.
# 		
# 
# The following steps may help to fix or alleviate common problems:
# * Check that the vector of sample times is correctly named and that the units are correct.
# * If passing a rooted tree, make sure that the root position was chosen correctly, or estimate the root position by passing an unrooted tree (e.g. pass ape::unroot(tree))
# * The root position may be poorly estimated. Try increasing the _searchRoot_ parameter in order to test more lineages as potential root positions.
# * The model may be fitted by a relaxed or strict molecular clock. Try changing the _clock_ parameter
# * A poor fit may be due to a small number of lineages with unusual / outlying branch lengths which can occur due to sequencing error or poor alignment. Try the *outlierTips* command to identify and remove these lineages.
# * Check that there is adequate variance in sample times in order to estimate a molecular clock by doing a root-to-tip regression. Try the *rootToTipRegressionPlot* command. If the clock rate can not be reliably estimated, you can fix the value to a range using the _meanRateLimits_ option which would estimate a time tree given the previous estimate of clock rates.
# 
# Phylogenetic tree with 126 tips and 125 internal nodes.
# 
# Tip labels:
#   PP_006XHKB.2, PP_00764EH.1, PP_00765UM.1, PP_0076617.1, PP_00765M1.1, PP_0076SYS.1, ...
# Node labels:
#   , 0/72, 75/73, 0/21, 0/4, 0/7, ...
# 
# Rooted; includes branch length(s).
# 
#  Time of common ancestor 
# 2026.19675611258 
# 
#  Time to common ancestor (before most recent sample) 
# 0.277216490158025 
# 
#  Weighted mean substitution rate (adjusted by branch lengths) 
# 0.0006622018554641 
# 
#  Unadjusted mean substitution rate 
# 0.000627985792730883 
# 
#  Clock model  
# additive 
# 
#  Coefficient of variation of rates 
# 1.56980990630568 

date_decimal( td1a$timeOfMRCA  )
# [1] "2026-03-13 19:35:00 UTC"

rootToTipRegressionPlot( td1a) 
# Root-to-tip mean rate: 0.000243214384349191 
# Root-to-tip p value: 0.447464470497682 
# Root-to-tip R squared (variance explained): 0.00466165798993011 
# Returning fitted linear model.
# 
# Call:
# lm(formula = dG[1:ape::Ntip(td)] ~ sts)
# 
# Coefficients:
# (Intercept)          sts  
#  -0.4925812    0.0002432  
# 


btd <- treedater::boot( td1a, btres, ncpu = 8, overrideTempConstraint = FALSE, quiet=FALSE  )
btd
#                            pseudo ML        2.5 %       97.5 %
# Time of common ancestor 2.026197e+03 2.026016e+03 2.026239e+03
# Mean substitution rate  6.279858e-04 3.678732e-04 1.471493e-03
# 
#  For more detailed output, $trees provides a list of each fit to each simulation 
date_decimal( btd$timeOf )
#                      2.5%                     97.5% 
# "2026-01-06 18:48:49 UTC" "2026-03-29 04:06:48 UTC" 


# bact dating 
trbd1 <- tr1
trbd1$edge.length <- trbd1$edge.length * aln_len
bd1a <- bactdate( trbd1, sts[tr1$tip.label],  initMu = aln_len*td1a$mean.rate , model = 'arc' , showProgress = TRUE)
bd1a
# Phylogenetic tree dated using BactDating
# 
# Phylogenetic tree with 126 tips and 125 internal nodes.
# 
# Tip labels:
#   PP_006XHKB.2, PP_00764EH.1, PP_00765UM.1, PP_0076617.1, PP_00765M1.1, PP_0076SYS.1, ...
# Node labels:
#   , 0/72, 75/73, 0/21, 0/4, 0/7, ...
# 
# Rooted; includes branch length(s).
# Probability of root branch=0.54
# likelihood=-2.65e+02 [-2.79e+02;-2.52e+02]
# prior=-2.30e+02 [-2.60e+02;-2.04e+02]
# mu=1.26e+01 [9.44e+00;1.59e+01]
# sigma=3.54e-01 [6.54e-02;7.34e-01]
# alpha=2.33e+00 [1.76e+00;3.04e+00]
# Root date=2025.86 [2025.60;2026.03]
# Root date for most likely root=2025.85 [2025.58;2026.03]

#EV mean rate is okay, but the root is not plausible 


# diagnodating 
roottotip( tr1, sts )

runTreeDater2=function(tre,dates,keepRoot=F,seqlen=NA, ...) {
	stopifnot( !is.na(seqlen ))
	if (!keepRoot)  tre=unroot(tre)
	sts=dates
	if ( is.null( names( sts) )) names(sts)=tre$tip.label
	o=capture.output(rtd<-suppressWarnings(treedater::dater(tre,sts, s = seqlen,...)))
	model <- 'poisson'
	relax = 0
	if ( 'clock' %in% ...names() ) model <- list(...)[['clock']]
	if ( model == 'additive') {
		model <- 'arc'
		 relax = rtd$sp 
	}
	## rescale 
	tre$edge.length<- tre$edge.length * seqlen 
	res=resDating(rtd,tre,algo='treedater',model=model,rate=rtd$mean.rate*seqlen,relax=relax,rootdate=rtd$timeOfMRCA)
	# return(list( rtd, res) )
	return(res)
}
ddtd <- runTreeDater2(tr1, sts,  keepRoot = TRUE, seqlen = aln_len,  omega0 = td1a$mean.rate
	, meanRateLimits = td1a$mean.rate*c(1,1+1e-6)
	, clock = 'additive')
plotLikBranches( ddtd )
plotResid( ddtd )
testResid( ddtd )
# 
# 	Anderson-Darling test of goodness-of-fit
# 	Null hypothesis: Normal distribution
# 	with parameters mean = 0.000, sd = 1.000
# 
# data:  n
# An = 1.915, p-value = 0.1024
# 
