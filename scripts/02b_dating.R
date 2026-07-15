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
rct <- relaxedClockTest(tr1, sts, aln_len, meanRateLimits = c(0.0005, 0.002), quiet=F, ncpu=8, overrideTempConstraint = FALSE, nreps = 1e3)

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


# bootstrap 
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


# parametric bootstrap 
# pbtd <- treedater::parboot( td1a, nreps=200, ncpu=8, quiet=F, overrideTempConstraint=F, overrideSearchRoot=F)


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
# NEEDS DiagnoDating >= 0.9.1, ie the fix-treedater-units branch (PR #1), which is where the
# additive clock, the seqlen argument and matching of named dates were fixed. Before that,
# treedater results always came back as model='poisson' with relax=0:
#   devtools::install_github('xavierdidelot/DiagnoDating', ref = 'fix-treedater-units')
roottotip( tr1, sts )

# DiagnoDating, like BactDating, wants branch lengths in substitutions rather than per site,
# so rescale the tree. seqlen then tells treedater how to convert back, and omega0 and
# meanRateLimits stay in treedater's own units of substitutions per site per year.
trdd <- tr1
trdd$edge.length <- trdd$edge.length * aln_len

ddtd <- runDating( trdd, sts[ tr1$tip.label ], algo = 'treedater', keepRoot = TRUE
	, seqlen = aln_len
	, clock = 'additive'
	, omega0 = td1a$mean.rate
	, meanRateLimits = td1a$mean.rate*c(1,1+1e-6) )
ddtd # rate is per genome per year, so comparable to mu from bd1a; relax comparable to sigma
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



# lik profile 
ratetologlik <- function(omega){
	print( omega )
	# dater(tr1, sts, s=aln_len,  clock="additive", omega0 = omega, meanRateLimits = omega*c(1,1+1e-1), maxit=1, quiet=TRUE, ncpu=4)$loglik
	td <- dater(tr1, sts, s=aln_len,  clock="additive", omega0 = omega, meanRateLimits = omega*c(.9,1.4), quiet=TRUE, ncpu=4)
	c( td$adjusted.mean.rate, td$loglik )
}
omegas <- seq( 0.0001, 0.002, length.out = 30)
# omegas <- seq( 0.0004, 0.0008, length.out = 30)
llprof <- sapply( omegas, ratetologlik )
omegas <- llprof[1,]
llprof <- llprof[2,]
# plot( omegas, llprof )
source( './profplot-02b.R')
prpl <- profplot( omegas, llprof )




# mlesky 
library(mlesky)
range_tr <- max(sts) - td1a$timeOfMRCA #
# [1] 0.2772165
range_tr_weeks <- ceiling(range_tr * 52.1775) # 
# [1] 15
mlsk <- mlskygrid(td1a, sampleTimes=sts, res=NULL, tau=NULL, tau_lower=0.001, tau_upper=100, 
                  ncross=10, ncpu=4, quiet=F, model=2)
plot(mlsk, logy=T)
mlsk$res #
# [1] 17
mlsk$tau 
# [1] 0.185064

mlsk_pboot <- mlesky::parboot(mlsk, nrep=1000, ncpu=8, dd=F)
# regular bootstrap NEEDS mlesky >= 0.1.9: #   devtools::install_github('emvolz-phylodynamics/mlesky')
# bmlsk <- mlesky::boot( mlsk, btd$trees, ncpu = 8, sampleTimes = sts )

# png(glue("{TR_DIR}/mlesky.png"), width=10, height=8, units="in", res=300)
plot(mlsk_pboot, logy=T, ggplot=T)
# dev.off()
#


############################################################################################
# Jazzy dated-tree figure: geography of tips + ancestral-locality reconstruction + TMRCA CI
#   - tips coloured by geoLocAdmin2 (top 5 localities + Other)
#   - branches & internal nodes painted by a parsimony reconstruction of locality
#     (node point size = reconstruction confidence; the root, being ambiguous, shows small)
#   - horizontal error bar at the root = 95% CI of the TMRCA from the treedater bootstrap (btd)
#   - dashed line at 2026.35 = the growth-rate change seen in the skygrid
# Needs phangorn (parsimony ASR); circular node pies would need ggimage (unavailable here),
# so used branch-painting + confidence-scaled node points to convey the same information.
# Produced twice: the full dated tree, and a version with unsupported (near-zero, < 0.5
# substitution) internal branches collapsed to polytomies, so the two can be compared side by side.
library(phangorn)

# Collapse internal branches supported by < min_subs substitutions into polytomies, keeping the
# treedater dates. Uses td$intree (substitutions, same edge order as the dated tree). NB re-dating
# a collapsed tree does not work -- dater() re-resolves polytomies -- so we collapse after dating.
# CAUTION on the threshold: di2multi collapses each branch independently, so a *staircase* of
# individually-short branches all collapses into one node. Here 86/124 internal branches are ~0
# (<0.1 sub) but the backbone rungs run up to ~1 sub; min_subs=0.5 removes the zero-length nodes
# and leaves a normal (degree-3) root, whereas min_subs>=1 also eats the ~1-sub root edge and
# chains the whole backbone onto the root (a spurious ~88-way basal polytomy). Keep it below 1.
collapse_dated <- function(td, aln_len, min_subs = 0.5) {
	phy  <- as.phylo(td)
	subs <- td$intree$edge.length * aln_len              # substitutions per edge (aligned to phy$edge)
	nd   <- node.depth.edgelength(phy)                   # node heights (time) from root
	tmp  <- phy; tmp$edge.length <- subs
	topo <- di2multi(tmp, tol = min_subs)                # collapse by substitutions
	Ntt  <- Ntip(topo); dts <- numeric(Ntt + topo$Nnode) # re-time from the original treedater dates
	dts[1:Ntt] <- (td$timeOfMRCA + nd[1:Ntip(phy)])[match(topo$tip.label, phy$tip.label)]
	desc <- Descendants(topo, (Ntt+1):(Ntt+topo$Nnode), "tips")
	for (i in seq_along(desc)) {
		tips <- topo$tip.label[desc[[i]]]
		mn   <- if (length(tips)==Ntip(phy)) Ntip(phy)+1 else getMRCA(phy, tips)
		dts[Ntt+i] <- td$timeOfMRCA + nd[mn]
	}
	topo$edge.length <- dts[topo$edge[,2]] - dts[topo$edge[,1]]
	topo$root.time   <- td$timeOfMRCA
	topo
}
.norm <- function(x){ x <- toupper(trimws(x)); x[x=="" | x=="NOT PROVIDED"] <- "UNKNOWN"
	dplyr::recode(x, "MONGWALU"="MONGBWALU","MUNGWALU"="MONGBWALU","NYAKUNDE"="NYANKUNDE") }
# shared locality set + palette (both trees have identical tips, so colours are comparable)
loc0 <- .norm(md_match$geoLocAdmin2[ match(as.phylo(td1a)$tip.label, md_match$accessionVersion) ])
top5 <- names(sort(table(loc0), decreasing=TRUE))[1:5]
LEV  <- c(top5, "Other")
pal  <- setNames(c("#0072B2","#D55E00","#009E73","#E69F00","#CC79A7","grey72"), LEV)
mrsd <- as.Date(date_decimal(max(sts))); ci <- btd$timeOfMRCA_CI; break_gr <- 2026.35
mb <- decimal_date(seq(as.Date("2026-01-01"), as.Date("2026-08-01"), by="month"))
ml <- format(seq(as.Date("2026-01-01"), as.Date("2026-08-01"), by="month"), "%b %Y")
# build the geography + parsimony-ASR figure for any dated tree (treedater object or phylo).
# Node points are sized by ASR confidence. On a collapsed (multifurcating) tree, UFBoot from tr1
# is also printed at internal nodes (skipped on the full binary tree, where it would be a wall of
# numbers). NB UFBoot is over-confident on near-zero branches, so read the printed values sceptically.
geo_tree_fig <- function(tree, root_date, title) {
	phy <- as.phylo(tree); nt <- Ntip(phy)
	loc <- .norm(md_match$geoLocAdmin2[ match(phy$tip.label, md_match$accessionVersion) ])
	loc[!loc %in% top5] <- "Other"
	dat <- phangorn::phyDat(as.matrix(data.frame(s=loc, row.names=phy$tip.label)), type="USER", levels=LEV)
	mat <- do.call(rbind, lapply(phangorn::ancestral.pars(phy, dat, "MPR"), function(z) z[1,])); colnames(mat) <- LEV
	map  <- LEV[apply(mat,1,which.max)]; conf <- apply(mat,1,max)
	ntrans <- sum(map[phy$edge[,1]] != map[phy$edge[,2]])
	lab <- rep(NA_character_, nt+phy$Nnode)               # UFBoot labels, only on collapsed trees
	if (!is.binary(phy)) {
		ufb <- as.numeric(sub(".*/", "", tr1$node.label))
		for (n in (nt+1):(nt+phy$Nnode)) {
			tips <- phy$tip.label[ unlist(phangorn::Descendants(phy, n, "tips")) ]
			mn <- if (length(tips)==nt) Ntip(tr1)+1 else getMRCA(tr1, tips)
			if (!is.null(mn) && !is.na(mn)) lab[n] <- as.character(round(ufb[mn - Ntip(tr1)]))
		}
	}
	nodedf <- data.frame(node=1:(nt+phy$Nnode), state=factor(map,levels=LEV), conf=conf, lab=lab)
	pg <- ggtree(tree, mrsd=mrsd, size=0.5) %<+% nodedf + aes(color=state)
	# ggtree rounds mrsd to whole days, shifting its x-axis; off re-aligns absolute-date annotations
	rootx <- pg$data$x[pg$data$node == nt+1]; rooty <- pg$data$y[pg$data$node == nt+1]
	off <- rootx - root_date; cih <- 0.4   # half-height of the thin TMRCA CI band, in tip-row units
	pg$layers <- c(annotate("rect", xmin=ci[1]+off, xmax=ci[2]+off, ymin=rooty-cih, ymax=rooty+cih, fill="grey30", alpha=0.45), pg$layers)  # thin TMRCA 95% CI band at the root, behind the tree
	pg +
		geom_vline(xintercept=break_gr+off, linetype="dashed", color="grey55", linewidth=0.4) +
		# annotate("label", x=break_gr+off, y=nt*0.965, label="growth-rate\nchange 2026.35", size=2.7,
		#          color="grey40", hjust=-0.04, vjust=1, lineheight=0.9, fill="white", alpha=0.7) +
		annotate("point", x=rootx, y=rooty, shape=18, size=2.8, color="grey20") +
		annotate("text", x=mean(ci)+off, y=rooty+cih+nt*0.02, label="TMRCA 95% CI", size=2.6, color="grey35", vjust=0) +
		geom_nodepoint(aes(color=state, size=conf), alpha=0.85) +
		geom_nodelab(aes(label=lab), geom="label", size=1.9, vjust=-0.6, color="grey20",
		             fill="white", label.size=0, label.padding=unit(0.04,"lines"), na.rm=TRUE) +
		geom_tippoint(aes(fill=state), shape=21, size=2, stroke=0.25, color="white") +
		scale_color_manual(values=pal, guide="none") +
		scale_fill_manual(values=pal, name="geoLocAdmin2") +
		scale_size_continuous(range=c(0.6,3), guide="none") +
		scale_x_continuous(breaks=mb+off, labels=ml, limits=c(min(ci[1]+off, rootx)-0.02, max(sts)+off+0.02)) +
		labs(title=title,
		     subtitle=glue("{nt} genomes | {phy$Nnode} internal nodes | MRCA {format(as.Date(date_decimal(root_date)),'%d %b %Y')}")) +
		theme_tree2() +
		theme(plot.title=element_text(face="bold",size=12), plot.subtitle=element_text(size=8.5,color="grey40"),
		      axis.text.x=element_text(size=9,color="grey40"), legend.position=c(0.12,0.75),
		      legend.background=element_rect(fill="white",color="grey85"), legend.key.size=unit(0.9,"lines"))
}
# (1) full-resolution dated tree
# pg_full <- geo_tree_fig(td1a, td1a$timeOfMRCA, "EBOV Bundibugyo: dated phylogeny (treedater, additive clock)")
# print(pg_full)
# ggsave(glue("{TR_DIR}/treedater_geo_tree.png"), pg_full, width=10, height=12, dpi=300)
# (2) same tree with near-zero internal branches (< 0.5 substitution) collapsed to polytomies
td1a_collapsed <- collapse_dated(td1a, aln_len, min_subs = 0.5)
pg_coll <- geo_tree_fig(td1a_collapsed, td1a_collapsed$root.time, "")
print(pg_coll)
# ggsave(glue("{TR_DIR}/treedater_geo_tree_collapsed.png"), pg_coll, width=10, height=12, dpi=300)


