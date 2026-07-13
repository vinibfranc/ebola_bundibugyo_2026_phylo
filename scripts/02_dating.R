library(ape)
library(phytools)
library(glue)
library(dplyr)
library(treedater) # install.github
library(lubridate)
library(ggtree) # from bioconductor
library(ggplot2)
#devtools::install_github("xavierdidelot/BactDating")
#devtools::install_github("xavierdidelot/DiagnoDating")
library(DiagnoDating,quietly=T)

#10 seqs: 2026_05_27_6pm
#16 seqs: 2026_06_01_2pm

TR_DIR <- "results/2026_07_09_5pm/"
tr <-read.tree(glue("{TR_DIR}/bdbv.treefile"))
plot(tr, show.tip.label = T)
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

tr_rooted <- root(tr, outgroup = "NC_014373", resolve.root = T)
plot(tr_rooted, show.tip.label = T)

tr_drop_outgroup <- drop.tip(tr_rooted, tip="NC_014373")
plot(tr_drop_outgroup, show.tip.label = T)

#tr_drop_outgroup2 <- unroot(tr_drop_outgroup)

# tr
td <- dater(tr_drop_outgroup, sts, s=aln_len,  clock="strict", 
            maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)
plot(td, show.tip.label = T)

td$mean.rate #  0.0003737813 (0.0004985478 vs 0.004114445 vs 0.00172937 10 seqs)
td$adjusted.mean.rate # same
# TODO too far back in past! Tried prev with Uganda samples and similar result!
# probably related to rooting to minimise residuals! Tried it but still giving end of 2025 estim
date_decimal(td$timeOfMRCA) # after rooting in outgroup and removing it: 2026-01-28
#initially: 2025-12-16 ("2026-04-26 08:03:30 UTC" vs "2026-03-21 05:05:54 UTC")

rootToTipRegressionPlot(td)
# Root-to-tip mean rate: 0.000967085925693231 
# Root-to-tip p value: 1.63614298888181e-05 
# Root-to-tip R squared (variance explained): 0.0046616579901538 vs 0.139611227042518

# suggestion 1: similar to Rambaut's fixed rates analysis, but bounds to rate rather than testing only two
# adjusted lower bound from 0.0012 to 0.0009
td_fixed_rate <- dater(tr_drop_outgroup, sts, s=aln_len,  clock="strict", meanRateLimits = c(0.0009, 0.0019),
                       maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)
# Fixing range of rate based on prev Ebola outbreaks give similar to their
# 0.001283867 vs 0.001261439 if 12 as lower bound
# 0.001136068 vs 0.001069195 if 9 as lower bound

td_fixed_rate$mean.rate # 0.001136068 (0.001069195 vs 0.001823019 vs 0.001632624 10 seqs)
date_decimal(td_fixed_rate$timeOfMRCA) # 2026-04-18 (2026-03-10 vs 2026-04-10 vs 2026-03-17) 
# using tr: matches their estimates (08 Mar and 15 Mar)
# using tr_drop_outgroup: 2026-04-18

# suggestion 2: additive clock with meanRateLimits instead of strict clock
td_additive_fixed_rate <- dater(tr_drop_outgroup, sts, s=aln_len,  clock="additive", meanRateLimits = c(0.0009, 0.0019),
                                maxit=10000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)
#NOTE: The estimated coefficient of variation of clock rates is high (>1).
td_additive_fixed_rate$mean.rate # 0.0009 (0.0009342773 vs 0.0019 vs 0.001605298)
date_decimal(td_additive_fixed_rate$timeOfMRCA) # 2026-04-07 vs (2026-03-06 vs 2026-04-14 vs 2026-03-18) 

rootToTipRegressionPlot(td_additive_fixed_rate, show.tip.labels = T, textopts = list(cex = 0.5)) #0.1396
tail_prob <- 0.05 # 0.01 giving 0 outliers
outliers <- outlierTips(td_additive_fixed_rate, alpha=tail_prob)
hist(outliers$q)
to_rm <- outliers[ outliers$q < tail_prob ,]
nrow(to_rm)
# 15 outliers
# rm outliers
tr2 <- drop.tip(tr_drop_outgroup, rownames(to_rm))

# relaxedClockTest after outlierTips
rct1 <- relaxedClockTest(tr2, sts, aln_len, quiet=F, ncpu=4)
rct1
rct1$strict_treedater
# Time of common ancestor 
# 2026.03838295281 
# 
# Time to common ancestor (before most recent sample) 
# 0.435589649925532 
# 
# Weighted mean substitution rate (adjusted by branch lengths) 
# 0.000338218149901144 
# 
# Unadjusted mean substitution rate 
# 0.000338218149901144 
rct1$relaxed_treedater
# Time of common ancestor 
# 2026.0367439475 
# 
# Time to common ancestor (before most recent sample) 
# 0.437228655238414 
# 
# Weighted mean substitution rate (adjusted by branch lengths) 
# 0.000345697301235511 
# 
# Unadjusted mean substitution rate 
# 0.000346071709937856 
# 
# Clock model  
# uncorrelated 
# 
# Coefficient of variation of rates 
# 0.000805696293406743
rct1$clock
# strict

rct2 <- relaxedClockTest(tr2, sts, aln_len, meanRateLimits = c(0.0009, 0.0019), quiet=F, ncpu=4)
rct2$strict_treedater
rct2$relaxed_treedater
rct2$clock
# "uncorrelated"
# so result of relaxedClockTest depends if specify meanRateLimits

# I guess probably ADAR sites not being masked and maybe over-masking at ends might be
# why tmrca so diff compared to Rambaut's analysis when not specifying rate limits?
# if decide to check later: https://github.com/artic-network/raccoon#typical-workflow

# also using diagnoDating (in paper Erik/Xavier use outlierLineages and plotLikBranches). Run here to understand https://github.com/xavierdidelot/DiagnoDating/blob/main/reproducibility/outlier.R 
# Not using diagnoDating to remove outliers so far... 
# Not sure how to get the units in the plotLikBranches right, i.e., same units as bactDating.

# re-estimate tree
td_additive_fixed_rate2 <- dater(tr2, sts, s=aln_len,  clock="additive", meanRateLimits = c(0.0009, 0.0019),
                                maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)
td_additive_fixed_rate2$mean.rate # 0.0009547757 (0.0009668545 vs 0.0019 vs 0.001605298)
date_decimal(td_additive_fixed_rate2$timeOfMRCA) # 2026-04-11 (vs 2026-03-03)
rootToTipRegressionPlot(td_additive_fixed_rate2) # 0.0032
# Very similar to before

# remove the ones Rambaut flagged as problematic
# PP_0075ZAY not in md
to_rm3 <- c("PP_0075Z66.1", "PP_00764TQ.1", "PP_00765FE.1")
tr3 <- drop.tip(tr_drop_outgroup, to_rm3)
td_additive_fixed_rate3 <- dater(tr3, sts, s=aln_len,  clock="additive", meanRateLimits = c(0.0009, 0.0019),
                                 maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)

td_additive_fixed_rate3$mean.rate # 0.0009827949 (0.0009668545)
date_decimal(td_additive_fixed_rate3$timeOfMRCA) # 2026-04-10 (2026-03-03)
rootToTipRegressionPlot(td_additive_fixed_rate3, show.tip.labels = T, textopts = list(cex=0.5)) # 0.0144 (0.2173 vs 0.16358)

pb_additive_fixed_rate <- treedater::parboot(td_additive_fixed_rate2, nreps=100, ncpu=4, quiet=F, overrideTempConstraint=F, overrideSearchRoot=F)
pb_additive_fixed_rate$timeOfMRCA_CI
date_decimal(as.numeric(pb_additive_fixed_rate$timeOfMRCA_CI))
# 2026-02-13 to 2026-04-17 (2026-01-12 to 2026-03-30 vs 2026-03-04 to 2026-04-18 vs 2023-02-06 to 2026-04-14)
pb_additive_fixed_rate$meanRate_CI # 0.0006664047 0.0013679325 (0.0007018393 0.0012866412 vs 0.001316129 0.002742892 vs 0.001148926 0.002369057)

# plot tree (meanRateLimits disabled and strict clock, first version)
mrsd_decimal <- max(td_additive_fixed_rate2$sts)
mrsd_date <- as.Date(date_decimal(mrsd_decimal))
mrsd_str <- format(mrsd_date, "%Y-%m-%d")
# build node date lookup for all nodes
get_node_dates <- function(td_obj) {
 n_tips  <- length(td_obj$tip.label)
 n_nodes <- td_obj$Nnode
 n_total <- n_tips + n_nodes
 
 node_dates <- rep(NA_real_, n_total)
 
 # Tips: dates are known from sampling times
 # sts names match tip.label order
 node_dates[1:n_tips] <- td_obj$sts[td_obj$tip.label]
 
 # edge matrix: col1 = parent, col2 = child
 # edge.length[i] = time elapsed from parent to child
 # child_date = parent_date + edge.length  (time flows forward)
 # => parent_date = child_date - edge.length
 
 edge <- td_obj$edge
 el <- td_obj$edge.length
 
 # Iteratively resolve: once a child date is known,
 # compute its parent date. Repeat until all internal nodes resolved.
 # Use a queue / iterate until convergence (handles any traversal order)
 
 max_iter <- n_total * 2
 iter <- 0
 
 while (any(is.na(node_dates[node_ids <- (n_tips + 1):n_total])) &&
        iter < max_iter) {
  iter <- iter + 1
  for (i in seq_len(nrow(edge))) {
   parent <- edge[i, 1]
   child <- edge[i, 2]
   
   if (!is.na(node_dates[child]) && is.na(node_dates[parent])) {
    node_dates[parent] <- node_dates[child] - el[i]
   }
  }
 }
 
 return(node_dates)  # length = n_tips + n_nodes
}

all_node_dates <- get_node_dates(td_additive_fixed_rate2)
n_tips <- length(td_additive_fixed_rate2$tip.label)
node_ids <- (n_tips + 1):(n_tips + td_additive_fixed_rate2$Nnode)

# data frame to attach via %<+%
node_date_df <- data.frame(
 node = node_ids,
 node_label = format(as.Date(date_decimal(all_node_dates[node_ids])), "%Y-%m-%d")
)

p_labeled <- ggtree(td_additive_fixed_rate2, mrsd = mrsd_str) %<+% node_date_df +
 #theme_tree2() +
 geom_tiplab(size = 3, align = TRUE, linesize = 0.3) +
 # date label on every internal node
 geom_label(
  aes(label = node_label),
  size = 2.5,
  fill = "white",
  label.size = 0.2,          # border thickness
  na.rm = TRUE
 ) +
 #scale_x_continuous(
 # labels = function(x) format(as.Date(date_decimal(x)), "%Y-%m-%d"),
 # breaks = scales::pretty_breaks(n = 6),
 # expand = expansion(mult = c(0.05, 0.3))
 #) +
 coord_cartesian(clip = "off") +
 theme(plot.margin = margin(t = 5, r = 120, b = 10, l = 15)) +
 labs(title = "Treedater estimates for the EBOV Bundibugyo outbreak") #x = "Date", 

print(p_labeled)
ggsave(glue("{TR_DIR}/treedater_tree.png"), p_labeled, width=12, height=15, dpi=300)

library(mlesky)

range_tr <- mrsd_decimal - td_additive_fixed_rate2$timeOfMRCA # 0.30
range_tr_weeks <- ceiling(range_tr * 52.1775) # 16

mlsk <- mlskygrid(td_additive_fixed_rate2, sampleTimes=sts, res=NULL, tau=NULL, tau_lower=0.001, tau_upper=100, 
                  ncross=10, ncpu=4, quiet=F, model=2)
plot(mlsk, logy=T)
mlsk$res # 14 (8 vs 5 vs 4)
mlsk$tau # 0.2020633 (0.2128471 vs 0.154 vs 0.6074089)

mlsk_pboot <- mlesky::parboot(mlsk, nrep=1000, ncpu=4, dd=F)

png(glue("{TR_DIR}/mlesky.png"), width=10, height=8, units="in", res=300)
plot(mlsk_pboot, logy=T, ggplot=T)
dev.off()


# BKP
# IMPORTANT: trees bigger now, so not running parboot for all configs
# normalApproxTMRCA useful?
pb <- treedater::parboot(td, nreps=1000, ncpu=4, quiet=F,
                         overrideTempConstraint=F, overrideSearchRoot=F)
plot(pb)
pb$timeOfMRCA_CI
date_decimal(as.numeric(pb$timeOfMRCA_CI))
# 2026-03-24 to 2026-04-30 (vs 2023-02-09 to 2026-04-24)
pb$meanRate_CI # 0.001647883 0.010259537 (vs 0.0001876817 0.0159350741)

pb_fixed_rate <- treedater::parboot(td_fixed_rate, nreps=1000, ncpu=4, quiet=F,
                                    overrideTempConstraint=F, overrideSearchRoot=F)
pb_fixed_rate$timeOfMRCA_CI
date_decimal(as.numeric(pb_fixed_rate$timeOfMRCA_CI))
# xx (2026-03-09 to 2026-04-17 vs 2023-02-06 to 2026-04-11)
pb_fixed_rate$meanRate_CI # xx (0.001350011 0.002461758 vs 0.001208905 0.002204855)