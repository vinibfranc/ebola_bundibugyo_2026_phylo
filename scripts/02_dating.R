library(ape)
library(glue)
library(dplyr)
library(treedater) # install.github
library(lubridate)
library(ggtree) # from bioconductor
library(ggplot2)

# IMPORTANT: added 50 Ns for all seqs at beginning of alignment and 10 Ns for all at the end 
# before estimating tree

TR_DIR <- "results/2026_05_27_6pm/"
tr <-read.tree(glue("{TR_DIR}/bdbv.treefile"))
plot(tr, show.tip.label = T)
aln_len <- 18940 - 50 - 10

DATA_DIR <- "data/2026_05_27_6pm/"
md <- read.csv(glue("{DATA_DIR}/ebola-bdbv_metadata_2026-05-27T1700.tsv"), sep="\t")
colnames(md)

length(intersect(tr$tip.label, md$accessionVersion))
md_match <- md %>% filter(accessionVersion %in% tr$tip.label)
md_match <- md_match %>% slice(match(tr$tip.label, accessionVersion))
all(tr$tip.label == md_match$accessionVersion)
md_match$sampleCollectionDate <- as.Date(md_match$sampleCollectionDate)
md_match$num_date <- decimal_date(md_match$sampleCollectionDate)

sts <- md_match$num_date
names(sts) <- md_match$accessionVersion

# any relevant to change default?
# omega0, minblen, estimateSampleTimes_densities, meanRateLimits

td <- dater(tr, sts, s=aln_len,  clock="strict", 
            maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)
plot(td, show.tip.label = T)

td$mean.rate # 0.00172937
td$adjusted.mean.rate # same
date_decimal(td$timeOfMRCA) # "2026-03-21 05:05:54 UTC"

# normalApproxTMRCA useful?
pb <- treedater::parboot(td, nreps=1000, ncpu=4, quiet=F,
              overrideTempConstraint=F, overrideSearchRoot=F)
plot(pb)
pb$timeOfMRCA_CI
date_decimal(as.numeric(pb$timeOfMRCA_CI))
# 2023-02-09 to 2026-04-24
pb$meanRate_CI # 0.0001876817 0.0159350741

# suggestion 1: similar to Rambaut's fixed rates analysis, but bounds to rate rather than testing only two
td_fixed_rate <- dater(tr, sts, s=aln_len,  clock="strict", meanRateLimits = c(0.0012, 0.0019),
            maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)

td_fixed_rate$mean.rate # 0.001632624 (vs 0.00172937)
date_decimal(td_fixed_rate$timeOfMRCA) # 2026-03-17 (vs 2026-03-21) 

pb_fixed_rate <- treedater::parboot(td_fixed_rate, nreps=1000, ncpu=4, quiet=F,
                                    overrideTempConstraint=F, overrideSearchRoot=F)
pb_fixed_rate$timeOfMRCA_CI
date_decimal(as.numeric(pb_fixed_rate$timeOfMRCA_CI))
# 2023-02-06 to 2026-04-11 (vs 2023-02-09 to 2026-04-24)
pb_fixed_rate$meanRate_CI # 0.001208905 0.002204855 (vs 0.0001876817 0.0159350741)

# suggestion 2: additive clock with meanRateLimits instead of strict clock
td_additive_fixed_rate <- dater(tr, sts, s=aln_len,  clock="additive", meanRateLimits = c(0.0012, 0.0019), 
                                maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)
#NOTE: The estimated coefficient of variation of clock rates is high (>1).

td_additive_fixed_rate$mean.rate # 0.001649809 (vs 0.001632624)
date_decimal(td_additive_fixed_rate$timeOfMRCA) # 2026-03-18 (vs 2026-03-17) 

pb_additive_fixed_rate <- treedater::parboot(td_additive_fixed_rate, nreps=1000, ncpu=4, quiet=F, overrideTempConstraint=F, overrideSearchRoot=F)
pb_additive_fixed_rate$timeOfMRCA_CI
date_decimal(as.numeric(pb_additive_fixed_rate$timeOfMRCA_CI))
# 2023-02-06 to 2026-04-14 (vs 2023-02-06 to 2026-04-11)
pb_additive_fixed_rate$meanRate_CI # 0.001148926 0.002369057 (vs 0.001208905 0.002204855)

# plot tree (meanRateLimits disabled and strict clock, first version)
mrsd_decimal <- max(td$sts)
mrsd_date <- as.Date(date_decimal(mrsd_decimal))
mrsd_str <- format(mrsd_date, "%Y-%m-%d")
# build node date lookup for all nodes
all_node_dates <- get_node_dates(td)
n_tips <- length(td$tip.label)
node_ids <- (n_tips + 1):(n_tips + td$Nnode)

# data frame to attach via %<+%
node_date_df <- data.frame(
 node = node_ids,
 node_label = format(as.Date(date_decimal(all_node_dates[node_ids])), "%Y-%m-%d")
)

p_labeled <- ggtree(td, mrsd = mrsd_str) %<+% node_date_df +
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
ggsave(glue("{TR_DIR}/treedater_tree.png"), p_labeled, width=10, height=8, dpi=300)

library(mlesky)

range_tr <- mrsd_decimal - td$timeOfMRCA
range_tr_weeks <- ceiling(range_tr * 52.1775) # using res=NULL for now

mlsk <- mlskygrid(td, sampleTimes=sts, res=NULL, tau=NULL, tau_lower=0.001, tau_upper=100, 
                  ncross=10, ncpu=4, quiet=F, model=2)
plot(mlsk)
mlsk$res # 4
mlsk$tau # 0.6074089

mlsk_pboot <- mlesky::parboot(mlsk, nrep=1000, ncpu=4, dd=F)

png(glue("{TR_DIR}/mlesky.png"), width=10, height=8, units="in", res=300)
plot(mlsk_pboot)
dev.off()