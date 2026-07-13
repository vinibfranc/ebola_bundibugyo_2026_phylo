# Copying some setup code from 02_dating.R, might want to refactor later!
library(ape)
library(glue)
library(dplyr)
library(treedater) # install.github
library(lubridate)
library(ggtree) # from bioconductor
library(ggplot2)
#devtools::install_github("xavierdidelot/BactDating")
#devtools::install_github("xavierdidelot/DiagnoDating")
library(DiagnoDating,quietly=T)
library(DescTools)

# 2026_05_27_6pm, 2026_06_01_2pm
TR_DIR <- "results/2026_07_09_5pm/"
#tr <-read.tree(glue("{TR_DIR}/bdbv.treefile"))
plot(tr2, show.tip.label = T)
# 18940 - 10
# 18940 - 50 - 36
aln_len <- 18900 - 65 - 432 

# 2026_05_27_6pm, 2026_06_01_2pm
DATA_DIR <- "data/2026_07_09_5pm/" 
# ebola-bdbv_metadata_2026-05-27T1700.tsv, ebola-bdbv_metadata_2026-06-01T1246.tsv
F_NAME <- "ebola-bdbv_metadata_2026-07-09T1550.tsv"
md <- read.csv(glue("{DATA_DIR}/{F_NAME}"), sep="\t")

length(intersect(tr2$tip.label, md$accessionVersion)) #15
md_match <- md %>% filter(accessionVersion %in% tr2$tip.label)
md_match <- md_match %>% dplyr::slice(match(tr2$tip.label, accessionVersion))
all(tr2$tip.label == md_match$accessionVersion)
md_match$sampleCollectionDate <- as.Date(md_match$sampleCollectionDate)
md_match$num_date <- decimal_date(md_match$sampleCollectionDate)

sts <- md_match$num_date
names(sts) <- md_match$accessionVersion

# rates <- c(0.0012, 0.0019)
summary(tr2$edge.length)
sum(tr2$edge.length)
max(tr2$edge.length)
hist(tr2$edge.length, breaks = 50)
tr4 <- tr2
tr4$edge.length <- tr4$edge.length * aln_len
hist(tr4$edge.length, breaks = 50)

#bactdate_rate12 <- runDating(tree=tr4, dates=sts, algo = "BactDating", rate=0.0012*aln_len, keepRoot=F) #... ?bactdate shows all additional params
bactdate_rate10 <- runDating(tree=tr4, dates=sts, algo = "BactDating", rate=0.0010*aln_len, keepRoot=F)
#bactdate_rate10 <- bactdate(tree=tr4, date=sts, initMu=0.0010*aln_len, updateRoot=T)
# Result from BactDating, model poisson, clock rate 22.08 subst/year, relaxation parameter 0.00, root date 2026.12 (mid-Feb)
bactdate_rate19 <- runDating(tree=tr4, dates=sts, algo = "BactDating", rate=0.0019*aln_len, keepRoot=F)
#bactdate_rate19 <- bactdate(tree=tr4, date=sts, initMu=0.0019*aln_len, updateRoot=T)
# Result from BactDating, model poisson, clock rate 34.97, relaxation parameter 0.00, root date 2026.20 (mid-March)

# tr4
# treedater_additive <- runDating(tree=tr2, dates=sts, algo = "treedater", keepRoot=F,
#                                 clock="additive", meanRateLimits = c(0.0009, 0.0019),
#                                 maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)
treedater_additive <- dater(tre=tr2, sts=sts, s=aln_len, clock="additive", meanRateLimits = c(0.0009, 0.0019),
                                maxit=1000, searchRoot=5, numStartConditions=10, quiet=F, ncpu=4)
# Result from treedater, model poisson, clock rate 1.20, relaxation parameter 0.00, root date 2026.33 (1st May)

#attributes(bactdate_rate10)
#attributes(treedater_additive)

# function to run diagno_dating diags on all three trees above

run_diagno_diagnostics <- function(
  dated_trees, sts, aln_len,
  out_dir = "results/diagno_diagnostics",
  width = 1200, height = 900, res = 300, ppcheck_nrep = 1000) {
 
 # Sanity checks
 stopifnot(is.list(dated_trees))
 stopifnot(!is.null(names(dated_trees)))
 #stopifnot(all(sapply(dated_trees, inherits, "resDating")))
 stopifnot(is.numeric(sts) && !is.null(names(sts)))
 stopifnot(is.numeric(aln_len) && length(aln_len) == 1 && aln_len > 0)
 
 # Create output director
 if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message("Created output directory: ", out_dir)
 }
 
 # open PNG device
 open_png <- function(label, step_num, step_name) {
  fname <- file.path(
   out_dir,
   sprintf("%s_%02d_%s.png", label, step_num, step_name)
  )
  png(filename = fname, width = width, height = height, res = res)
  invisible(fname)
 }
 
 # does this resDating have an MCMC posterior?
 # BactDating stores MCMC iterations in r$record (a matrix).
 # treedater is a point-estimate method, so r$record is NULL.
 has_posterior <- function(r) {
  !is.null(r$record) && is.matrix(r$record) && nrow(r$record) > 1
 }
 
 # safe inputtree for roottotip
 # roottotip() needs branch lengths in NUMBER OF SUBSTITUTIONS.
 # BactDating inputtree has distances in substitutions (ready)
 # treedater inputtree has per-site lengths -> multiply by aln_len
 prepare_tree_for_roottotip <- function(r, sts, aln_len) {
  tr_in <- r$inputtree
  
  # Scale per-site → substitution counts ONLY for treedater.
  # BactDating inputtree is already in substitution counts
  # (tr4 = tr * aln_len was passed to runDating()).
  if (r$algo == "treedater") {
   tr_in$edge.length <- tr_in$edge.length #* aln_len
  }
  
  # Root the tree to maximise root-to-tip correlation before roottotip()
  # initRoot() is from BactDating (re-exported by DiagnoDating)
  tr_rooted <- tryCatch(
   initRoot(tr_in, sts),
   error = function(e) {
    message("initRoot() failed (", conditionMessage(e),
            "), falling back to unrooted tree.")
    tr_in
   }
  )
  tr_rooted
 }
 
 # Loop over each dated tree
 for (nm in names(dated_trees)) {
  
  if("treedater" %in% class(dated_trees[[nm]])) {
   r <- resDating(dated_trees[[nm]], tr2)
   print(attributes(r))
   print("a")
  } else {
   #r <- resDating(dated_trees[[nm]], tr4)
   r <- dated_trees[[nm]]
   print("b")
  }
  
  #r <- dated_trees[[nm]]
  is_bayesian <- has_posterior(r)
  
  message("\n========================================")
  message("Processing: ", nm,"  [algo = ", r$algo, " | posterior = ", is_bayesian, "]")
  message("========================================")
  
  # Step 1: Root-to-tip regression
  message("[1/6] Root-to-tip regression ...")
  fname1 <- open_png(nm, 1, "roottotip")
  tryCatch({
   tr_rtt <- prepare_tree_for_roottotip(r, sts, aln_len)
   roottotip(
    tree = tr_rtt, date = sts,
    showFig = TRUE, showTree = TRUE, permTest = 10000
   )
  }, error = function(e) {
   message("WARNING – roottotip() failed: ", conditionMessage(e))
   plot.new()
   title(main = paste("roottotip FAILED\n", conditionMessage(e)), col.main = "red")
  })
  dev.off()
  message("Saved: ", fname1)
  
  # Step 2: Plot the dated tree
  message("[2/6] Plotting dated tree ...")
  fname2 <- open_png(nm, 2, "dated_tree")
  tryCatch({
   # Do NOT modify r$tree — plot() uses it as-is
   # For treedater: r$tree branch lengths are in years (calendar time)
   # For BactDating: r$tree branch lengths are in substitutions
   # Both axes show the correct units for their respective method
   plot(r)
   
   # Add a subtitle clarifying the x-axis units per method
   units_label <- if (r$algo == "treedater") {
    "x-axis: decimal year (calendar time) | branch lengths in years"
   } else {
    "x-axis: decimal year | branch lengths in substitutions"
   }
   mtext(units_label, side = 1, line = 3, cex = 0.75, col = "grey40")
   
  }, error = function(e) {
   message("    WARNING – plot(r) failed: ", conditionMessage(e))
   plot.new()
   title(main = paste("plot(r) FAILED\n", conditionMessage(e)), col.main = "red")
  })
  dev.off()
  message("Saved: ", fname2)
  
  # Step 3: Posterior predictive check
  # Requires MCMC posterior, only available for BactDating results.
  message("[3/6] Posterior predictive check ...")
  fname3 <- open_png(nm, 3, "ppcheck")
  if (!is_bayesian) {
   message("SKIPPED – ppcheck() requires an MCMC posterior ",
           "(not available for algo = '", r$algo, "').")
   plot.new()
   title(
    main = paste0("ppcheck: SKIPPED\n", "algo '", r$algo, "' produces no MCMC posterior.\n",
                  "ppcheck() is only applicable to BactDating results."
    ), col.main = "darkorange", cex.main = 0.9
   )
  } else {
   tryCatch({
    ppcheck(
     x = r, nrep = ppcheck_nrep, showProgress = FALSE, showPlot = TRUE
    )
   }, error = function(e) {
    message("WARNING – ppcheck() failed: ", conditionMessage(e))
    plot.new()
    title(main = paste("ppcheck FAILED\n", conditionMessage(e)), col.main = "red")
   })
  }
  dev.off()
  message("    Saved: ", fname3)
  
  # Step 4: Likelihood of branches
  message("[4/6] Plotting branch likelihoods ...")
  fname4 <- open_png(nm, 4, "lik_branches")
  tryCatch({
   
   # Do NOT modify r — plotLikBranches() uses r$tree and r$resid as-is.
   # For treedater: x-axis (branch duration) is in years, y-axis rate
   # is in subst/site/year — these are the native treedater units.
   # For BactDating: x-axis is in years, y-axis rate is in subst/year.
   # Scaling would corrupt the internal consistency of the object.
   plotLikBranches(r)
   
   # Annotate with units so the reader knows what they are looking at
   rate_units <- if (r$algo == "treedater") {
    "Rate units: subst/site/year | Duration units: years"
   } else {
    "Rate units: subst/year | Duration units: years"
   }
   mtext(rate_units, side = 1, line = 0.5, outer = TRUE,
         cex = 0.75, col = "grey40")
   
  }, error = function(e) {
   message("WARNING – plotLikBranches() failed: ", conditionMessage(e))
   plot.new()
   title(main = paste("plotLikBranches FAILED\n", conditionMessage(e)),
         col.main = "red")
  })
  dev.off()
  message("Saved: ", fname4)
  
  
  # Step 5: Residuals
  message("[5/6] Plotting residuals ...")
  fname5 <- open_png(nm, 5, "residuals")
  tryCatch({
   plotResid(r)
  }, error = function(e) {
   message("WARNING – plotResid() failed: ", conditionMessage(e))
   plot.new()
   title(main = paste("plotResid FAILED\n", conditionMessage(e)),
         col.main = "red")
  })
  dev.off()
  message("Saved: ", fname5)
  
  # Step 6: Posterior distribution of p-values
  # Requires MCMC posterior, only available for BactDating results.
  message("[6/6] Posterior distribution of p-values ...")
  fname6 <- open_png(nm, 6, "postdist_pvals")
  if (!is_bayesian) {
   message("SKIPPED – postdistpvals() requires an MCMC posterior ",
           "(not available for algo = '", r$algo, "').")
   plot.new()
   title(
    main = paste0(
     "postdistpvals: SKIPPED\n",
     "algo '", r$algo, "' produces no MCMC posterior.\n",
     "postdistpvals() is only applicable to BactDating results."
    ), col.main = "darkorange", cex.main = 0.9
   )
  } else {
   tryCatch({
    postdistpvals(r, showPlot = TRUE)
   }, error = function(e) {
    message("    WARNING – postdistpvals() failed: ", conditionMessage(e))
    plot.new()
    title(main = paste("postdistpvals FAILED\n", conditionMessage(e)), col.main = "red")
   })
  }
  dev.off()
  message("Saved: ", fname6)
  
  message("Done with: ", nm)
 }
 
 message("\nAll diagnostics complete. Files written to: ", out_dir)
 invisible(out_dir)
}

dated_trees <- list(
 bactdate_rate10 = bactdate_rate10, bactdate_rate19 = bactdate_rate19,
 treedater_additive = treedater_additive
)

run_diagno_diagnostics(
 dated_trees = dated_trees, sts = sts, aln_len = aln_len,
 out_dir = "results/diagno_diagnostics",
 width = 1400, height = 1000, res = 150,
 ppcheck_nrep = 1000
)

# TODO issues:
# treedater_additive_02_dated_tree.png has axis in subst/site? (0 to 0.05) instead of numeric year
# treedater_additive_04_lik_branches.png also has substitutions in y-axis as subst/site/year (0 to 0.0025) 
# and x-axis in tree as subst per site while for bactdating those are in # substitutions and 
# numeric year units (2026.24 to 2026.36)