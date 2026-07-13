library(cod)
library(ggnewscale)
library(ggplot2)
library(ggtree)

NCPU <- 4
PROF_CONTROL_LIST <- list(logtaulb = 0, logtauub = 50, res = 51, ncpu = NCPU)

cod_fit_bdbv <- codls( td_additive_fixed_rate2, logtau = NULL, 
                       profcontrol = PROF_CONTROL_LIST)

COD_DIR <- "results/cod"
if(!dir.exists(COD_DIR)) dir.create(COD_DIR)

pc <- plot(cod_fit_bdbv)
pc$layers[[length(pc$layers)]] <- NULL  # Remove tip labels
ggsave(file.path(COD_DIR, "bdbv_cod_fit_no_covar.png"), plot = pc, bg="white", dpi=300)

# Find optimal clustering threshold
chis_bdbv <- chindices(cod_fit_bdbv, clths = seq(0.1, 1.5, length = 20))
# Compute clusters
clusterdf_bdbv <- computeclusters(cod_fit_bdbv, clth = chis_bdbv$threshold[which.max(chis_bdbv$CH)])

md_match$sampleCollectionDate <- as.Date(md_match$sampleCollectionDate)
md_match$decimal_date <- as.numeric(decimal_date(md_match$sampleCollectionDate))
mrsd_s <- max(md_match$sampleCollectionDate)

# Plot clusters without tip labels
plotclusters_no_tips <- function(f, clusterdf, decimal_d, mrsd_, ...) {
 
 mrsd <- max(decimal_d)
 
 cl <- sort(unique(clusterdf$clusterid))
 cmat <- sapply(cl, function(x) (clusterdf$clusterid == x))
 colnames(cmat) <- cl
 rownames(cmat) <- clusterdf$tip.label
 
 # Debug: Check if tip labels match
 print("Checking tip label matching:")
 print(paste("Tree tips:", length(f$data$tip.label)))
 print(paste("Cluster tips:", length(clusterdf$tip.label)))
 print(paste("Matching tips:", sum(f$data$tip.label %in% clusterdf$tip.label)))
 
 # Only subset if tip labels match, otherwise reorder clusterdf to match tree
 if (all(f$data$tip.label %in% clusterdf$tip.label)) {
  cmat <- cmat[f$data$tip.label, , drop = FALSE]
 } else {
  # Alternative approach: reorder clusterdf to match tree tip order
  warning("Tip labels don't match perfectly. Attempting to reorder...")
  
  # Find common tip labels
  common_tips <- intersect(f$data$tip.label, clusterdf$tip.label)
  if (length(common_tips) == 0) {
   stop("No matching tip labels found between tree and cluster data")
  }
  
  # Subset both to common tips
  tree_order <- match(common_tips, f$data$tip.label)
  cluster_order <- match(common_tips, clusterdf$tip.label)
  
  cmat <- cmat[cluster_order, , drop = FALSE]
  rownames(cmat) <- common_tips
 }
 
 plot_no_tips <- function(f, mrsd) {
  f2beta <- f$coef
  tr1 <- f$data
  class(tr1) <- 'phylo'
  
  tr1$nodetimes <- mrsd - tr1$nodetimes
  # key to work!
  tr1$edge.length <- tr1$nodetimes[tr1$edge[,1]] - tr1$nodetimes[tr1$edge[,2]]
  #tr1$edge.length <- mrsd - tr1$nodetimes
  #print(paste("calendar_time time range:", min(calendar_time), "to", max(calendar_time)))
  print(paste("Node time range:", min(tr1$nodetimes), "to", max(tr1$nodetimes)))
  print(paste("MRSD:", mrsd))
  
  print("summaries")
  print(summary(tr1$edge.length))
  print(summary(tr1$nodetimes))
  
  print("heights")
  tree_height <- max(node.depth.edgelength(tr1))
  print(tree_height)
  node_height <- max(tr1$nodetimes)
  print(node_height)
  
  
  fdf <- data.frame(node = 1:length(tr1$nodetimes), theta = f2beta, date = tr1$nodetimes)
  #gtr1 <- ggtree(tr1, mrsd = as.Date(mrsd)) %<+% fdf
  gtr1 <- ggtree(tr1, mrsd = mrsd_) %<+% fdf #paste0(floor(mrsd), "-01-01")
  #gtr1 <- revts(gtr1)
  
  # Calculate breaks for x-axis
  date_range <- range(tr1$nodetimes)
  start_year <- floor(date_range[1] / 5) * 5 # Round down to nearest 5
  print(start_year)
  end_year <- ceiling(date_range[2] / 5) * 5 # Round up to nearest 5
  print(end_year)
  
  #tree_xlim <- range(mrsd - f$data$nodetimes)
  
  gtr1 <- gtr1 + aes(color = theta) + #x = date, 
   scale_x_continuous(breaks = seq(start_year, end_year, by = 5)) +
   scale_color_gradient2(low = 'blue', mid = 'lightblue', high = 'red', 
                         midpoint = 0, limits = range(fdf$theta), name = "psi") + 
   #coord_cartesian(xlim = c(tree_xlim[1], tree_xlim[2] + 1)) +
   theme_tree2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(plot = gtr1, breaks = seq(start_year, end_year, by = 5)))
 }
 
 plot_data <- plot_no_tips(f, mrsd)
 base_plot <- plot_data$plot
 #base_plot <- plot_no_tips(f, mrsd)
 
 # Additional debug info
 print(paste("cmat dimensions:", paste(dim(cmat), collapse = " x ")))
 #print("cmat row names (first 5):")
 #print(head(rownames(cmat), 5))
 print(paste("Number of clusters:", ncol(cmat)))
 #print("Cluster summary:")
 #print(table(clusterdf$clusterid))
 
 # Check if we have meaningful clusters
 if (ncol(cmat) == 1) {
  warning("Only one cluster found. This may indicate clustering threshold is too high.")
  # For single cluster, we can still create a basic plot
  result_plot <- base_plot + 
   ggtitle("Single cluster - all sequences grouped together")
 } else {
  # Try to create the heatmap with error handling
  tryCatch({
   result_plot <- gheatmap(base_plot, cmat, color=NULL, ...) + scale_x_ggtree() # key to work! + scale_y_continuous(expand=c(0, 0.3)) #offset = 0.02, width = 0.3
  }, error = function(e) {
   warning(paste("gheatmap failed:", e$message))
   # Return base plot with annotation
   result_plot <<- base_plot + 
    ggtitle(paste("Clustering completed -", ncol(cmat), "clusters found"))
  })
 }
 
 return(result_plot)
}

pc2 <- plotclusters_no_tips(cod_fit_bdbv, clusterdf_bdbv, decimal_d = md_match$decimal_date, mrsd_=mrsd_s)
ggsave(file.path(COD_DIR, "bdbv_cod_clusters.png"), plot = pc2, bg="white", dpi=300)

plot_cod_with_annotation <- function(cod_res, dataset, cluster_df = NULL, date_col, mrsd_ , 
                                     annotation_cols = c("geoLocAdmin2"), annotation_cats = c("BUNIA","RWAMPARA","NIZI", "OTHER"), logtau_label = NULL, offset = 0.3, width = 0.5) {
 
 mrsd <- max(as.numeric(dataset[[date_col]]), na.rm = TRUE)
 
 tree_tips <- cod_res$data$tip.label
 
 # psi values
 fdf <- data.frame(
  node = 1:length(cod_res$data$nodetimes), 
  theta = cod_res$coef
 )
 
 # Align dataset to tree
 hm_df <- dataset[match(tree_tips, dataset$tip.label), , drop = FALSE]
 
 # Add CH cluster assignments
 if (!is.null(cluster_df)) {
  clust_match <- cluster_df[match(tree_tips, cluster_df$tip.label), , drop = FALSE]
  hm_df$CH_cluster <- as.character(clust_match$clusterid)
  hm_df$CH_cluster[is.na(hm_df$CH_cluster)] <- "Unclustered"
  hm_df$CH_cluster <- factor(hm_df$CH_cluster)
 }
 
 #hm_df$
 
 print(annotation_cols)
 View(hm_df)
 
 # Fill NAs in annotation cols
 for (col in annotation_cols) {
  
  hm_df[[col]] <- factor(
   hm_df[[col]], levels = annotation_cats
  )
  
  hm_df[[col]][is.na(hm_df[[col]])] <- "OTHER"
  print(table(hm_df[[col]]) )
  
 }
 
 rownames(hm_df) <- tree_tips
 
 cod_res$data$nodetimes <- mrsd - cod_res$data$nodetimes
 cod_res$data$edge.length <- cod_res$data$nodetimes[cod_res$data$edge[,1]] - cod_res$data$nodetimes[cod_res$data$edge[,2]]
 
 date_range <- range(cod_res$data$nodetimes)
 start_year <- floor(date_range[1] / 5) * 5 # Round down to nearest 5
 print(start_year)
 end_year <- ceiling(date_range[2] / 5) * 5 # Round up to nearest 5
 print(end_year)
 
 # Base plot
 print("a")
 
 tr_plot <- ggtree(cod_res$data, mrsd = mrsd_) %<+% fdf #paste0(floor(mrsd), "-01-01")
 
 tr_plot <- tr_plot + aes(color = theta) + #x = date,
  scale_x_continuous(breaks = seq(start_year, end_year, by = 5)) +
  scale_color_gradient2(low = 'blue', mid = 'lightblue', high = 'red',
                        midpoint = 0, limits = range(fdf$theta), name = "psi") +
  #coord_cartesian(xlim = c(tree_xlim[1], tree_xlim[2] + 1)) +
  theme_tree2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
 #print(tr_plot)
 
 print("b")
 
 #col_width <- width / 2
 x_max <- max(tr_plot$data$x, na.rm = TRUE)  # rightmost tip in x-axis units
 x_range <- diff(range(tr_plot$data$x, na.rm = TRUE))
 # Absolute column width and gap in x-axis units (years)
 col_width_abs <- x_range * 0.04 # e.g. ~4 years wide per column
 gap_abs <- x_range * 0.001
 total_extra <- gap_abs + 3 * col_width_abs + 2 * gap_abs # NEW: 1 gap + 3 cols + 2 inter-col gaps
 total_width <- x_max + total_extra + x_range * 0.01 # NEW: matches xlim upper bound
 # gheatmap 'width' is relative to x_range, so convert:
 col_width_rel <- col_width_abs / x_range
 # Offsets: each is absolute distance from rightmost tip
 off1 <- gap_abs
 off2 <- off1 + col_width_abs + gap_abs
 
 hm1_mat <- hm_df[, annotation_cols[1], drop = FALSE]
 
 tr_plot <- gheatmap(
  tr_plot, hm1_mat,
  offset = off1, width = col_width_rel, #width / 2,
  colnames = F, color = NULL
 ) + #scale_fill_brewer( palette = "Set1", name = annotation_cols[1], na.value = "grey90")
  scale_fill_manual(
   values = c( BUNIA = "#E41A1C", RWAMPARA = "#377EB8", NIZI = "#4DAF4A", OTHER = "grey80"),
   drop = FALSE, name = annotation_cols[1]
  )
 #print(tr_plot)
 
 # NEW scale, then heatmap layer 2: CH clusters
 if (!is.null(cluster_df)) {
  hm2_mat <- hm_df[, "CH_cluster", drop = FALSE] # data.frame, NOT as.matrix()
  tr_plot <- tr_plot + new_scale_fill() # reset fill BEFORE gheatmap
  tr_plot <- gheatmap(
   tr_plot, hm2_mat,
   #offset = offset + (width / 2) + 0.05, # small gap between columns
   offset = off2, #offset + col_width,
   #width = width / 2, 
   width = col_width_rel, #col_width,
   colnames = F, color = NULL
  ) +
   scale_fill_brewer(palette = "Set2", name = "CH cluster", na.value = "grey90") +
   scale_x_ggtree() # restore x-axis ONCE at end
 }
 
 tr_plot$data$label <- ''
 
 if (!is.null(logtau_label)) {
  tr_plot <- tr_plot + ggtitle(logtau_label)
 }
 
 return(tr_plot)
}

md_match$geoLocAdmin2 <- toupper(md_match$geoLocAdmin2)
md_match_red <- md_match
md_match_red <- md_match_red %>% dplyr::select(accessionVersion, geoLocAdmin2, decimal_date) %>%
 dplyr::rename(tip.label = accessionVersion)
#rownames(md_match_red) <- md_match_red$accessionVersion

pcod_annot_bdbv <- plot_cod_with_annotation(
 cod_fit_bdbv, md_match_red, clusterdf_bdbv, date_col = "decimal_date",	mrsd_ = mrsd_s, 
 annotation_cols = c("geoLocAdmin2")
)
ggsave(file.path(COD_DIR, "bdbv_cod_fit_ANNOT.png"), plot = pcod_annot_bdbv, bg="white", width=10, height=8, dpi=300)
