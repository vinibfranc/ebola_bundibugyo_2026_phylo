
## ---- profile likelihood plot -------------------------------------------------
profplot <- function( omegas, llprof){
	library(scales)
	prof <- data.frame( omega = omegas, ll = llprof )
	prof <- prof[ is.finite(prof$ll), ]
	# smooth in log10(rate) space -- the axis the profile is read on. The same fit
	# supplies the drawn curve AND the support interval, so band and line can't disagree.
	lfit <- loess( ll ~ log10(omega), data = prof, span = 0.6 )
	grid <- data.frame( omega = 10^seq( log10(min(prof$omega)), log10(max(prof$omega)),
                                    	length.out = 1000 ) )
	grid$ll <- predict( lfit, newdata = grid )
	llmax <- max( grid$ll, na.rm = TRUE )
	peak  <- grid$omega[ which.max(grid$ll) ]
	thr   <- llmax - 1.96                      # support interval: within 1.96 loglik units
	inside <- grid$omega[ which(grid$ll >= thr) ]
	lo <- min(inside); hi <- max(inside)
	# flag if the profile never falls 1.96 units below the max inside the grid
	open_lo <- lo <= min(prof$omega) * 1.0001
	open_hi <- hi >= max(prof$omega) * 0.9999
	if (open_lo || open_hi)
		warning( glue("Support interval is open at the {ifelse(open_lo,'lower','upper')} end -- widen the range of *omegas*.") )
	ink <- "#39424E"; muted <- "#7A8794"; accent <- "#2E6DB4"
	p_prof <- ggplot( prof, aes(omega, ll) ) +
		annotate( "rect", xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf,
	          	fill = accent, alpha = 0.10 ) +
		geom_hline( yintercept = thr,  colour = muted, linetype = "dashed", linewidth = 0.35 ) +
		geom_vline( xintercept = peak, colour = muted, linetype = "dotted", linewidth = 0.35 ) +
		geom_line( data = grid, aes(omega, ll), colour = accent, linewidth = 0.9 ) +
		geom_point( shape = 21, size = 2.1, stroke = 0.6, colour = muted, fill = "white" ) +
		scale_x_log10( breaks = c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3),
	               	labels = label_number( scale = 1e4, accuracy = 0.1 ) ) +
		labs(
			title    = "Profile likelihood for the clock rate (additive clock)",
			subtitle = glue("shaded: within 1.96 log-lik units of the max, {signif(lo,3)} to {signif(hi,3)}"),
			x        = expression(paste("Clock rate (", 10^-4, " subs/site/year, log scale)")),
			y        = "Log likelihood"
		) +
		theme_minimal( base_size = 12 ) +
		theme(
			panel.grid.minor = element_blank(),
			panel.grid.major = element_line( linewidth = 0.25, colour = "grey90" ),
			plot.title       = element_text( face = "bold", colour = ink ),
			plot.subtitle    = element_text( colour = muted, size = 10 ),
			axis.title       = element_text( colour = muted, size = 10 ),
			axis.text        = element_text( colour = muted )
		)
	# ggsave( glue("{TR_DIR}/loglik_profile_rate.png"), p_prof, width = 7.5, height = 5, dpi = 300 )
	p_prof 
}

