# Does the post-2026.35 growth-rate slowdown in the skygrid coincide with a shift in the
# GEOGRAPHIC composition of sampling? Concentrated sampling of a local transmission cluster
# oversamples closely-related lineages -> excess coalescence that can mimic a slowdown.
# Crosstab sample date x geoLocAdmin2 at weekly and twice-weekly resolution.
suppressMessages({ library(ape); library(dplyr); library(lubridate); library(glue)
                   library(tidyr); library(ggplot2); library(patchwork) })

TR_DIR   <- "results/2026_07_09_5pm/"
DATA_DIR <- "data/2026_07_09_5pm/"
F_NAME   <- "ebola-bdbv_metadata_2026-07-09T1550.tsv"
BREAK    <- 2026.35   # growth-rate change (~2026-05-08)

# --- tree samples that drive the skygrid (tr1: rooted, outgroup + outliers dropped) ---
tr <- read.tree(glue("{TR_DIR}/bdbv.treefile"))
.rootanddrop <- function(t0) drop.tip(root(t0, outgroup="NC_014373", resolve.root=T), tip="NC_014373")
tr1 <- .rootanddrop(tr)

norm_admin2 <- function(x){ x <- toupper(trimws(x)); x[x=="" | x=="NOT PROVIDED"] <- "UNKNOWN"
	dplyr::recode(x, "MONGWALU"="MONGBWALU","MUNGWALU"="MONGBWALU","NYAKUNDE"="NYANKUNDE") }

md <- read.csv(glue("{DATA_DIR}/{F_NAME}"), sep="\t") %>%
	filter(accessionVersion %in% tr1$tip.label) %>%
	mutate(date    = as.Date(sampleCollectionDate),
	       admin2  = norm_admin2(geoLocAdmin2),
	       week    = floor_date(date, "week", week_start = 1),
	       # twice-weekly: split each Monday-week into Mon-Thu and Fri-Sun halves
	       halfwk  = floor_date(date, "week", week_start = 1) +
	                 ifelse(wday(date, week_start = 1) <= 4, 0, 4))

top5 <- md %>% count(admin2, sort=TRUE) %>% slice_head(n=5) %>% pull(admin2)
LEV  <- c(top5, "Other")
md$loc <- factor(ifelse(md$admin2 %in% top5, md$admin2, "Other"), levels = LEV)

cat(glue("Tree samples: {nrow(md)}   {min(md$date)} to {max(md$date)}   break {as.Date(date_decimal(BREAK))} (={BREAK})\n"))
cat(glue("Top localities: {paste(top5, collapse=', ')}; all others -> Other ({sum(md$loc=='Other')} samples, {n_distinct(md$admin2[md$loc=='Other'])} districts)\n\n"))

shannon <- function(v){ p<-v/sum(v); -sum(p*log(p)) }

# crosstab (rows = time bin) + per-bin concentration over the FULL admin2 set
xtab <- function(binvar){
	comp <- md %>% count(bin=.data[[binvar]], loc) %>%
		pivot_wider(names_from=loc, values_from=n, values_fill=0) %>% arrange(bin)
	conc <- md %>% group_by(bin=.data[[binvar]]) %>%
		summarise(dec=round(decimal_date(min(bin)),3), n=n(),
		          nDistrict=n_distinct(admin2),
		          topShare=round(max(table(admin2))/n(),2),
		          H=round(shannon(as.numeric(table(admin2))),2), .groups="drop")
	out <- comp %>% left_join(conc %>% select(bin,dec,nDistrict,topShare,H), by="bin") %>%
		mutate(cross=ifelse(dec>=BREAK,"|","")) %>%            # marks bins at/after the break
		relocate(dec, .after=bin) %>% relocate(cross, .after=bin)
	for (L in LEV) if (!L %in% names(out)) out[[L]] <- 0L
	out %>% select(bin, cross, dec, all_of(LEV), nDistrict, topShare, H)
}

cat("================ WEEKLY  (cross '|' = bin at/after 2026.35) ================\n")
print(as.data.frame(xtab("week")), row.names=FALSE)
cat("\n================ TWICE-WEEKLY (Mon-Thu / Fri-Sun halves) ================\n")
print(as.data.frame(xtab("halfwk")), row.names=FALSE)

# ---- figure at twice-weekly resolution --------------------------------------------------
pal <- setNames(c("#0072B2","#D55E00","#009E73","#E69F00","#CC79A7","grey72"), LEV)
break_date <- as.Date(date_decimal(BREAK)); ink<-"#39424E"; muted<-"#7A8794"
conc_h <- md %>% group_by(halfwk) %>% summarise(n=n(), H=shannon(as.numeric(table(admin2))), .groups="drop")

pA <- ggplot(md, aes(halfwk, fill=loc)) +
	geom_vline(xintercept=break_date, linetype="dashed", colour=muted, linewidth=0.4) +
	geom_bar(width=3, colour="white", linewidth=0.15) +
	scale_fill_manual(values=pal, name="geoLocAdmin2") +
	labs(title="Twice-weekly sample composition by locality",
	     subtitle=glue("Dashed: growth-rate change {BREAK} ({break_date})"), x=NULL, y="samples") +
	theme_minimal(base_size=12) +
	theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
	      plot.title=element_text(face="bold",colour=ink), plot.subtitle=element_text(colour=muted,size=10),
	      axis.text=element_text(colour=muted), axis.title=element_text(colour=muted,size=10))
pB <- ggplot(conc_h, aes(halfwk, H)) +
	geom_vline(xintercept=break_date, linetype="dashed", colour=muted, linewidth=0.4) +
	geom_line(colour="#0072B2", linewidth=0.7) +
	geom_point(aes(size=n), shape=21, fill="#0072B2", colour="white", stroke=0.5) +
	scale_size_area(max_size=6, name="samples/bin") +
	labs(title="Locality diversity per bin (higher = less concentrated)",
	     x="sample collection (twice-weekly bins)", y="Shannon H") +
	theme_minimal(base_size=12) +
	theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold",colour=ink,size=12),
	      axis.text=element_text(colour=muted), axis.title=element_text(colour=muted,size=10))
ggsave(glue("{TR_DIR}/sampling_geography.png"), pA/pB + plot_layout(heights=c(1.3,1)), width=9, height=8, dpi=300)
cat(glue("\nSaved figure: {TR_DIR}sampling_geography.png\n"))
