# Keep in fasta file only focalAccessions that are listed in metadata tsv (seqset) file 
# (each accession separated by comma)

# Read adjusted alignment file
library(Biostrings)

F_NAME <- "ebola-bdbv_aligned-ADJUSTED-nuc_2026-07-09T1551.fasta"
SEQSET_NAME <- "BDBV_DRC_2026-07-07.tsv"

DATA_FOLDER <- "data/2026_07_09_5pm"
#RES_FOLDER <- "results/2026_07_09_5pm"

metadata <- read.delim(file.path(DATA_FOLDER, SEQSET_NAME), sep = "\t", stringsAsFactors = FALSE)
accessions <- unique(trimws(unlist(strsplit(metadata$focalAccessions, ","))))
fasta <- readDNAStringSet(file.path(DATA_FOLDER, F_NAME))
# Keep only sequences whose names match the accessions
fasta_filtered <- fasta[names(fasta) %in% accessions]
F_OUT_NAME <- "ebola-bdbv_aligned-ADJUSTED-MATCH-nuc_2026-07-09T1551.fasta"
writeXStringSet(fasta_filtered, file.path(DATA_FOLDER, F_OUT_NAME))
