DATA="data/2026_07_09_5pm" 
#RES_FOLDER="results/2026_07_09_5pm"
F_IN_NAME="ebola-bdbv_aligned-ADJUSTED-MATCH-nuc_2026-07-09T1551.fasta"
F_OUT_NAME="ebola-bdbv_aligned-ADJUSTED-MATCH-REF-nuc_2026-07-09T1551.fasta"

# OLD
#mafft --thread 4 --addfragments "$DATA/$F_IN_NAME" "$DATA/refseq_bdbv_NC_014373.fasta" > "$DATA/$F_OUT_NAME" #2>/dev/null

mafft --thread 4 --keeplength --add "$DATA/refseq_bdbv_NC_014373.fasta" "$DATA/$F_IN_NAME" > "$DATA/$F_OUT_NAME"
echo "Raccoon quality control running:"
raccoon aln-qc "$DATA/$F_OUT_NAME" -d "$DATA/alignment_qc_results" --reference-id NC_014373