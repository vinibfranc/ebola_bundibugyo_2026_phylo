# 2026_06_01_2pm, 2026_05_27_6pm
DATA="data/2026_07_09_5pm" 
RES_FOLDER="results/2026_07_09_5pm"
#2026-06-01T1246, 2026-05-27T1700
F_NAME="ebola-bdbv_aligned-ADJUSTED-MATCH-REF-nuc_2026-07-09T1551.fasta"
MD_FILE="$DATA/ebola-bdbv_metadata_2026-07-09T1550.tsv"

mkdir -p $RES_FOLDER
iqtree2 -s "$DATA/$F_NAME" -m MFP -B 1000 --bnni -alrt 1000 -T AUTO --prefix "$RES_FOLDER/bdbv"

echo "Creating intermediate file to make sure raccoon doesn't break"
#sed 's/[0-9]*\.[0-9]*\/\([0-9]*\)/\1/g' "$RES_FOLDER/bdbv.treefile" > "$RES_FOLDER/bdbv_clean.treefile"
sed 's/)[0-9.]*\/[0-9.]*/)/g' "$RES_FOLDER/bdbv.treefile" > "$RES_FOLDER/bdbv_clean.treefile"

raccoon tree-qc --tree "$RES_FOLDER/bdbv_clean.treefile" --alignment "$DATA/$F_NAME" \
  -d "$RES_FOLDER/tree_qc_results" --outgroup-ids "$DATA/outgroup_id.txt" \
  --run-adar --adar-window 300 --adar-min-count 3 --run-apobec
  
  
#--tip-fields $MD_FILE --tip-field-delimiter "\t" --tip-date-field sampleCollectionDate \
#--metadata-id-field accessionVersion --metadata-location-field geoLocAdmin2 --metadata-date-field sampleCollectionDate \
#--asr-state alignment.masked.fasta.state \ (needs --ancestral flag in iqtree above)