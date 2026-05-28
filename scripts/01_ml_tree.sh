
DATA="data/2026_05_27_6pm"
RES_FOLDER="results/2026_05_27_6pm"
mkdir -p $RES_FOLDER
iqtree2 \
  -s "$DATA/ebola-bdbv_aligned-ADJUSTED-nuc_2026-05-27T1700.fasta" \
  -m MFP \
  -B 1000 \
  --bnni \
  -alrt 1000 \
  -T AUTO \
  --prefix "$RES_FOLDER/bdbv"
