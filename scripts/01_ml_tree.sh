# 2026_06_01_2pm, 2026_05_27_6pm
DATA="data/2026_07_09_5pm" 
RES_FOLDER="results/2026_07_09_5pm"
#2026-06-01T1246, 2026-05-27T1700
F_NAME="ebola-bdbv_aligned-ADJUSTED-MATCH-nuc_2026-07-09T1551.fasta"

mkdir -p $RES_FOLDER
iqtree2 \
  -s "$DATA/$F_NAME" \
  -m MFP \
  -B 1000 \
  --bnni \
  -alrt 1000 \
  -T AUTO \
  --prefix "$RES_FOLDER/bdbv"