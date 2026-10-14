#!/bin/bash
# Download EcoMAC (Carrera et al., Mol Syst Biol 2014, doi:10.15252/msb.20145108):
# 2,198 RMA-normalized E. coli Affymetrix arrays x 4,189 genes, with per-array
# metadata, as redistributed in the CovertLab wcEcoli repository. Then look up
# each array's Affymetrix platform in GEO (series platform; per-sample
# platform for series that mix platforms). Writes data/ecomac/.
# Usage (repo root):  tools/get_ecomac.sh
set -u
cd "$(dirname "$0")/.." || exit 1
D=data/ecomac; mkdir -p "$D"
B=https://raw.githubusercontent.com/CovertLab/wcEcoli/master/reconstruction/ecoli/scripts/nca/data/compendium
for f in README.md samples.tsv gene_names.tsv gene_synonyms.tsv experimental_shift.tsv seq/EcoMAC.tsv; do
  [ -s "$D/$(basename $f)" ] || curl -sfL -o "$D/$(basename $f)" "$B/$f" || { echo "download failed: $f"; exit 1; }
done
geo() { curl -s -m 30 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$1&targ=self&form=text&view=brief" |
        grep "^!$2_platform_id" | cut -d' ' -f3 | paste -sd,; sleep 0.4; }
if [ ! -s "$D/gse_platform.tsv" ]; then
  for g in $(cut -f3 "$D/samples.tsv" | tail -n +2 | grep '^GSE' | sort -u); do
    printf '%s\t%s\n' "$g" "$(geo "$g" Series)"; done > "$D/gse_platform.tsv"
fi
if [ ! -s "$D/gsm_platform.tsv" ]; then
  mixed=$(awk -F'\t' '$2 ~ /,/ {print $1}' "$D/gse_platform.tsv" | paste -sd'|')
  for s in $(awk -F'\t' -v m="^($mixed)$" 'NR>1 && $3 ~ m {print $2}' "$D/samples.tsv" | sort -u); do
    printf '%s\t%s\n' "$s" "$(geo "$s" Sample)"; done > "$D/gsm_platform.tsv"
fi
echo "EcoMAC in $D: $(awk 'END{print NR}' "$D/EcoMAC.tsv") genes x $(head -1 "$D/EcoMAC.tsv" | awk -F'\t' '{print NF}') arrays"
