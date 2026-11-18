. jobs_pc/common.sh
# mem: 2
# cpus: 1
# Arm B expression for the long-method tools: header rename (orf to gene) of the frozen Kemmeren matrix; body checked byte-identical; hash of the original checked.
set -e

[ "$(sha256sum data/yeast/deleteome_all_1484.tsv | cut -d' ' -f1)" = "$(cut -d' ' -f1 docs/tcaclr_confirm_selections/armB_expression.sha256)" ]
mkdir -p results/tcaclr_confirm/armB
{ printf 'gene\t'; head -1 data/yeast/deleteome_all_1484.tsv | cut -f2-; tail -n +2 data/yeast/deleteome_all_1484.tsv; } > results/tcaclr_confirm/armB/expression.tsv
cmp <(tail -n +2 data/yeast/deleteome_all_1484.tsv) <(tail -n +2 results/tcaclr_confirm/armB/expression.tsv)
