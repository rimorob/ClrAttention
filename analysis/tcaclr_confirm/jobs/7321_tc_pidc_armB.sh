. jobs_pc/common.sh
# mem: 8
# cpus: 32
# after: 7301_tc_prep_armB
# Arm B: PIDC by tools/pidc_chunked.jl, one process, 32 threads (equivalence gates of the scoreboard competitor study checked first).
set -e
python3 - <<'PY'
import csv,sys
for f in ["gate_m3d_vs_saved","gate_sub300_vs_pidc_jl","gate_sub600_vs_pidc_jl"]:
    r=list(csv.DictReader(open(f"results/scoreboard/complex_competitors/pidc_gate/{f}.csv")))[0]
    if r["pass"]!="True": sys.exit("PIDC equivalence gate failed: "+f)
PY
env -i HOME="$HOME" PATH=/usr/bin:/bin "$HOME/.juliaup/bin/julia" +1.6 -t 32 tools/pidc_chunked.jl results/tcaclr_confirm/armB/expression.tsv results/tcaclr_confirm/armB/pidc.f64 results/tcaclr_confirm/armB/pidc.log
