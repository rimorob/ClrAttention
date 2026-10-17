#!/bin/bash
# Text progress bars for both job lanes (jobs_pc, jobs_pc2): running job,
# queue with cumulative ETA, memory and load.
# A job script may declare, in comment lines:
#   # expect: <minutes>               expected run time (time-based bar, marked ~)
#   # progress: <total> <regex>       count log lines matching regex (log-based bar)
# Usage (repo root):  tools/progress.sh        or   watch -n 30 tools/progress.sh
cd "$(dirname "$0")/.." || exit 1
now=$(date +%s)
bar() {  # fraction(0..1) width
  awk -v f="$1" -v w="$2" 'BEGIN { if (f > 1) f = 1; n = int(f * w + 0.5); s = "";
    for (i = 0; i < w; i++) s = s (i < n ? "#" : "."); printf "[%s] %3d%%", s, 100 * f }'
}
hm() { local m=$1; [ "$m" -lt 0 ] && m=0; printf "%dh%02dm" $((m / 60)) $((m % 60)); }
# expect/progress from the job script, else from <root>/estimates.tsv (job<TAB>minutes),
# which also works for scripts that are already running (never edit those)
field() { local v; v=$(grep -m1 -E "^# $2:" "$1" 2>/dev/null | sed -E "s/^# $2: *//")
  if [ -z "$v" ] && [ "$2" = expect ]; then
    v=$(awk -F'\t' -v j="$(basename "$1" .sh)" '$1 == j { print $2 }' "$(dirname "$(dirname "$1")")/estimates.tsv" 2>/dev/null); fi
  echo "$v"; }

echo "ClrAttention jobs  $(date '+%a %F %H:%M')   load $(cut -d' ' -f1 /proc/loadavg)   mem $(free -g | awk '/^Mem/ {print $3 " used / " $7 " avail GB"}')"
for root in jobs_pc jobs_pc2; do
  unit=$([ $root = jobs_pc ] && echo clr-jobs || echo clr-jobs2)
  mem=$(systemctl --user show $unit -p MemoryCurrent --value 2>/dev/null)
  [[ "$mem" =~ ^[0-9]+$ ]] && mem="$((mem / 1073741824)) GB" || mem="-"
  echo
  echo "== $root ($unit: $(systemctl --user is-active $unit 2>/dev/null), $mem)"
  t_free=$now
  for f in "$root"/running/*.sh; do
    [ -e "$f" ] || { echo "   (nothing running)"; continue; }
    b=$(basename "$f" .sh)
    st=$(grep " START $b\$" "$root/runner.log" | tail -1 | cut -d' ' -f1,2)
    s=$(date -d "$st" +%s 2>/dev/null || echo $now); el=$(( (now - s) / 60 ))
    exp=$(field "$f" expect); prog=$(field "$f" progress)
    if [ -n "$prog" ]; then
      tot=${prog%% *}; re=${prog#* }
      cnt=$(grep -c -E "$re" "$root/done/$b.log" 2>/dev/null); cnt=${cnt:-0}
      fr=$(awk -v c="$cnt" -v t="$tot" 'BEGIN { print c / t }')
      eta=$(awk -v f="$fr" -v e="$el" 'BEGIN { print (f > 0.02) ? int(e / f - e) : -1 }')
      lab="$cnt/$tot steps"
    elif [ -n "$exp" ]; then
      fr=$(awk -v e="$el" -v x="$exp" 'BEGIN { f = e / x; print (f > 0.99) ? 0.99 : f }')
      eta=$(( exp - el )); [ $eta -lt 5 ] && eta=5; lab="~time-based"
    else fr=0; eta=-1; lab="no estimate"; fi
    printf "   %-34s %s  %s elapsed  %s  %s\n" "$b" "$(bar "$fr" 30)" "$(hm $el)" \
      "$([ "$eta" -ge 0 ] && echo "ETA $(date -d "@$((now + eta * 60))" '+%a %H:%M')" || echo 'ETA ?')" "($lab)"
    [ "$eta" -ge 0 ] && t_free=$((now + eta * 60))
  done
  q=$(ls "$root"/queue/*.sh 2>/dev/null)
  [ -z "$q" ] && { echo "   queue empty"; continue; }
  echo "   queued:"
  unknown=0
  for f in $q; do
    b=$(basename "$f" .sh); exp=$(field "$f" expect)
    if [ -n "$exp" ] && [ $unknown -eq 0 ]; then t_free=$((t_free + exp * 60))
      printf "     %-34s ~%-7s done by ~%s\n" "$b" "$(hm $exp)" "$(date -d "@$t_free" '+%a %H:%M')"
    else [ -z "$exp" ] && unknown=1
      printf "     %-34s %s\n" "$b" "$([ -n "$exp" ] && echo "~$(hm $exp)" || echo 'no estimate')"; fi
  done
done
