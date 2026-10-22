#!/bin/bash
# Thermal and fan watch while jobs run (the two GPUs sit in adjacent slots and
# one card's fan can rub the other). Logs to docs/status/thermal.log:
# every 5 min for the first hour, then every 20 min, until --until (HH:MM).
#   WARN  (logged, then exit 2): CPU package >= 90 C; GPU >= 82 C; a processor
#         fan below 300 rpm with a CPU >= 70 C; a GPU fan at 0% with that GPU >= 65 C.
#   CRIT  (stop the scheduler and every job, then exit 3): CPU package >= 95 C or GPU >= 88 C.
# Usage (repo root):  systemd-run --user --unit=clr-thermal bash tools/thermal_watch.sh 16:00
cd "$(dirname "$0")/.." || exit 1
until_s=$(date -d "${1:-16:00}" +%s); start=$(date +%s)
log=docs/status/thermal.log; mkdir -p docs/status
note() { echo "$(date '+%F %T') $*" >> "$log"; }
pkg_temps() { for h in /sys/class/hwmon/hwmon*; do [ "$(cat $h/name)" = coretemp ] || continue
  for t in $h/temp*_input; do case "$(cat ${t%_input}_label 2>/dev/null)" in Package*) echo $(( $(cat $t) / 1000 ));; esac; done; done; }
dell=$(grep -l dell_smm /sys/class/hwmon/hwmon*/name | head -1 | xargs dirname)
note "thermal watch started, until $(date -d @$until_s '+%H:%M')"
while [ "$(date +%s)" -lt "$until_s" ]; do
  cpus=$(pkg_temps | paste -sd/); cmax=$(pkg_temps | sort -n | tail -1)
  fans=$(for f in $dell/fan*_input; do cat $f; done | paste -sd/)
  pfan_min=$(for f in $dell/fan*_input; do [ "$(cat ${f%_input}_label)" = "Processor Fan" ] && cat $f; done | sort -n | head -1)
  gpus=$(nvidia-smi --query-gpu=index,temperature.gpu,fan.speed,power.draw,utilization.gpu --format=csv,noheader,nounits 2>/dev/null)
  gline=$(echo "$gpus" | awk -F', ' '{ printf "GPU%s %sC fan %s%% %sW util %s%%; ", $1, $2, $3, int($4), $5 }')
  note "CPU pkg ${cpus}C | fans ${fans} rpm | ${gline}load $(cut -d' ' -f1 /proc/loadavg)"
  gmax=$(echo "$gpus" | awk -F', ' 'BEGIN { m = 0 } { if ($2 + 0 > m) m = $2 + 0 } END { print m }')
  gstall=$(echo "$gpus" | awk -F', ' '$3 + 0 == 0 && $2 + 0 >= 65 { print "GPU" $1 }')
  if [ "${cmax:-0}" -ge 95 ] || [ "${gmax:-0}" -ge 88 ]; then
    note "CRIT: CPU ${cmax}C / GPU ${gmax}C: stopping the scheduler and all jobs"
    systemctl --user stop clr-sched
    for u in $(systemctl --user list-units 'clrjob-*' --no-legend --plain | awk '{ print $1 }'); do systemctl --user stop "$u"; done
    exit 3
  fi
  w=""
  [ "${cmax:-0}" -ge 90 ] && w="$w CPU package ${cmax}C;"
  [ "${gmax:-0}" -ge 82 ] && w="$w GPU ${gmax}C;"
  [ -n "$pfan_min" ] && [ "$pfan_min" -lt 300 ] && [ "${cmax:-0}" -ge 70 ] && w="$w processor fan ${pfan_min} rpm at CPU ${cmax}C;"
  [ -n "$gstall" ] && w="$w fan at 0% on hot ${gstall};"
  if [ -n "$w" ]; then note "WARN:$w"; exit 2; fi
  [ $(( $(date +%s) - start )) -lt 3600 ] && sleep 300 || sleep 1200
done
note "thermal watch ended normally"
