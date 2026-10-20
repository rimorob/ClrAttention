#!/bin/bash
# Resource-aware job scheduler (replaces the one-at-a-time tools/job_runner.sh).
# Runs several jobs from <root>/<qdir>/*.sh at once, in name order (= priority),
# while their declared needs fit the budget. A job declares, in comment lines:
#   # mem: <GB>          expected peak memory (default 30)
#   # cpus: <N>          threads it keeps busy (default 8)
#   # after: a,b         job names that must have ended with rc=0 first
# Each job runs in its own systemd user scope (clrjob-<name>.scope) with
# MemoryMax = 1.3 x mem (no swap), CPUWeight 10 and nice 19: a job that
# outgrows its cap is killed alone, never the lane or the desktop, and
# stopping the scope stops everything it started (R workers included):
#     systemctl --user stop clrjob-<name>.scope
# Admission: running declared mem + mem <= BUDGET_MEM, cpus + cpus <= BUDGET_CPU,
# and MemAvailable - mem >= MIN_AVAIL (GB). A job too big for the budget still
# runs when nothing else does. Later jobs may start ahead of an earlier one
# that does not fit yet (backfill), unless it is held by "# after:".
# Jobs started by the old runner (running/*.sh without a .res file) count with
# the defaults. Logs: <root>/runner.log (START/END lines, END with the
# measured peak, sampled every 15 s), <root>/done/<job>.log. A job file runs once unchanged for 30 s.
# Usage (repo root):  tools/scheduler.sh <root> [qdir]
#   e.g. systemd-run --user --unit=clr-sched -p CPUWeight=10 bash tools/scheduler.sh jobs_pc sched
cd "$(dirname "$0")/.." || exit 1
ROOT="${1:-jobs_pc}"; QDIR="$ROOT/${2:-queue}"
BUDGET_MEM=${BUDGET_MEM:-100}; BUDGET_CPU=${BUDGET_CPU:-50}; MIN_AVAIL=${MIN_AVAIL:-15}
DEF_MEM=30; DEF_CPU=8
export JOBS_ROOT="$ROOT"
mkdir -p "$QDIR" "$ROOT/running" "$ROOT/done" "$ROOT/control"
log() { echo "$(date '+%F %T') $*" >> "$ROOT/runner.log"; }
log "scheduler started (pid $$) on $(hostname), queue $QDIR, budget ${BUDGET_MEM} GB / ${BUDGET_CPU} threads, min available ${MIN_AVAIL} GB"
pgrep -f "tools/control_runner.sh $ROOT\$" >/dev/null || nohup tools/control_runner.sh "$ROOT" >/dev/null 2>&1 &

decl() { local v; v=$(grep -m1 -E "^# $2:" "$1" | sed -E "s/^# $2: *//; s/ .*//"); echo "${v:-$3}"; }
last_rc() { grep -E " END $1 rc=" "$ROOT/runner.log" | tail -1 | sed -E 's/.* rc=([0-9]+).*/\1/'; }
avail_gb() { awk '/^MemAvailable/ { print int($2 / 1048576) }' /proc/meminfo; }
declare -A held

start_job() {   # name mem cpus
  local b=$1 mem=$2 cpu=$3 cap
  cap=$(awk -v m="$mem" 'BEGIN { printf "%d", m * 1.3 + 0.999 }')
  mv "$QDIR/$b.sh" "$ROOT/running/$b.sh" || return
  echo "$mem $cpu" > "$ROOT/running/$b.res"; echo 0 > "$ROOT/running/$b.peak"
  log "START $b (mem ${mem} GB, cap ${cap} GB, cpus ${cpu})"
  systemctl --user reset-failed "clrjob-$b.scope" 2>/dev/null   # a failed scope of the same name blocks the unit name
  (
    systemd-run --user --scope --quiet --unit="clrjob-$b" --nice=19 \
      -p MemoryMax="${cap}G" -p MemorySwapMax=0 -p CPUWeight=10 -p IOWeight=10 \
      bash "$ROOT/running/$b.sh" > "$ROOT/done/$b.log" 2>&1
    rc=$?
    peak=$(awk '{ printf "%.1f", $1 / 1024 }' "$ROOT/running/$b.peak" 2>/dev/null)
    systemctl --user reset-failed "clrjob-$b.scope" 2>/dev/null
    mv "$ROOT/running/$b.sh" "$ROOT/done/$b.sh"; rm -f "$ROOT/running/$b.res" "$ROOT/running/$b.peak"
    log "END $b rc=$rc peak=${peak:-?}GB declared=${mem}GB"
  ) &
}

while true; do
  # memory accounting for running jobs: declared needs, and measured peaks
  used_mem=0; used_cpu=0; nrun=0
  for f in "$ROOT"/running/*.sh; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .sh); nrun=$((nrun + 1))
    if [ -f "$ROOT/running/$b.res" ]; then read m c < "$ROOT/running/$b.res"
      cur=$(systemctl --user show "clrjob-$b.scope" -p MemoryCurrent --value 2>/dev/null)
      if [[ "$cur" =~ ^[0-9]+$ ]]; then g=$((cur / 1048576)); p=$(cat "$ROOT/running/$b.peak" 2>/dev/null || echo 0)   # MB
        [ "$g" -gt "${p:-0}" ] && echo "$g" > "$ROOT/running/$b.peak"; fi
    else m=$DEF_MEM; c=$DEF_CPU; fi            # started by the old runner
    used_mem=$((used_mem + m)); used_cpu=$((used_cpu + c))
  done
  # admission in name order
  for f in "$QDIR"/*.sh; do
    [ -e "$f" ] || continue
    [ $(( $(date +%s) - $(stat -c %Y "$f") )) -ge 30 ] || continue
    b=$(basename "$f" .sh)
    mem=$(decl "$f" mem $DEF_MEM); cpu=$(decl "$f" cpus $DEF_CPU); deps=$(decl "$f" after "")
    ok=1
    for d in ${deps//,/ }; do
      rc=$(last_rc "$d")
      if [ -z "$rc" ]; then ok=0
      elif [ "$rc" != 0 ]; then ok=0
        [ -z "${held[$b]}" ] && { log "HOLD $b: dependency $d ended rc=$rc"; held[$b]=1; }
      fi
    done
    [ $ok = 1 ] || continue
    if [ $nrun -eq 0 ] || { [ $((used_mem + mem)) -le $BUDGET_MEM ] && [ $((used_cpu + cpu)) -le $BUDGET_CPU ] && \
                            [ $(( $(avail_gb) - mem )) -ge $MIN_AVAIL ]; }; then
      start_job "$b" "$mem" "$cpu"
      used_mem=$((used_mem + mem)); used_cpu=$((used_cpu + cpu)); nrun=$((nrun + 1))
      sleep 2
    fi
  done
  date '+%F %T' > "$ROOT/heartbeat"
  sleep 15
done
