#!/bin/bash
# Resource-aware job scheduler (replaces the one-at-a-time tools/job_runner.sh).
# Standing project rule: keep the machine at ~70-80% utilization, reliably.
# Runs several jobs from <root>/<qdir>/*.sh at once, in name order (= priority),
# while their needs fit the budget. A job declares, in comment lines:
#   # mem: <GB>          expected peak memory (default 30; decimals round up)
#   # cpus: <N>          threads it keeps busy (default 8)
#   # after: a,b         job names that must have ended with rc=0 first
# Each job runs in its own systemd user scope clrjob-<name>.scope inside
# clrjobs.slice (~/.config/systemd/user/clrjobs.slice: MemoryMax 95G, no swap,
# so all jobs together can never take the desktop's memory). Per-job cap:
# MemoryMax = min(1.3 x mem, slice cap), no swap, OOMPolicy=kill (one dead
# worker ends the whole job), CPUWeight 10, nice 19. Stopping the scope stops
# everything the job started:   systemctl --user stop clrjob-<name>.scope
# Jobs survive a scheduler restart: the scope writes the exit code to
# running/<job>.rc, and every loop finalizes jobs whose scope has ended
# (END line with rc, systemd result and measured peak).
# Admission: sum over running jobs of max(declared, current) (after 30 min:
# max(1.1 x measured peak, current)) + mem <= BUDGET_MEM,
# declared cpus + cpus <= BUDGET_CPU (or 1-min load + cpus <= BUDGET_CPU), and MemAvailable - mem >= MIN_AVAIL (always).
# Backfill: later jobs may start ahead of an earlier one that does not fit,
# unless the earliest eligible job has waited > RESERVE_MIN minutes (then the
# queue waits for it). Jobs started by the old runner (running/*.sh without a
# .res file) count with their measured memory (clr-jobs service), at least the
# default. A job file runs once unchanged for 30 s (copy as .tmp, then mv).
# Usage (repo root):  tools/scheduler.sh <root> [qdir]
#   systemd-run --user --unit=clr-sched -p Restart=on-failure -p CPUWeight=10 bash tools/scheduler.sh jobs_pc sched
cd "$(dirname "$0")/.." || exit 1
ROOT="${1:-jobs_pc}"; QDIR="$ROOT/${2:-queue}"
# CPU over-commit only slows jobs; memory is the hard limit
BUDGET_MEM=${BUDGET_MEM:-90}; BUDGET_CPU=${BUDGET_CPU:-56}; MIN_AVAIL=${MIN_AVAIL:-15}
SLICE=${SLICE:-clrjobs.slice}; SLICE_CAP=${SLICE_CAP:-95}; RESERVE_MIN=${RESERVE_MIN:-30}
DEF_MEM=30; DEF_CPU=8
export JOBS_ROOT="$ROOT"
mkdir -p "$QDIR" "$ROOT/running" "$ROOT/done" "$ROOT/control"
log() { echo "$(date '+%F %T') $*" >> "$ROOT/runner.log"; }
log "scheduler started (pid $$) on $(hostname), queue $QDIR, slice $SLICE, budget ${BUDGET_MEM} GB / ${BUDGET_CPU} threads, min available ${MIN_AVAIL} GB"
pgrep -f "tools/control_runner.sh $ROOT\$" >/dev/null || nohup tools/control_runner.sh "$ROOT" >/dev/null 2>&1 &

declare -A warned eligible_since
num() {   # positive number -> integer (decimals round up); else the default (warned in the main loop)
  if [[ "$1" =~ ^[0-9]+([.][0-9]+)?$ ]]; then awk -v v="$1" 'BEGIN { i = int(v); print (v > i) ? i + 1 : i }'; else echo "$2"; fi
}
check() {  # warn once about an invalid declaration (runs in the main shell, so warned[] persists)
  [ -n "$1" ] && ! [[ "$1" =~ ^[0-9]+([.][0-9]+)?$ ]] && [ -z "${warned[$2$1]}" ] && { log "WARN $2: invalid value '$1', using the default"; warned[$2$1]=1; }
}
decl() { grep -m1 -E "^# $2:" "$1" | sed -E "s/^# $2: *//; s/[[:space:]].*//"; }
avail_gb() { awk '/^MemAvailable/ { print int($2 / 1048576) }' /proc/meminfo; }
scope_mem_gb() { local v; v=$(systemctl --user show "$1" -p MemoryCurrent --value 2>/dev/null)
  [[ "$v" =~ ^[0-9]+$ ]] && echo $(( (v + 1073741823) / 1073741824 )) || echo 0; }
# a dependency is satisfied when its newest END is rc=0 and it is not queued or running again
# (no comparison with the dependent's file time: a dependency may well have finished before
# the dependent was queued)
dep_state() {   # dep dependent_file -> ok | wait | failed | unknown
  local d=$1 f=$2 line ts
  [ -e "$QDIR/$d.sh" ] || [ -e "$ROOT/running/$d.sh" ] && { echo wait; return; }
  line=$(grep -F " END $d rc=" "$ROOT/runner.log" | tail -1)
  if [ -z "$line" ]; then [ -e "$ROOT/done/$d.sh" ] && echo wait || echo unknown; return; fi
  case "$line" in *" rc=0 "*|*" rc=0") echo ok ;; *) echo failed ;; esac
}

start_job() {   # name mem cpus
  local b=$1 mem=$2 cpu=$3 cap
  cap=$(awk -v m="$mem" -v s="$SLICE_CAP" 'BEGIN { c = int(m * 1.3 + 0.999); print (c > s) ? s : c }')
  mv "$QDIR/$b.sh" "$ROOT/running/$b.sh" || return
  rm -f "$ROOT/running/$b.rc"
  echo "$mem $cpu" > "$ROOT/running/$b.res"; echo 0 > "$ROOT/running/$b.peak"
  systemctl --user reset-failed "clrjob-$b.scope" 2>/dev/null   # a failed scope of the same name blocks the unit name
  log "START $b (mem ${mem} GB, cap ${cap} GB, cpus ${cpu})"
  setsid systemd-run --user --scope --quiet --unit="clrjob-$b" --slice="$SLICE" --nice=19 \
    -p MemoryMax="${cap}G" -p MemorySwapMax=0 -p OOMPolicy=kill -p CPUWeight=10 -p IOWeight=10 \
    bash -c 'bash "$1" > "$2" 2>&1; echo $? > "$3"' _ "$ROOT/running/$b.sh" "$ROOT/done/$b.log" "$ROOT/running/$b.rc" \
    >/dev/null 2>&1 < /dev/null &
  sleep 2
}
finalize() {    # job whose scope has ended
  local b=$1 rc res peak mem
  res=$(systemctl --user show "clrjob-$b.scope" -p Result --value 2>/dev/null)
  if [ -s "$ROOT/running/$b.rc" ]; then rc=$(cat "$ROOT/running/$b.rc")
  elif [ "$res" = success ]; then rc=0
  elif [ "$res" = oom-kill ]; then rc=137
  else rc=255; fi
  [ "$res" = oom-kill ] && [ "$rc" = 0 ] && rc=137
  peak=$(awk '{ printf "%.1f", $1 / 1024 }' "$ROOT/running/$b.peak" 2>/dev/null)
  read mem _ < "$ROOT/running/$b.res"
  systemctl --user reset-failed "clrjob-$b.scope" 2>/dev/null
  mv "$ROOT/running/$b.sh" "$ROOT/done/$b.sh"; rm -f "$ROOT/running/$b.res" "$ROOT/running/$b.peak" "$ROOT/running/$b.rc"
  log "END $b rc=$rc result=${res:-?} peak=${peak:-?}GB declared=${mem}GB"
}

while true; do
  used_mem=0; used_cpu=0; nrun=0
  for f in "$ROOT"/running/*.sh; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .sh)
    if [ -f "$ROOT/running/$b.res" ]; then
      # ended: the job wrote its exit code, or its scope is gone > 60 s after the launch
      # (killed, e.g. out of memory, or the scope never started)
      if [ -e "$ROOT/running/$b.rc" ]; then finalize "$b"; continue; fi
      state=$(systemctl --user is-active "clrjob-$b.scope" 2>/dev/null)
      if [ "$state" = failed ] || { [ "$state" != active ] && [ "$state" != activating ] && [ "$state" != deactivating ] && \
         [ $(( $(date +%s) - $(stat -c %Y "$ROOT/running/$b.res") )) -gt 60 ]; }; then finalize "$b"; continue; fi
      read m c < "$ROOT/running/$b.res"
      cur=$(systemctl --user show "clrjob-$b.scope" -p MemoryCurrent --value 2>/dev/null)
      if [[ "$cur" =~ ^[0-9]+$ ]]; then mb=$((cur / 1048576)); p=$(cat "$ROOT/running/$b.peak" 2>/dev/null || echo 0)
        [ "$mb" -gt "${p:-0}" ] && echo "$mb" > "$ROOT/running/$b.peak"; fi
      # after 30 min, count the job at 1.1 x its measured peak instead of the declaration
      if [ $(( $(date +%s) - $(stat -c %Y "$ROOT/running/$b.res") )) -gt 1800 ]; then
        m=$(awk '{ printf "%d", $1 * 1.1 / 1024 + 0.999 }' "$ROOT/running/$b.peak" 2>/dev/null); m=${m:-$DEF_MEM}; fi
      g=$(scope_mem_gb "clrjob-$b.scope"); [ "$g" -gt "$m" ] && m=$g
    else                                        # started by the old one-at-a-time runner
      g=$(scope_mem_gb clr-jobs.service); m=$DEF_MEM; [ "$g" -gt "$m" ] && m=$g; c=$DEF_CPU
    fi
    nrun=$((nrun + 1)); used_mem=$((used_mem + m)); used_cpu=$((used_cpu + c))
  done
  first_waiting=""
  for f in "$QDIR"/*.sh; do
    [ -e "$f" ] || continue
    [ $(( $(date +%s) - $(stat -c %Y "$f") )) -ge 30 ] || continue
    b=$(basename "$f" .sh)
    rm_=$(decl "$f" mem); rc_=$(decl "$f" cpus); check "$rm_" "$b"; check "$rc_" "$b"
    mem=$(num "$rm_" $DEF_MEM); cpu=$(num "$rc_" $DEF_CPU); deps=$(decl "$f" after)
    ok=1
    for d in ${deps//,/ }; do
      s=$(dep_state "$d" "$f")
      [ "$s" = ok ] && continue
      ok=0
      { [ "$s" = failed ] || [ "$s" = unknown ]; } && [ -z "${warned[$b$d$s]}" ] && { log "HOLD $b: dependency $d is $s"; warned[$b$d$s]=1; }
    done
    [ $ok = 1 ] || continue
    [ -n "${eligible_since[$b]}" ] || eligible_since[$b]=$(date +%s)
    load=$(awk '{ printf "%d", $1 + 0.999 }' /proc/loadavg)      # declared threads are guesses too: actual load also counts
    if [ $((used_mem + mem)) -le $BUDGET_MEM ] && { [ $((used_cpu + cpu)) -le $BUDGET_CPU ] || [ $((load + cpu)) -le $BUDGET_CPU ]; } && \
       [ $(( $(avail_gb) - mem )) -ge $MIN_AVAIL ] || \
       { [ $nrun -eq 0 ] && [ $(( $(avail_gb) - (mem < SLICE_CAP ? mem : SLICE_CAP) )) -ge $MIN_AVAIL ]; }; then
      [ -n "$first_waiting" ] && [ $(( $(date +%s) - ${eligible_since[$first_waiting]} )) -gt $((RESERVE_MIN * 60)) ] && break
      start_job "$b" "$mem" "$cpu"; unset "eligible_since[$b]"
      used_mem=$((used_mem + mem)); used_cpu=$((used_cpu + cpu)); nrun=$((nrun + 1))
    else
      [ -z "$first_waiting" ] && first_waiting=$b
    fi
  done
  date '+%F %T' > "$ROOT/heartbeat"
  sleep 15
done
