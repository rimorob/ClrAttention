#!/bin/bash
# Local job runner: executes shell jobs dropped into <root>/queue/ one at a
# time, in name order, logging to <root>/done/<job>.log and <root>/runner.log.
# <root> defaults to "jobs" (the Mac); the Linux workstation uses "jobs_pc",
# so the two machines never pick up each other's jobs from the shared
# (iCloud) folder. Start it once from the repo root and leave it running:
#     macOS:  caffeinate -i tools/job_runner.sh
#     Linux:  nohup systemd-inhibit --what=sleep:idle tools/job_runner.sh jobs_pc >/dev/null 2>&1 &
# It also starts the control channel for the same root (tools/control_runner.sh),
# which runs <root>/control/*.sh immediately, e.g. to stop the current job:
#     cp tools/kill_running_job.sh <root>/control/
# The queue is rescanned after every job (name order = priority, including
# jobs added while another was running). A job file is run only once it has
# been unchanged for 30 s, so a file still
# arriving through cloud sync is not run half-written.
# SECURITY: anything placed in <root>/queue/*.sh is executed as you. Only this
# repo folder is watched.
cd "$(dirname "$0")/.." || exit 1
ROOT="${1:-jobs}"
export JOBS_ROOT="$ROOT"
mkdir -p "$ROOT/queue" "$ROOT/running" "$ROOT/done" "$ROOT/control"
mtime() { stat -c %Y "$1" 2>/dev/null || stat -f %m "$1"; }
echo "$(date '+%F %T') runner started (pid $$) on $(hostname) in $PWD, root $ROOT" >> "$ROOT/runner.log"
pgrep -f "tools/control_runner.sh $ROOT\$" >/dev/null || \
  nohup tools/control_runner.sh "$ROOT" >/dev/null 2>&1 &
while true; do
  for f in "$ROOT"/queue/*.sh; do
    [ -e "$f" ] || continue
    [ $(( $(date +%s) - $(mtime "$f") )) -ge 30 ] || continue
    b=$(basename "$f" .sh)
    mv "$f" "$ROOT/running/$b.sh" || continue
    echo "$(date '+%F %T') START $b" >> "$ROOT/runner.log"
    bash "$ROOT/running/$b.sh" > "$ROOT/done/$b.log" 2>&1
    rc=$?
    mv "$ROOT/running/$b.sh" "$ROOT/done/$b.sh"
    echo "$(date '+%F %T') END $b rc=$rc" >> "$ROOT/runner.log"
    break   # rescan the queue after every job, so jobs added meanwhile run in name order
  done
  date '+%F %T' > "$ROOT/heartbeat"
  sleep 20
done
