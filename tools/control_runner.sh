#!/bin/bash
# Control channel next to the job runner: executes <root>/control/*.sh
# immediately (not queued behind the running job), e.g. to stop a job.
# Started by tools/job_runner.sh, or by hand:
#     nohup tools/control_runner.sh [root] >/dev/null 2>&1 &
# <root> defaults to "jobs". Logs to <root>/control/control.log.
# Stop it with: touch <root>/control/quit
# SECURITY: like the job runner, anything placed in <root>/control/*.sh is
# executed as you. Only this repo folder is watched.
cd "$(dirname "$0")/.." || exit 1
ROOT="${1:-jobs}"
export JOBS_ROOT="$ROOT"
C="$ROOT/control"
mkdir -p "$C/done"
echo "$(date '+%F %T') control runner started (pid $$) on $(hostname)" >> "$C/control.log"
echo $$ > "$C/control_runner.pid"
while [ ! -e "$C/quit" ]; do
  for f in "$C"/*.sh; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .sh)
    mv "$f" "$C/done/$b.sh" || continue
    echo "$(date '+%F %T') RUN $b" >> "$C/control.log"
    bash "$C/done/$b.sh" > "$C/done/$b.log" 2>&1
    echo "$(date '+%F %T') END $b rc=$?" >> "$C/control.log"
  done
  date '+%F %T' > "$C/heartbeat"
  sleep 5
done
rm -f "$C/quit" "$C/control_runner.pid"
echo "$(date '+%F %T') control runner stopped" >> "$C/control.log"
