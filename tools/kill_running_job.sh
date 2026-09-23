#!/bin/bash
# Kill the job the runner is currently executing (the "bash jobs/running/*.sh"
# process and everything it started, e.g. Rscript and its PSOCK workers).
# The runner then logs the job as ended (non-zero rc) and moves on to the
# next queued job. Finished bootstrap draws are kept on disk and skipped on
# restart. Use from the control channel: cp tools/kill_running_job.sh jobs/control/
cd "$(dirname "$0")/.." 2>/dev/null || cd "$(pwd)"
kill_tree() {
  local p=$1 c
  for c in $(pgrep -P "$p"); do kill_tree "$c"; done
  kill -TERM "$p" 2>/dev/null
}
pids=$(pgrep -f "bash jobs/running/")
[ -z "$pids" ] && { echo "no running job"; exit 0; }
for p in $pids; do echo "killing job tree $p: $(ps -o args= -p "$p")"; kill_tree "$p"; done
sleep 3
# PSOCK workers are separate R processes; kill any orphans still serving.
pkill -f "parallel:::.workRSOCK|parallel:::.slaveRSOCK" 2>/dev/null && echo "killed leftover PSOCK workers"
exit 0
