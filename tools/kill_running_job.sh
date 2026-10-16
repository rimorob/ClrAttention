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
pids=$(pgrep -f "bash ${JOBS_ROOT:-jobs}/running/")
[ -z "$pids" ] && { echo "no running job"; exit 0; }
# PSOCK workers are separate R processes, not children of the job. Note the
# ports their masters (R processes in this job's tree) listen on, and later
# kill only the workers connected to those ports: other jobs' workers (a
# second job lane) must survive.
tree() { local p=$1 c; echo "$p"; for c in $(pgrep -P "$p"); do tree "$c"; done; }
ports=""
for p in $pids; do for q in $(tree "$p"); do
  ports="$ports $(ss -ltnpH 2>/dev/null | grep "pid=$q," | awk '{n = split($4, a, ":"); print a[n]}')"
done; done
for p in $pids; do echo "killing job tree $p: $(ps -o args= -p "$p")"; kill_tree "$p"; done
sleep 3
for port in $ports; do
  pkill -f "parallel:::.workRSOCK.*PORT=$port( |$)" 2>/dev/null && echo "killed leftover PSOCK workers of port $port"
done
exit 0
