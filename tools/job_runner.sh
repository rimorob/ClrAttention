#!/bin/bash
# Local job runner: executes shell jobs dropped into jobs/queue/ one at a
# time, in name order, logging to jobs/done/<job>.log and jobs/runner.log.
# Start it once in Terminal from the repo root and leave it running:
#     caffeinate -i tools/job_runner.sh
# (caffeinate keeps the Mac awake while it runs). Stop with Ctrl-C.
# To stop the *current* job without stopping the runner, use the control
# channel (tools/control_runner.sh): cp tools/kill_running_job.sh jobs/control/
# SECURITY: anything placed in jobs/queue/*.sh is executed as you. Only this
# repo folder is watched.
cd "$(dirname "$0")/.." || exit 1
mkdir -p jobs/queue jobs/running jobs/done
echo "$(date '+%F %T') runner started (pid $$) in $PWD" >> jobs/runner.log
while true; do
  for f in jobs/queue/*.sh; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .sh)
    mv "$f" "jobs/running/$b.sh" || continue
    echo "$(date '+%F %T') START $b" >> jobs/runner.log
    bash "jobs/running/$b.sh" > "jobs/done/$b.log" 2>&1
    rc=$?
    mv "jobs/running/$b.sh" "jobs/done/$b.sh"
    echo "$(date '+%F %T') END $b rc=$rc" >> jobs/runner.log
  done
  date '+%F %T' > jobs/heartbeat
  sleep 20
done
