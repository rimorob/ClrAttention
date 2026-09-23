#!/bin/bash
# Control channel next to the job runner: executes jobs/control/*.sh
# immediately (not queued behind the running job), e.g. to stop a job.
# Started in the background by a queued job, or by hand:
#     nohup tools/control_runner.sh >/dev/null 2>&1 &
# Logs to jobs/control/control.log. Stop it with: touch jobs/control/quit
# SECURITY: like the job runner, anything placed in jobs/control/*.sh is
# executed as you. Only this repo folder is watched.
cd "$(dirname "$0")/.." || exit 1
mkdir -p jobs/control/done
echo "$(date '+%F %T') control runner started (pid $$)" >> jobs/control/control.log
echo $$ > jobs/control/control_runner.pid
while [ ! -e jobs/control/quit ]; do
  for f in jobs/control/*.sh; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .sh)
    mv "$f" "jobs/control/done/$b.sh" || continue
    echo "$(date '+%F %T') RUN $b" >> jobs/control/control.log
    bash "jobs/control/done/$b.sh" > "jobs/control/done/$b.log" 2>&1
    echo "$(date '+%F %T') END $b rc=$?" >> jobs/control/control.log
  done
  date '+%F %T' > jobs/control/heartbeat
  sleep 5
done
rm -f jobs/control/quit jobs/control/control_runner.pid
echo "$(date '+%F %T') control runner stopped" >> jobs/control/control.log
