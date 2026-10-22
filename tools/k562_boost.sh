#!/bin/bash
# One-shot: when K562 keff splits 3 and 4 are saved, restart the selection job with 4 workers
# (resumes from cached splits; nothing in progress is lost at that moment).
cd "$(dirname "$0")/.." || exit 1
S=results/keff_test/k562ess/splits
until [ -s $S/split_03.csv ] && [ -s $S/split_04.csv ]; do sleep 60; done
f=jobs_pc/running/0893b_keff_select_k562.sh
[ -e "$f" ] || { echo "$(date '+%F %T') k562 boost: job not running, nothing to do" >> jobs_pc/runner.log; exit 0; }
sed -e 's/^# mem: .*/# mem: 60/' -e 's/^# cpus: .*/# cpus: 8/' -e 's/--workers 2 /--workers 4 /' "$f" > jobs_pc/sched_hold/0893c_keff_select_k562.sh
systemctl --user stop clrjob-0893b_keff_select_k562.scope; sleep 30
touch -d '-1 min' jobs_pc/sched_hold/0893c_keff_select_k562.sh; mv jobs_pc/sched_hold/0893c_keff_select_k562.sh jobs_pc/sched/
echo "$(date '+%F %T') k562 boost: 0893b stopped after splits 3-4; requeued as 0893c with 4 workers" >> jobs_pc/runner.log
