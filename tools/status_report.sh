#!/bin/bash
# Write a status report of the jobs_pc queue as a Markdown table:
#   docs/status/status_<YYYY-MM-DD>.md  (and docs/status/latest.md)
# Jobs since --since (default: 36 h ago) with start, end, duration and exit
# code from jobs_pc/runner.log; running and queued jobs; memory and load; and
# docs/status/notes.md (interpretation written by the analysis session) if present.
# Usage (repo root):  tools/status_report.sh [--since "2026-10-15 00:00"]
cd "$(dirname "$0")/.." || exit 1
since=$(date -d '36 hours ago' '+%F %T')
[ "${1:-}" = "--since" ] && since=$2
out=docs/status/status_$(date +%F).md
mkdir -p docs/status
{
  echo "# ClrAttention job status, $(date '+%A %F %H:%M %Z')"
  echo
  echo "Runner service: \`$(systemctl --user is-active clr-jobs 2>/dev/null)\`; memory: $(free -g | awk '/^Mem/ {print $3 " GB used, " $7 " GB available of " $2 " GB"}'); load: $(cut -d' ' -f1-3 /proc/loadavg)."
  echo
  echo "## Jobs since $since"
  echo
  echo "| Job | Started | Ended | Duration | Result |"
  echo "|---|---|---|---|---|"
  running=" $(ls jobs_pc/running 2>/dev/null | sed 's/\.sh$//' | paste -sd' ') "
  awk -v since="$since" -v running="$running" '
    { ts = $1 " " $2 }
    $3 == "START" && ts >= since { n++; job[n] = $4; st[n] = ts; open[$4] = n }
    $3 == "END" && ($4 in open) { k = open[$4]; en[k] = ts; split($5, a, "="); rc[k] = a[2]; delete open[$4] }
    END {
      "date +%s" | getline now
      for (k = 1; k <= n; k++) {
        cmd = "date -d \"" st[k] "\" +%s"; cmd | getline s; close(cmd)
        if (k in en) { cmd = "date -d \"" en[k] "\" +%s"; cmd | getline e; close(cmd)
                       res = (rc[k] == "0") ? "done" : "FAILED (rc=" rc[k] ")" }
        else if (index(running, " " job[k] " ") && open[job[k]] == k) { e = now; res = "running" }
        else { e = s; res = "lost (no END logged; killed)" }
        d = e - s; dur = (e > s) ? sprintf("%dh %02dm", int(d / 3600), int((d % 3600) / 60)) : "-"
        printf "| %s | %s | %s | %s | %s |\n", job[k], substr(st[k], 6, 11), (k in en) ? substr(en[k], 6, 11) : "-", dur, res }
    }' jobs_pc/runner.log
  echo
  echo "**Queued (in order):** $(ls jobs_pc/queue 2>/dev/null | sed 's/\.sh$//' | paste -sd, | sed 's/,/, /g')"
  echo
  if [ -s docs/status/notes.md ]; then echo "## Notes"; echo; cat docs/status/notes.md; fi
} > "$out"
cp "$out" docs/status/latest.md
echo "wrote $out"
