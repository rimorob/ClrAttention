#!/bin/bash
# Write tools/progress.sh output to docs/status/progress.txt and an auto-refreshing
# docs/status/progress.html (open in a browser; reloads every 30 s). Run by the
# clr-progress.timer user unit.
cd "$(dirname "$0")/.." || exit 1
mkdir -p docs/status
out=$(tools/progress.sh 2>&1)
printf '%s\n' "$out" > docs/status/progress.txt.new && mv docs/status/progress.txt.new docs/status/progress.txt
esc=$(printf '%s' "$out" | sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g')
cat > docs/status/progress.html.new <<HTML
<!doctype html><html><head><meta charset="utf-8"><meta http-equiv="refresh" content="30">
<title>ClrAttention jobs</title>
<style>body{background:#111;color:#ddd;font:14px/1.45 ui-monospace,Menlo,Consolas,monospace;margin:16px}
pre{white-space:pre}</style></head><body><pre>$esc</pre></body></html>
HTML
mv docs/status/progress.html.new docs/status/progress.html
