#!/bin/bash
# usage: run.sh <panel> [--continue-run ID --feedback "..."]   (key read from desktop config, never printed)
set -e
cd /home/boris/Projects/ClrAttention
P=$1; shift
export GOOGLE_API_KEY=$(python3 -c "import json,os;d=json.load(open(os.path.expanduser('~/.config/Claude/claude_desktop_config.json')));print(list(d['mcpServers'].values())[0]['env']['GOOGLE_API_KEY'])" )
export GOOGLE_VLM_MODEL=gemini-3.8-flash GOOGLE_IMAGE_MODEL=gemini-nano-banana-2.1
D=docs/whitepapers/figures/pb_runs/panels_20261030
if [ "$1" = "--continue-run" ]; then
  exec nice -n 19 ~/.local/bin/paperbanana generate --output-dir $D/runs --vlm-model gemini-3.8-flash --image-model gemini-nano-banana-2.1 "$@"
fi
exec nice -n 19 ~/.local/bin/paperbanana generate --output-dir $D/runs \
  --vlm-model gemini-3.8-flash --image-model gemini-nano-banana-2.1 \
  -i $D/in/panel$P.txt -c "$(cat $D/in/panel$P.caption)" \
  --image $D/ref/panel$P.png --image $D/in/style_2007_panelA.png -n 1 -ar $(cat $D/in/panel$P.ar)
