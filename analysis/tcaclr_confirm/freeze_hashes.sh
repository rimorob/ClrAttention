#!/bin/bash
# Writes docs/tcaclr_confirm_selections/<arm>.sha256 (SHA-256 of the truth members file, the excluded-pairs file and the expression matrix
# of one arm; names without directories, the matrix named by its basename). The score stage refuses to run on a mismatch.
# Usage: freeze_hashes.sh armA|armB <expression file>
cd "$(dirname "$0")/../.." || exit 1
D=docs/tcaclr_confirm_selections; arm=$1; expr=$2
case $arm in armA) pre=armA_complex_portal ;; armB) pre=armB_cyc2008 ;; *) echo "arm?"; exit 1 ;; esac
out="$D/$arm.sha256"; : > "$out"
for f in "$D/${pre}_members.tsv" "$D/${pre}_excluded_pairs.tsv" "$expr"; do echo "$(sha256sum "$f" | cut -d' ' -f1)  $(basename "$f")" >> "$out"; done
cat "$out"
