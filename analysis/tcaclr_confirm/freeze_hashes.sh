#!/bin/bash
# Writes docs/tcaclr_confirm_selections/<arm>.sha256 (SHA-256 of the truth members file and the excluded-pairs file, names relative to that
# directory: the format of the score stage's hash guard) and <arm>_expression.sha256 (SHA-256 of the expression matrix, by basename).
# Usage: freeze_hashes.sh armA|armB <expression file>
cd "$(dirname "$0")/../.." || exit 1
D=docs/tcaclr_confirm_selections; arm=$1; expr=$2
case $arm in armA) pre=armA_complex_portal ;; armB) pre=armB_cyc2008 ;; *) echo "arm?"; exit 1 ;; esac
: > "$D/$arm.sha256"
for f in "$D/${pre}_members.tsv" "$D/${pre}_excluded_pairs.tsv"; do echo "$(sha256sum "$f" | cut -d' ' -f1)  $(basename "$f")" >> "$D/$arm.sha256"; done
echo "$(sha256sum "$expr" | cut -d' ' -f1)  $(basename "$expr")" > "$D/${arm}_expression.sha256"
cat "$D/$arm.sha256" "$D/${arm}_expression.sha256"
