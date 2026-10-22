#!/bin/bash
# Download BEELINE experimental scRNA-seq inputs and ground-truth networks
# (Zenodo record 3701939) into data/beeline.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/beeline && cd data/beeline
[ -f BEELINE-data.zip ] || curl -L --fail -o BEELINE-data.zip https://zenodo.org/api/records/3701939/files/BEELINE-data.zip/content
[ -f BEELINE-Networks.zip ] || curl -L --fail -o BEELINE-Networks.zip https://zenodo.org/api/records/3701939/files/BEELINE-Networks.zip/content
unzip -qo BEELINE-data.zip 'BEELINE-data/inputs/scRNA-Seq/*'
unzip -qo BEELINE-Networks.zip
ls BEELINE-data/inputs/scRNA-Seq Networks
