#!/usr/bin/env python
"""Step 3.2(d)/3.5 R1 setup: finish populating the SEPARATE lamin instance ~/scfm_models/lamin_r1 (LAMIN_SETTINGS_DIR=~/scfm_models/lamin_r1/home)
so that scprint 2.2.0's GNInfer can run. No truth file is read.

scdataloader.utils.populate_my_ontology() (v1.0.1) filled cell types, organisms, sex, ethnicity, assays, tissues and developmental stages and then
raised TypeError at its mouse DevelopmentalStage step (`records[-4][0]`); the Disease step and the Gene step were
then run separately (this script); the mouse developmental-stage step was skipped (not needed by GNInfer, which reads only organism). The Gene step differs in one way forced by sqlite ("too many SQL variables" in bt.Gene.from_values on ~60k ids): the Gene
registry is filled in chunks of 400 and only with the Ensembl ids of the scPRINT vocabulary (44,756). The gene table is bionty's public Ensembl
table; ids of the vocabulary that the public table lacks are listed in the output json (they would shift the Collator's positions, so the script
reports whether Collator's gene order equals model.genes).
  LAMIN_SETTINGS_DIR=$HOME/scfm_models/lamin_r1/home python -I analysis/scfm_r1_setup.py --out OUT.json
"""
import argparse, json, os, time
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
a = ap.parse_args()
assert "lamin_r1" in os.environ.get("LAMIN_SETTINGS_DIR", ""), "must run in the separate lamin instance"
import bionty as bt, lamindb as ln
t0 = time.time()
R = {}
if bt.Disease.filter().count() == 0:          # the Disease step of populate_my_ontology (it follows the failing mouse-stage step)
    ln.save(bt.Disease.from_values(bt.Disease.public().df().index, field="ontology_id"), parents=False)
    bt.Disease(name="normal", ontology_id="PATO:0000461").save()
    bt.Disease(name="unknown", ontology_id="unknown").save()
vocab = list(torch.load(os.path.expanduser("~/scfm_models/scprint/medium-v1.5.ckpt"), map_location="cpu", weights_only=False)["hyper_parameters"]["genes"])
R["vocab"] = len(vocab)
have = bt.Gene.filter().count()
R["genes_already_in_registry"] = have
for name, prefix in (("human", "ENSG"), ("mouse", "ENSMUSG")):
    ids = [g for g in vocab if g.startswith(prefix)]
    pub = bt.Gene.public(organism=name).df()
    pub_ids = set(pub["ensembl_gene_id"])
    present = [g for g in ids if g in pub_ids]
    R[f"{name}_vocab"], R[f"{name}_in_public_table"] = len(ids), len(present)
    R[f"{name}_missing_from_public_examples"] = [g for g in ids if g not in pub_ids][:10]
    R[f"{name}_n_missing_from_public"] = len(ids) - len(present)
    for i in range(0, len(present), 400):
        recs = bt.Gene.from_values(present[i:i + 400], field="ensembl_gene_id", organism=name)
        ln.save(recs)
    print(name, len(present), time.time() - t0, flush=True)
R["genes_in_registry_after"] = bt.Gene.filter().count()
R["s"] = time.time() - t0
json.dump(R, open(a.out, "w"), indent=1)
print(json.dumps(R, indent=1))
