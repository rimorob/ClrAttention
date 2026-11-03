#!/usr/bin/env python
"""Truth-free QC table of the fine-tuning arm (E17): training summary (loss, steps, weight change) from the per-case training json, and the distance
between the fine-tuned and the zero-shot R0 (relative L1 and Spearman on the off-diagonal), split-half stability and organism-head accuracy of both.
Writes results/scfm_ft/qc_<arm>.csv. NO TRUTH FILE IS READ.
  python -I analysis/scfm_ft_qc.py scprint|scgpt
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy.stats import rankdata
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FT = os.path.expanduser("~/scfm_models/step5_ft"); S3 = os.path.expanduser("~/scfm_models/step3")
arm = sys.argv[1]


def off(a, b):
    m = ~np.eye(a.shape[0], dtype=bool)
    x, y = a[m].astype(np.float64), b[m].astype(np.float64)
    return float(np.abs(x - y).sum() / (0.5 * (x.sum() + y.sum()))), float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


rows = []
if arm == "scprint":
    cases = [f"{d}_{n}" for d in ("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L") for n in (500, 1000)]
    for c in cases:
        t = json.load(open(f"{FT}/scprint_ckpt/{c}.json")); e = json.load(open(f"{FT}/scprint_attn/{c}.json")); z = json.load(open(f"{S3}/scprint_medium_A/{c}.json"))
        a = np.load(f"{FT}/scprint_attn/{c}.nets.npz", allow_pickle=True); b = np.load(f"{S3}/scprint_medium_A/{c}.nets.npz", allow_pickle=True)
        assert list(a["genes"]) == list(b["genes"])
        rel, sp = off(a["R0"], b["R0"])
        relp, spp = off(a["R0p"], b["R0p"]); rel2, sp2 = off(a["R2"], b["R2"])
        rows.append(dict(case=c, n_cells=t["n_cells"], G=t["n_genes"], optimizer_steps=t["optimizer_steps"], micro_bs=t["micro_bs"], loss_epoch1=t["epoch_mean_loss"][0],
                         loss_epoch10=t["epoch_mean_loss"][-1], probe_before=t["probe_loss_before"], probe_after=t["probe_loss_after"],
                         rel_weight_change=t["rel_weight_change_all"], train_s=t["train_s"], state_sha256=t["state_sha256"],
                         ft_vs_zs_R0_relL1=rel, ft_vs_zs_R0_spearman=sp, ft_vs_zs_R0p_relL1=relp, ft_vs_zs_R2_relL1=rel2,
                         splithalf_R0_ft=e["split_half"]["R0"]["relL1"], splithalf_R0_zs=z["split_half"]["R0"]["relL1"],
                         orgacc_ft=e["organism_accuracy"], orgacc_zs=z["organism_accuracy"], r0p_heads_ft=e["r0p_n_heads"], r0p_heads_zs=z["r0p_n_heads"],
                         row_sum_min=e["row_sum_min"], row_sum_max=e["row_sum_max"], nan=e["attn_nan"] + e["attn_inf"]))
else:
    cases = ["hESC_500", "hESC_1000", "hHep_500", "hHep_1000"]
    for c in cases:
        t = json.load(open(f"{FT}/scgpt_ckpt/{c}.json")); e = json.load(open(f"{FT}/scgpt_attn/{c}.json")); z = json.load(open(f"{S3}/scgpt_A/{c}.json"))
        a = np.load(f"{FT}/scgpt_attn/{c}.nets.npz", allow_pickle=True); b = np.load(f"{S3}/scgpt_A/{c}.nets.npz", allow_pickle=True)
        assert list(a["genes"]) == list(b["genes"])
        rel, sp = off(a["attn_last_mean"], b["attn_last_mean"])
        rt, st = off(a["tok_cos"], b["tok_cos"]) if False else (float(np.abs(a["tok_cos"] - b["tok_cos"]).mean()), float(np.corrcoef(a["tok_cos"].ravel(), b["tok_cos"].ravel())[0, 1]))
        rh, sh = (float(np.abs(a["hid_cos"] - b["hid_cos"]).mean()), float(np.corrcoef(a["hid_cos"].ravel(), b["hid_cos"].ravel())[0, 1]))
        rows.append(dict(case=c, n_cells=t["n_cells"], G=t["n_genes"], optimizer_steps=t["optimizer_steps"], micro_bs=t["micro_bs"], loss_epoch1=t["epoch_mean_loss"][0],
                         loss_epoch10=t["epoch_mean_loss"][-1], probe_before=t["probe_loss_before"], probe_after=t["probe_loss_after"],
                         rel_weight_change=t["rel_weight_change_all"], train_s=t["train_s"], state_sha256=t["state_sha256"],
                         ft_vs_zs_attn_relL1=rel, ft_vs_zs_attn_spearman=sp, ft_vs_zs_tokcos_meanabs=rt, ft_vs_zs_tokcos_corr=st,
                         ft_vs_zs_hidcos_meanabs=rh, ft_vs_zs_hidcos_corr=sh, cls_mass_last_max_ft=max(e["cls_token_mass"][-1]), cls_mass_last_max_zs=max(z["cls_token_mass"][-1]),
                         row_sum_min=e["row_sum_min"], row_sum_max=e["row_sum_max"], nan=e["attn_nan"] + e["attn_inf"] + e["hid_cos_nan"]))
df = pd.DataFrame(rows)
os.makedirs(f"{ROOT}/results/scfm_ft", exist_ok=True)
df.to_csv(f"{ROOT}/results/scfm_ft/qc_{arm}.csv", index=False)
pd.set_option("display.width", 250, "display.max_columns", 40)
print(df.drop(columns=["state_sha256"]).round(4).to_string(index=False))
