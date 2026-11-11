#!/usr/bin/env python3
"""DREAM5 scorer gate inputs: tidy converter for the organisers' published results table and the flat prediction directory.
Run with `python3 -I analysis/dream5_published_table.py <command> ...` (the inputs are untrusted downloads: read only).

  tidy <r-survey4.xls.csv> <out.csv>
      The sheet 'survey4.xls' of DREAM5_NetworkInference_MethodScoringResults.xls (syn2787254), exported to CSV by LibreOffice, is
      a messy sheet: category header rows ('Regression', 'Mutual information', 'Correlation', 'Bayesian networks', 'Other', 'Meta',
      'Community') followed by rows whose column 2 is the within-category number (1..8) or 'RankSum' (Community). 0-based columns:
      2 overall score, 3 AUPR score, 4 AUROC score (the organisers' -log10 p scores averaged over networks 1, 3, 4; overall = mean of
      the two), 5-7 AUPR nets 1, 3, 4, 8-10 AUROC nets 1, 3, 4, 11-12 network 2 (not scored, ignored). There are NO per-network
      p-values in the table. Writes team (e.g. Regression3, MI4, Community), score_overall, score_aupr, score_auroc, aupr1/3/4,
      auroc1/3/4. Asserts the expected row counts (Regression 8, MI 5, Correlation 3, Bayesian 6, Other 8, Meta 5, Community 1 = 36)
      and that every numeric cell parses; fails loudly otherwise.

  flat-pred <Network_predictions dir> <flat_dir>
      Builds <flat_dir>/<team>/DREAM5_NetworkInference_<team>_Network<n>.txt (n = 1..4) as symlinks into the extracted archive, for
      every team directory under 'Challenge participants' and 'Off-the-shelf methods' and for 'Community integration' (files
      DREAM5_NetworkInference_Community_Network<n>.txt -> flat_dir/Community/). The extracted files are never modified. Asserts that
      each team has all four files and that team names are unique across the three folders.
"""
import csv, os, re, sys

CATEGORY = {"Regression": "Regression", "Mutual information": "MI", "Correlation": "Correlation", "Bayesian networks": "Bayesian",
            "Other": "Other", "Meta": "Meta", "Community": "Community"}
EXPECTED = {"Regression": 8, "MI": 5, "Correlation": 3, "Bayesian": 6, "Other": 8, "Meta": 5, "Community": 1}
OUT_COLS = ["team", "score_overall", "score_aupr", "score_auroc", "aupr1", "aupr3", "aupr4", "auroc1", "auroc3", "auroc4"]
SRC_COLS = {"score_overall": 2, "score_aupr": 3, "score_auroc": 4, "aupr1": 5, "aupr3": 6, "aupr4": 7, "auroc1": 8, "auroc3": 9, "auroc4": 10}


def tidy_rows(rows):
    """rows: csv rows (lists of str) of the exported sheet. Returns a list of dicts (OUT_COLS)."""
    out, cat = [], None
    for k, r in enumerate(rows):
        c1 = r[1].strip() if len(r) > 1 else ""
        if not c1:
            continue
        if c1 in CATEGORY:                       # category header row (no numbers in it)
            assert not any(x.strip() for x in r[2:13]), "row %d: category header %r carries numbers" % (k, c1)
            cat = CATEGORY[c1]; continue
        if c1 in ("Method",):                    # sheet title rows above the data
            continue
        if cat is None:
            continue
        if cat == "Community":
            assert c1 == "RankSum", "row %d: Community row is %r, expected 'RankSum'" % (k, c1)
            team = "Community"
        else:
            assert re.fullmatch(r"[1-9][0-9]*", c1), "row %d: category %s has row label %r, expected an integer" % (k, cat, c1)
            team = "%s%s" % (cat, c1)
        d = {"team": team}
        for name, j in SRC_COLS.items():
            assert j < len(r) and r[j].strip() != "", "row %d (%s): column %d (%s) is empty" % (k, team, j, name)
            try: d[name] = float(r[j])
            except ValueError: raise AssertionError("row %d (%s): column %d (%s) is not numeric: %r" % (k, team, j, name, r[j]))
        out.append(d)
    cnt = {}
    for d in out: cnt[re.sub(r"[0-9]+$", "", d["team"])] = cnt.get(re.sub(r"[0-9]+$", "", d["team"]), 0) + 1
    assert cnt == EXPECTED, "unexpected row counts %r, expected %r" % (cnt, EXPECTED)
    assert len(out) == 36 and len({d["team"] for d in out}) == 36, "expected 36 distinct teams, got %d" % len(out)
    for d in out:   # the table's own definition (overall = mean of the two scores) is checked against its rounding (3 decimals; 4 in some cells)
        assert abs(d["score_overall"] - (d["score_aupr"] + d["score_auroc"]) / 2) <= 2e-3, "%s: overall is not the mean of the two scores" % d["team"]
    return out


def tidy(src, dst):
    with open(src, newline="") as f: rows = list(csv.reader(f))
    out = tidy_rows(rows)
    with open(dst, "w", newline="") as f:
        w = csv.writer(f); w.writerow(OUT_COLS)
        for d in out: w.writerow([d["team"]] + [repr(d[c]) for c in OUT_COLS[1:]])
    print("wrote %s: %d teams" % (dst, len(out)))


def flat_pred(root, flat):
    root = os.path.abspath(root); flat = os.path.abspath(flat)
    jobs = {}   # team -> {n: source path}
    for sub in ("Challenge participants", "Off-the-shelf methods"):
        base = os.path.join(root, sub)
        for t in sorted(os.listdir(base)):
            if not os.path.isdir(os.path.join(base, t)): continue
            assert t not in jobs, "team %s appears twice" % t
            jobs[t] = {n: os.path.join(base, t, "DREAM5_NetworkInference_%s_Network%d.txt" % (t, n)) for n in (1, 2, 3, 4)}
    jobs["Community"] = {n: os.path.join(root, "Community integration", "DREAM5_NetworkInference_Community_Network%d.txt" % n) for n in (1, 2, 3, 4)}
    for t, fs in jobs.items():
        for n, p in fs.items(): assert os.path.isfile(p), "missing prediction file %s" % p
    for t, fs in jobs.items():
        d = os.path.join(flat, t); os.makedirs(d, exist_ok=True)
        for n, p in fs.items():
            link = os.path.join(d, "DREAM5_NetworkInference_%s_Network%d.txt" % (t, n))
            if os.path.islink(link): os.remove(link)
            assert not os.path.exists(link), "refusing to replace the non-link %s" % link
            os.symlink(p, link)
    print("wrote %s: %d teams x 4 symlinks" % (flat, len(jobs)))


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == "tidy": tidy(a[1], a[2])
    elif len(a) == 3 and a[0] == "flat-pred": flat_pred(a[1], a[2])
    else: sys.exit(__doc__)
