#!/usr/bin/env python3
"""DREAM5 network-inference challenge: staged, logged download (design: docs/plans/dream5_design.md).

Run with the dedicated environment:   ~/venvs/dream5/bin/python analysis/dream5_download.py <stage>

Stages (a stage downloads only its own whitelist of Synapse IDs; nothing else can be fetched):
  listing   anonymous check of the Synapse tree against the manifest below (public metadata,
            no login, no download).
  public    organiser-hosted, login-free copy of the EXPRESSION matrices only (GeneNetWeaver site;
            transposed and replicate-averaged formats; no TF lists, no gene-ID maps, no truth).
  data      Synapse, training data of networks 1, 3, 4: expression, chip features, gene IDs, TF
            lists (+ the small 'original data' name/experiment tables for E. coli and yeast).
            NO gold standard. Needs the PI's own Synapse login.
  scoring   Synapse, INSTRUCTIONS.txt and 'Evaluation scripts.zip' (scoring code only). The zip
            is listed and refused if any member name looks like a gold standard.
  truth     Synapse, the five gold-standard files and the published challenge results/predictions.
            REFUSED unless the design is fixed and pushed and the data log is committed (see gate()).

Credentials. This script never asks for, reads, prints or stores a token or password. For the
stages that need Synapse it calls synapseclient's own login(), which uses the PI's existing
SYNAPSE_AUTH_TOKEN environment variable or ~/.synapseConfig. Nothing credential-like is written
by this script. The download cache is redirected into the output directory (data/, git-ignored).

Every file is verified against the expected name, then its size, SHA-256 and the Synapse-reported
MD5 are appended to docs/dream5_selections/download_log.tsv (commit that log; it is the record).
"""
import argparse, csv, datetime, hashlib, json, os, subprocess, sys, urllib.request, zipfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN = "docs/plans/dream5_design.md"
LOG = os.path.join(REPO, "docs/dream5_selections/download_log.tsv")
SYN_REPO = "https://repo-prod.prod.sagebase.org/repo/v1/"
GNW_URL = "https://gnw.sourceforge.net/resources/DREAM5_NetworkInferenceChallenge_AlternativeDataFormats.zip"

# (synapse id, expected file name, subdirectory under <out>/synapse)
NETS = {"net1": ("Network 1 - in silico", 2787224), "net3": ("Network 3 - E. coli", 2787232),
        "net4": ("Network 4 - S. cerevisiae", 2787236)}
DATA = [  # training data (anonymised IDs), networks 1, 3, 4 only; S. aureus (net2) has no usable gold standard
    ("syn2787224", "net1_gene_ids.tsv", "net1"), ("syn2787225", "net1_chip_features.tsv", "net1"),
    ("syn2787226", "net1_expression_data.tsv", "net1"), ("syn2787227", "net1_transcription_factors.tsv", "net1"),
    ("syn2787232", "net3_gene_ids.tsv", "net3"), ("syn2787233", "net3_chip_features.tsv", "net3"),
    ("syn2787234", "net3_expression_data.tsv", "net3"), ("syn2787235", "net3_transcription_factors.tsv", "net3"),
    ("syn2787236", "net4_gene_ids.tsv", "net4"), ("syn2787237", "net4_chip_features.tsv", "net4"),
    ("syn2787238", "net4_expression_data.tsv", "net4"), ("syn2787239", "net4_transcription_factors.tsv", "net4"),
    # 'original data': real gene / TF names and experiment tables only (not the expression matrices)
    ("syn3130841", "README.txt", "original"),
    ("syn3130846", "ecoli_experiments.tsv", "original/Ecoli"), ("syn3130847", "ecoli_gene_names.tsv", "original/Ecoli"),
    ("syn3130848", "ecoli_tf_names.tsv", "original/Ecoli"),
    ("syn3130854", "yeast_experiments.tsv", "original/Scerevisiae"), ("syn3130855", "yeast_gene_names.tsv", "original/Scerevisiae"),
    ("syn3130856", "yeast_tf_names.tsv", "original/Scerevisiae"),
]
SCORING = [("syn2787218", "INSTRUCTIONS.txt", "scoring"), ("syn2787219", "Evaluation scripts.zip", "scoring")]
TRUTH = [
    ("syn2787240", "DREAM5_NetworkInference_GoldStandard_Network1 - in silico.tsv", "truth"),
    ("syn2787241", "DREAM5_NetworkInference_GoldStandardSigned_Network1 - in silico.tsv", "truth"),
    ("syn2787243", "DREAM5_NetworkInference_GoldStandard_Network3 - E. coli.tsv", "truth"),
    ("syn2787244", "DREAM5_NetworkInference_GoldStandard_Network4 - S. cerevisiae.tsv", "truth"),
    ("syn2787254", "DREAM5_NetworkInference_MethodScoringResults.xls", "published"),
    ("syn2787251", "Network_predictions.zip", "published"),
]  # network 2 (S. aureus) gold standard (syn2787242) is deliberately not fetched
MANIFEST = {"data": DATA, "scoring": SCORING, "truth": TRUTH}


def sh(*a):
    return subprocess.run(a, cwd=REPO, capture_output=True, text=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def log_row(stage, sid, name, path, version="", md5=""):
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    new = not os.path.exists(LOG)
    with open(LOG, "a", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        if new:
            w.writerow(["utc", "stage", "source_id", "name", "version", "bytes", "sha256", "synapse_md5", "local_path"])
        w.writerow([datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), stage, sid, name,
                    version, os.path.getsize(path), sha256(path), md5, os.path.relpath(path, REPO)])


def check_out(out):
    r = sh("git", "check-ignore", "-q", out)
    if r.returncode != 0:
        sys.exit(f"refusing: {out} is not git-ignored (downloads must stay under data/)")


def gate():
    """Truth stage only: design FIXED and pushed, data log committed and pushed (house style)."""
    problems = []
    if sh("git", "ls-files", "--error-unmatch", DESIGN).returncode != 0:
        problems.append(f"{DESIGN} is not tracked")
    if sh("git", "diff", "--quiet", "HEAD", "--", DESIGN).returncode != 0:
        problems.append(f"{DESIGN} differs from HEAD")
    txt = open(os.path.join(REPO, DESIGN)).read() if os.path.exists(os.path.join(REPO, DESIGN)) else ""
    if not any(l.startswith("Status: FIXED ") for l in txt.splitlines()):
        problems.append("no 'Status: FIXED <date>' line at the start of a line in the design (it is still a draft)")
    if sh("git", "rev-parse", "--is-shallow-repository").stdout.strip() == "true":
        problems.append("shallow clone")
    up = sh("git", "rev-parse", "--abbrev-ref", "@{u}")
    if up.returncode != 0:
        problems.append("no upstream branch")
    for p in (DESIGN, "docs/dream5_selections/download_log.tsv"):
        c = sh("git", "log", "-1", "--format=%H", "--", p).stdout.strip()
        if not c:
            problems.append(f"{p} has no commit")
        elif up.returncode == 0 and sh("git", "merge-base", "--is-ancestor", c, "@{u}").returncode != 0:
            problems.append(f"last commit of {p} ({c[:7]}) is not on the upstream branch (push first)")
    # the log (data stage, no truth) must be committed and pushed before any gold standard is fetched
    if sh("git", "diff", "--quiet", "HEAD", "--", "docs/dream5_selections/download_log.tsv").returncode != 0:
        problems.append("download_log.tsv has uncommitted changes")
    if problems:
        sys.exit("truth stage REFUSED:\n  - " + "\n  - ".join(problems))
    c = sh("git", "log", "-1", "--format=%H", "--", DESIGN).stdout.strip()
    print(f"gate passed; design commit {c}")
    return c


def listing():
    def post(p, d):
        r = urllib.request.Request(SYN_REPO + p, json.dumps(d).encode(), {"Content-Type": "application/json"})
        return json.load(urllib.request.urlopen(r, timeout=60))
    found = {}

    def walk(i):
        for c in post("entity/children", {"parentId": i, "includeTypes": ["folder", "file"]})["page"]:
            found[c["id"]] = c["name"]
            if c["type"].endswith("Folder"):
                walk(c["id"])
    walk("syn2787209"); walk("syn3130840")
    bad = 0
    for st, items in MANIFEST.items():
        for sid, name, _ in items:
            ok = found.get(sid) == name
            bad += not ok
            print(("ok   " if ok else "MISMATCH ") + f"{st:8s} {sid} {name}" + ("" if ok else f"  (Synapse: {found.get(sid)})"))
    sys.exit(1 if bad else 0)


def public(out):
    d = os.path.join(out, "public_gnw"); os.makedirs(d, exist_ok=True)
    z = os.path.join(d, "DREAM5_NetworkInferenceChallenge_AlternativeDataFormats.zip")
    if not os.path.exists(z):
        urllib.request.urlretrieve(GNW_URL, z)
    with zipfile.ZipFile(z) as zf:
        names = zf.namelist()
        if any("gold" in n.lower() for n in names):
            sys.exit("unexpected gold-standard-like member in the public zip; not extracted")
        zf.extractall(d)
    log_row("public", GNW_URL, os.path.basename(z), z)
    print("public expression archive:", z, os.path.getsize(z), "bytes", sha256(z))


def synapse_stage(stage, out, dry):
    items = MANIFEST[stage]
    if dry:
        for sid, name, sub in items:
            print(f"would download {sid}  {name}  ->  {out}/synapse/{sub}/")
        return
    import synapseclient
    root = os.path.join(out, "synapse")
    os.makedirs(root, exist_ok=True)
    syn = synapseclient.Synapse(cache_root_dir=os.path.join(out, "_synapse_cache"), skip_checks=True)
    try:
        syn.login(silent=True)      # PI's own SYNAPSE_AUTH_TOKEN or ~/.synapseConfig; nothing is passed or stored here
    except Exception as e:
        sys.exit(f"Synapse login failed ({type(e).__name__}). Set SYNAPSE_AUTH_TOKEN in your own shell or create "
                 f"~/.synapseConfig yourself (see Synapse docs); this script never handles credentials.")
    for sid, name, sub in items:
        dest = os.path.join(root, sub); os.makedirs(dest, exist_ok=True)
        try:
            ent = syn.get(sid, downloadLocation=dest, ifcollision="overwrite.local")
        except Exception as e:
            sys.exit(f"download of {sid} ({name}) failed: {type(e).__name__}: {str(e)[:300]}\n"
                     f"If this is a 403, open https://www.synapse.org/Synapse:syn2787209 once in your browser and accept the "
                     f"challenge terms, then rerun.")
        path = os.path.join(dest, ent.name)
        if ent.name != name:
            sys.exit(f"{sid}: expected '{name}', Synapse gave '{ent.name}'; stopping")
        fh = getattr(ent, "_file_handle", {}) or {}
        log_row(stage, sid, name, path, version=str(getattr(ent, "versionNumber", "")), md5=fh.get("contentMd5", ""))
        print(f"{sid} {name}  {os.path.getsize(path)} bytes  sha256 {sha256(path)[:16]}...")
        if stage == "scoring" and name.endswith(".zip"):
            with zipfile.ZipFile(path) as zf:
                bad = [n for n in zf.namelist() if "gold" in n.lower() or "standard" in n.lower()]
                print("  zip members:", len(zf.namelist()), "; gold-like members:", bad or "none")
                if bad:
                    sys.exit("scoring zip contains gold-standard-like members; not extracting. Tell the lead.")
    with open(os.path.join(root, "SHA256SUMS"), "a") as f:
        for sid, name, sub in items:
            p = os.path.join(root, sub, name)
            f.write(f"{sha256(p)}  {sub}/{name}\n")
    sh("rm", "-rf", os.path.join(out, "_synapse_cache"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["listing", "public", "data", "scoring", "truth"])
    ap.add_argument("--out", default=os.path.join(REPO, "data/dream5"))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.stage == "listing":
        listing()
    check_out(a.out)
    if a.stage == "public":
        public(a.out)
    else:
        if a.stage == "truth" and not a.dry_run:
            gate()
        synapse_stage(a.stage, a.out, a.dry_run)
