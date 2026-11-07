"""DREAM5 Challenge 4 official-scorer helpers (Python side).  Run with `python3 -I`.  No network access, no Synapse.

  convert-densities <scripts_root> <out_dir>
      Reads <scripts_root>/INPUT/probability_densities/Network{1,3,4}_{AUPR,AUROC}.mat (variables X, Y; scipy) and writes plain
      TSV (columns X, Y; 17 significant digits) Network{n}_{AUPR,AUROC}.tsv into <out_dir>, so R can read them without MATLAB.
  port-run <port_py> <gold_dir> <pdf_dir> <pred_dir> <out_tsv> [team ...]
      Runs the faithful Python port (data/dream5/truth_extracted_tools/d5_official_port.py, path given by <port_py>) on
      <pred_dir>/<team>/DREAM5_NetworkInference_<team>_Network{1,3,4}.txt for every team directory (or the teams named), with the
      gold edge files in <gold_dir> (DREAM5_NetworkInference_Edges_Network{n}.tsv) and the .mat densities in <pdf_dir>.  Writes one
      TSV row per team and network: team net P N T L TPL AUROC AUPR P_AUROC P_AUPR.
"""
import sys, os, importlib.util
import numpy as np
import scipy.io as sio

NETS = (1, 3, 4)


def convert_densities(root, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for n in NETS:
        for m in ("AUPR", "AUROC"):
            d = sio.loadmat(os.path.join(root, "INPUT", "probability_densities", f"Network{n}_{m}.mat"))
            X = d["X"].ravel().astype(float); Y = d["Y"].ravel().astype(float)
            assert X.shape == Y.shape and len(X) > 2
            with open(os.path.join(out_dir, f"Network{n}_{m}.tsv"), "w") as f:
                f.write("X\tY\n")
                for x, y in zip(X, Y):
                    f.write("%s\t%s\n" % (repr(float(x)), repr(float(y))))
            print("wrote", n, m, len(X))


def load_port(path):
    spec = importlib.util.spec_from_file_location("d5_official_port", path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod


def port_run(port_py, gold_dir, pdf_dir, pred_dir, out_tsv, teams):
    port = load_port(port_py)
    if not teams:
        teams = sorted(t for t in os.listdir(pred_dir) if os.path.isdir(os.path.join(pred_dir, t)))
    with open(out_tsv, "w") as f:
        f.write("team\tnet\tP\tN\tT\tL\tTPL\tAUROC\tAUPR\tP_AUROC\tP_AUPR\n")
        for t in teams:
            for n in NETS:
                gold = port.load_dream_network(os.path.join(gold_dir, f"DREAM5_NetworkInference_Edges_Network{n}.tsv"))
                pf = os.path.join(pred_dir, t, f"DREAM5_NetworkInference_{t}_Network{n}.txt")
                pred = port.load_dream_network(pf)
                pa = sio.loadmat(os.path.join(pdf_dir, f"Network{n}_AUPR.mat")); po = sio.loadmat(os.path.join(pdf_dir, f"Network{n}_AUROC.mat"))
                r = port.evaluate(gold, pred, pa, po)
                f.write("\t".join([t, str(n)] + [repr(float(r[k])) for k in ("P", "N", "T", "L", "TPL", "AUROC", "AUPR", "P_AUROC", "P_AUPR")]) + "\n")
    print("wrote", out_tsv)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "convert-densities" and len(a) == 3: convert_densities(a[1], a[2])
    elif a and a[0] == "port-run" and len(a) >= 6: port_run(a[1], a[2], a[3], a[4], a[5], a[6:])
    else: sys.exit(__doc__)
