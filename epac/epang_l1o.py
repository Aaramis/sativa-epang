#!/usr/bin/env python3
"""Leave-one-out via EPA-ng, replacing RAxML `-f O` in SATIVA.

Produces a placement list in the SAME format as EpaJsonParser.get_placement(),
with edge numbers = the refjson B= numbering (the one classify_seq expects via
bid_taxonomy_map). The classification/decision stays 100% SATIVA.

Approach: k-fold. For each fold, its leaves are removed from the reference tree +
alignment and the fold's sequences are placed with EPA-ng; the EPA-ng edges are
remapped to B= by leaf bipartition. All placements are kept
(--filter-acc-lwr 0.99999 --filter-max 100000) to recover the full LWR mass.
"""
import os, sys, re, json, glob, subprocess, shutil
sys.setrecursionlimit(200000)
from ete3 import Tree

# epa-ng resolved from PATH (provided by the sativa.yaml conda env); overridable
# via SATIVA_EPANG_BIN.
EPANG = os.environ.get("SATIVA_EPANG_BIN") or shutil.which("epa-ng") or "epa-ng"

def _bip_map(nhx_tree_str, tag):
    if tag == "EDGE":
        nhx_tree_str = re.sub(r"\{(\d+)\}", r"[&&NHX:EDGE=\1]", nhx_tree_str)
    t = Tree(nhx_tree_str, format=1)
    allL = frozenset(t.get_leaf_names())
    m = {}
    for n in t.traverse():
        e = getattr(n, tag, None)
        if e is None:
            continue
        desc = frozenset(n.get_leaf_names())
        side = desc if len(desc) <= len(allL) - len(desc) else allL - desc
        m[frozenset(side)] = str(e)
    return m, allL

def _read_fasta(path):
    seqs, cur = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                cur = line[1:].split()[0]; seqs[cur] = []
            else:
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

def _write_fasta(d, path):
    with open(path, "w") as f:
        for k, v in d.items():
            f.write(">%s\n%s\n" % (k, v))

# --- matching RAxML's placement settings -------------------------------------------
# SATIVA's RAxML leave-one-out and EPA-ng do not place under the same rules out of the box.
# On a reference tree of <=1000 taxa SATIVA runs `raxmlHPC -f O` with NO preplacement
# heuristic, RAxML-style branch-length optimisation, and keeps placements up to an
# accumulated LWR of 0.999. EPA-ng defaults to a two-phase heuristic (--dyn-heur 0.99999),
# a faster "sliding" branch-length optimisation, and here we kept an accumulated LWR of
# 0.99999. Each of these can be lined up with RAxML through an environment variable, so the
# defaults stay fast and the strict settings are available for concordance work:
#
#   SATIVA_EPANG_HEUR=off      -> --no-heur      (evaluate every branch, as RAxML does)
#   SATIVA_EPANG_BLO=raxml     -> --raxml-blo    (RAxML-style branch-length optimisation)
#   SATIVA_EPANG_ACC_LWR=0.999 -> --filter-acc-lwr 0.999 (RAxML's threshold)
#   SATIVA_EPANG_FOLDS=<N>     -> one sequence per fold = the strict leave-one-out
#
# Setting all four makes EPA-ng answer the same question as RAxML `-f O`; what remains is
# the likelihood implementation itself.
def epang_placement_flags():
    flags = ["--filter-acc-lwr", os.environ.get("SATIVA_EPANG_ACC_LWR", "0.99999"),
             "--filter-max", "100000"]
    if os.environ.get("SATIVA_EPANG_HEUR", "on").lower() in ("off", "no", "0", "false"):
        flags.append("--no-heur")
    if os.environ.get("SATIVA_EPANG_BLO", "sliding").lower() == "raxml":
        flags.append("--raxml-blo")
    return flags


def run_epang_l1o(refjson_tree_str, refaln_path, reftree_path, raxml_outdir,
                  workdir, folds=5, threads=1, log=None):
    def _log(m):
        if log: log.info("[epang-l1o] " + m)
        else: sys.stderr.write("[epang-l1o] " + m + "\n")

    # EPA-ng model: RAxML_info.mfresolv if present, otherwise GTR+G (EPA-ng re-evaluates)
    info = glob.glob(os.path.join(raxml_outdir, "RAxML_info.mfresolv*"))
    model = info[0] if info else "GTR+G"
    _log("EPA-ng model: %s" % model)

    bidmap, allL = _bip_map(refjson_tree_str, "B")        # bipartition -> B-id
    aln = _read_fasta(refaln_path)
    # align the alignment keys onto the leaf names
    leaves = list(allL)
    aln_by_leaf = {}
    for lf in leaves:
        for cand in (lf, lf[2:] if lf.startswith("r_") else "r_"+lf):
            if cand in aln:
                aln_by_leaf[lf] = aln[cand]; break
    missing = [l for l in leaves if l not in aln_by_leaf]
    if missing:
        raise RuntimeError("epang-l1o: %d leaves without a sequence (e.g. %s)" % (len(missing), missing[:3]))

    full_tree = Tree(reftree_path, format=1)
    ordered = sorted(leaves)
    K = min(folds, len(ordered))
    folds_list = [ordered[i::K] for i in range(K)]

    os.makedirs(workdir, exist_ok=True)
    placements = []
    for fi, fold in enumerate(folds_list):
        fold_set = set(fold)
        ref_leaves = [l for l in leaves if l not in fold_set]
        if len(ref_leaves) < 4 or not fold:
            continue
        wd = os.path.join(workdir, "fold_%03d" % fi); os.makedirs(wd, exist_ok=True)
        tpr = full_tree.copy(method="newick")
        tpr.prune(ref_leaves, preserve_branch_length=True)
        tpr.write(outfile=os.path.join(wd, "ref.nwk"), format=5)
        _write_fasta({l: aln_by_leaf[l] for l in ref_leaves}, os.path.join(wd, "ref.fasta"))
        _write_fasta({l: aln_by_leaf[l] for l in fold},       os.path.join(wd, "query.fasta"))
        cmd = [EPANG, "-t", os.path.join(wd, "ref.nwk"), "-s", os.path.join(wd, "ref.fasta"),
               "-q", os.path.join(wd, "query.fasta"), "-m", model, "--outdir", wd, "--redo",
               "-T", str(threads)] + epang_placement_flags()
        r = subprocess.run(cmd, capture_output=True, text=True)
        if os.environ.get("SATIVA_EPANG_DEBUG"):
            # What EPA-ng makes of the model file it was handed. Above 500 taxa SATIVA
            # builds the reference tree under GTRCAT, and a CAT RAxML_info carries
            # "alpha: 1.000000" -- a placeholder, since CAT fits no gamma shape.
            for line in (r.stdout or "").splitlines():
                if any(k in line.lower() for k in ("model", "alpha", "rate")):
                    _log("epa-ng says: " + line.strip())
        jpf = os.path.join(wd, "epa_result.jplace")
        if r.returncode != 0 or not os.path.isfile(jpf):
            _log("EPA-ng FAIL fold %d: %s" % (fi, r.stderr[-300:])); continue
        d = json.load(open(jpf))
        ie = d["fields"].index("edge_num")
        e2side, foldL = _bip_map(d["tree"], "EDGE")
        # table bipartition-restreinte -> B
        restricted = {}
        for side, b in ((s, bidmap[s]) for s in bidmap):
            s2 = side & foldL
            key = s2 if len(s2) <= len(foldL)-len(s2) else (foldL - s2)
            restricted.setdefault(frozenset(key), b)
        epa2b = {}
        for side, e in e2side.items():
            key = side if len(side) <= len(foldL)-len(side) else (foldL - side)
            b = restricted.get(frozenset(key))
            if b is not None:
                epa2b[int(e)] = b
        # Edges that fail to map back to SATIVA's B= numbering are dropped, and with them
        # their likelihood weight -- which shifts every confidence classify_seq computes.
        # SATIVA_EPANG_DEBUG reports how much mass that is.
        ilwr = d["fields"].index("like_weight_ratio") if "like_weight_ratio" in d["fields"] else None
        kept_mass = dropped_mass = 0.0
        kept_edges = dropped_edges = 0
        for pl in d["placements"]:
            name = (pl.get("n") or pl.get("nm"))[0]
            if isinstance(name, list): name = name[0]
            newp = []
            for row in pl["p"]:
                b = epa2b.get(int(row[ie]))
                if b is None:
                    dropped_edges += 1
                    if ilwr is not None: dropped_mass += float(row[ilwr])
                    continue
                kept_edges += 1
                if ilwr is not None: kept_mass += float(row[ilwr])
                rr = list(row); rr[ie] = int(b)
                newp.append(rr)
            if newp:
                placements.append({"p": newp, "n": [name]})
        if os.environ.get("SATIVA_EPANG_DEBUG"):
            total = kept_mass + dropped_mass
            _log("fold %d: %d edges kept, %d dropped; LWR mass dropped %.4f%%"
                 % (fi, kept_edges, dropped_edges,
                    100.0 * dropped_mass / total if total else 0.0))
    _log("placements produits: %d (K=%d folds)" % (len(placements), K))
    return placements


def run_epang_final(reftree_path, refaln_path, raxml_outdir, workdir, threads=1, log=None):
    """Pass 2 (confirmation) via EPA-ng, replacing the RAxML `-f v` of run_epa_once.

    Places the 'suspect' sequences (those pruned from the tree = reference alignment MINUS
    the pruned tree's leaves) onto the mislabel-free reference tree, and writes a jplace
    directly consumable by EpaJsonParser (tree {N} + self-consistent placements: the
    bid_tax_map is rebuilt from that tree by SATIVA).
    Returns the jplace path, or None if there is no suspect sequence.
    """
    def _log(m):
        (log.info if log else (lambda x: sys.stderr.write(x + "\n")))("[epang-final] " + m)

    os.makedirs(workdir, exist_ok=True)
    info = glob.glob(os.path.join(raxml_outdir, "RAxML_info.mfresolv*"))
    model = info[0] if info else "GTR+G"
    # SATIVA's own pass 2 deliberately does NOT reuse the reference model here
    # ("don't load the model, since it's invalid for the pruned tree", run_epa_once), while
    # we hand EPA-ng the full-tree RAxML_info. SATIVA_EPANG_FINAL_MODEL overrides it so
    # that choice can be measured rather than assumed.
    model = os.environ.get("SATIVA_EPANG_FINAL_MODEL", model)

    tree = Tree(reftree_path, format=1)
    ref_leaves = set(tree.get_leaf_names())
    aln = _read_fasta(refaln_path)
    ref_seqs = {n: s for n, s in aln.items() if n in ref_leaves}
    query_seqs = {n: s for n, s in aln.items() if n not in ref_leaves}
    if not query_seqs:
        _log("no suspect sequence to re-place")
        return None
    _write_fasta(ref_seqs, os.path.join(workdir, "ref.fasta"))
    _write_fasta(query_seqs, os.path.join(workdir, "query.fasta"))

    cmd = [EPANG, "-t", reftree_path, "-s", os.path.join(workdir, "ref.fasta"),
           "-q", os.path.join(workdir, "query.fasta"), "-m", model,
           "--outdir", workdir, "--redo", "-T", str(threads)] + epang_placement_flags()
    r = subprocess.run(cmd, capture_output=True, text=True)
    jpf = os.path.join(workdir, "epa_result.jplace")
    if r.returncode != 0 or not os.path.isfile(jpf):
        raise RuntimeError("EPA-ng final placement failed: " + (r.stderr[-300:] if r.stderr else "no jplace"))
    _log("%d suspects re-placed (model %s)" % (len(query_seqs), os.path.basename(model)))
    return jpf
