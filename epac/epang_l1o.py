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
import os, sys, re, json, glob, subprocess, shutil, time
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

# --- linear-time bipartition mapping -----------------------------------------------
# _bip_map above is the reference for what follows: it asks ete3 for the leaf names under
# every node and takes the complement when that side is the larger one, both O(n) per node,
# so the whole map is O(n^2) and the leave-one-out pays it once per fold. Below, a clade is
# summarised by the XOR of one random 128-bit value per leaf: a node's summary is the XOR of
# its children's, so a tree costs O(n), and the complement of a side is just `all ^ side`,
# which makes the canonical key min(side, all ^ side) free. Two distinct clades collide only
# if 128 random bits collide.
#
# The de-duplication order is reproduced as it stands: inside one tree the later node wins
# (dict overwrite in level order, the order ete3's traverse() uses), and when several
# reference bipartitions collapse onto the same one inside a fold, the first wins.
import random

_LEAF_HASH = {}
_HASH_RNG = random.Random(20160104)

def _leaf_hash(name):
    h = _LEAF_HASH.get(name)
    if h is None:
        h = _HASH_RNG.getrandbits(128)
        _LEAF_HASH[name] = h
    return h

_B_RE = re.compile(r"B=(\d+)")
_STOP = "(),:;[]{}"

def _read_annot(s, i, node, name, tag, tag_style):
    """Consume a node's name, its comments, its {edge} label and its branch length."""
    n = len(s)
    if i < n and s[i] == "'":
        j = s.index("'", i + 1)
        name[node] = s[i + 1:j]
        i = j + 1
    else:
        start = i
        while i < n and s[i] not in _STOP:
            i += 1
        if i > start:
            name[node] = s[start:i]
    while i < n:
        c = s[i]
        if c == "[":
            j = s.index("]", i)
            if tag_style == "B":
                m = _B_RE.search(s, i, j)
                if m:
                    tag[node] = m.group(1)
            i = j + 1
        elif c == "{":
            j = s.index("}", i)
            if tag_style == "EDGE":
                tag[node] = s[i + 1:j]
            i = j + 1
        elif c == ":":
            i += 1
            while i < n and s[i] not in _STOP:
                i += 1
        else:
            break
    return i

def _parse_newick_arrays(s, tag_style):
    """Newick -> (children, name, tag, root): only what the bipartition map needs."""
    children, name, tag = [], [], []
    def _new(parent):
        children.append([]); name.append(None); tag.append(None)
        idx = len(children) - 1
        if parent is not None:
            children[parent].append(idx)
        return idx
    stack, i, n, root = [], 0, len(s), None
    while i < n:
        c = s[i]
        if c == "(":
            node = _new(stack[-1] if stack else None)
            stack.append(node)
            i += 1
        elif c in ", \t\r\n":
            i += 1
        elif c == ";":
            break
        elif c == ")":
            node = stack.pop()
            i = _read_annot(s, i + 1, node, name, tag, tag_style)
            if not stack:
                root = node
        else:
            node = _new(stack[-1] if stack else None)
            i = _read_annot(s, i, node, name, tag, tag_style)
            if not stack:
                root = node
    if root is None:
        raise RuntimeError("epang-l1o: could not parse the newick string")
    return children, name, tag, root

def _tree_orders(children, root):
    """(postorder, levelorder)."""
    level, head = [], 0
    level.append(root)
    while head < len(level):
        level.extend(children[level[head]])
        head += 1
    post, stack = [], [root]
    while stack:
        i = stack.pop()
        post.append(i)
        stack.extend(children[i])
    post.reverse()
    return post, level

def _clade_hashes(children, name, post, keep=None):
    """XOR summary of the leaves under each node, counting only those in `keep`."""
    h = [0] * len(children)
    for i in post:
        kids = children[i]
        if not kids:
            nm = name[i]
            if nm is not None and (keep is None or nm in keep):
                h[i] = _leaf_hash(nm)
        else:
            v = 0
            for c in kids:
                v ^= h[c]
            h[i] = v
    return h

class _RefBipartitions:
    """The refjson tree, parsed once, ready to be restricted to any fold's leaves."""

    def __init__(self, tree_str):
        self.children, self.name, self.tag, self.root = _parse_newick_arrays(tree_str, "B")
        self.post, self.level = _tree_orders(self.children, self.root)
        self.leaves = [self.name[i] for i in self.post
                       if not self.children[i] and self.name[i] is not None]
        full = _clade_hashes(self.children, self.name, self.post)
        allx = full[self.root]
        dedup = {}
        for i in self.level:
            if self.tag[i] is None:
                continue
            key = full[i]
            other = allx ^ key
            dedup[key if key < other else other] = (i, self.tag[i])
        self.tagged = list(dedup.values())

    def leaf_neighbourhood(self):
        """For each leaf, the B= ids that stop existing once that leaf is pruned away.

        Pruning a leaf removes its pendant branch and merges its sister branch with its
        parent branch, so those three ids are the ones a placement of that sequence must
        not be allowed to use.
        """
        parent = [-1] * len(self.children)
        for i, kids in enumerate(self.children):
            for c in kids:
                parent[c] = i
        out = {}
        for i in self.post:
            if self.children[i] or self.name[i] is None:
                continue
            ids = set()
            if self.tag[i] is not None:
                ids.add(self.tag[i])
            p = parent[i]
            if p >= 0:
                if self.tag[p] is not None:
                    ids.add(self.tag[p])
                for c in self.children[p]:
                    if c != i and self.tag[c] is not None:
                        ids.add(self.tag[c])
            out[self.name[i]] = ids
        return out

    def leaf_pendant(self):
        """For each leaf, only its own pendant branch."""
        return {self.name[i]: ({self.tag[i]} if self.tag[i] is not None else set())
                for i in self.post if not self.children[i] and self.name[i] is not None}

    def restrict(self, keep):
        """{bipartition of the pruned tree -> B=}, first reference branch wins."""
        h = _clade_hashes(self.children, self.name, self.post, keep)
        allx = h[self.root]
        out = {}
        for i, b in self.tagged:
            key = h[i]
            other = allx ^ key
            out.setdefault(key if key < other else other, b)
        return out

def _epa_edge_map(tree_str):
    """{bipartition -> EPA-ng edge number} for a jplace tree, later node wins."""
    children, name, tag, root = _parse_newick_arrays(tree_str, "EDGE")
    post, level = _tree_orders(children, root)
    h = _clade_hashes(children, name, post)
    allx = h[root]
    out = {}
    for i in level:
        if tag[i] is None:
            continue
        key = h[i]
        other = allx ^ key
        out[key if key < other else other] = tag[i]
    return out


def _flatten_tree(tree):
    """Freeze an ete3 tree into arrays, so pruning does not need ete3 at all.

    Returns (children, name, dist, postorder). ete3's own prune() is the single most
    expensive Python step of the leave-one-out: it calls get_common_ancestor for every
    kept node, and with K folds that cost is paid K times. The arrays here are built once.
    """
    index = {}
    for node in tree.traverse("postorder"):
        index[id(node)] = len(index)
    size = len(index)
    children = [[] for _ in range(size)]
    name = [""] * size
    dist = [0.0] * size
    postorder = []
    for node in tree.traverse("postorder"):
        i = index[id(node)]
        postorder.append(i)
        name[i] = node.name or ""
        dist[i] = float(node.dist)
        children[i] = [index[id(c)] for c in node.children]
    return children, name, dist, postorder


def _prune_to_newick(flat, keep):
    """Newick of the tree restricted to `keep`, with branch lengths preserved.

    Same contract as ete3's prune(..., preserve_branch_length=True): a node left with a
    single child is suppressed and its branch length added to that child. O(n) per fold.
    """
    children, name, dist, postorder = flat
    keep = set(keep)
    piece = [None] * len(name)     # Newick of the kept subtree below each node, or None
    extra = [0.0] * len(name)      # branch length inherited from suppressed ancestors

    for i in postorder:
        kids = [c for c in children[i] if piece[c] is not None]
        if not children[i]:
            piece[i] = name[i] if name[i] in keep else None
            continue
        if not kids:
            piece[i] = None
        elif len(kids) == 1:
            c = kids[0]
            piece[i] = piece[c]
            extra[i] = extra[c] + dist[c]      # collapse: the child carries both branches
        else:
            parts = []
            for c in kids:
                # ete3's writer formats branch lengths with %0.6g; matching it exactly is
                # what makes the pruned tree the same tree, not merely the same topology.
                # %0.6f, six decimals rather than six significant digits, silently rounded
                # short branches and moved the lengths by up to 5e-7.
                parts.append("%s:%0.6g" % (piece[c], dist[c] + extra[c]))
            piece[i] = "(" + ",".join(parts) + ")"
            extra[i] = 0.0

    root = postorder[-1]
    if piece[root] is None:
        raise RuntimeError("epang-l1o: pruning removed every leaf")
    body = piece[root]
    return (body if body.startswith("(") else "(%s)" % body) + ";"


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
#
# Two more knobs go the other way, towards less work rather than more:
#
#   SATIVA_EPANG_DYN_HEUR=<f>  -> --dyn-heur <f>   thorough insertion on the branches
#                                 accumulating <f> of the preplacement weight (EPA-ng
#                                 default 0.99999)
#   SATIVA_EPANG_FIX_HEUR=<f>  -> --fix-heur <f>   thorough insertion on that share of the
#                                 branches, by preplacement rank. This is the same kind of
#                                 shortcut SATIVA itself asks RAxML for above 1000 taxa
#                                 (-G 500/n: 31% of the branches at n=1600, 9% at n=5402).
#   SATIVA_EPANG_PRECISION=<n> -> --precision <n>  decimals in the jplace (default 10); the
#                                 file is read back by json.load, so its size is our cost.
def epang_placement_flags():
    flags = ["--filter-acc-lwr", os.environ.get("SATIVA_EPANG_ACC_LWR", "0.99999"),
             "--filter-max", "100000"]
    fix_heur = os.environ.get("SATIVA_EPANG_FIX_HEUR")
    dyn_heur = os.environ.get("SATIVA_EPANG_DYN_HEUR")
    if os.environ.get("SATIVA_EPANG_HEUR", "on").lower() in ("off", "no", "0", "false"):
        flags.append("--no-heur")
    elif fix_heur:
        flags += ["--fix-heur", fix_heur]
    elif dyn_heur:
        flags += ["--dyn-heur", dyn_heur]
    if os.environ.get("SATIVA_EPANG_BLO", "sliding").lower() == "raxml":
        flags.append("--raxml-blo")
    precision = os.environ.get("SATIVA_EPANG_PRECISION")
    if precision:
        flags += ["--precision", precision]
    return flags


def _slow_edge_map(refjson_tree_str, epa_tree_str, keep):
    """The ete3 remap, kept as the reference SATIVA_EPANG_MAP_CHECK compares against."""
    bidmap, _ = _bip_map(refjson_tree_str, "B")
    e2side, foldL = _bip_map(epa_tree_str, "EDGE")
    keepS = frozenset(keep)
    restricted = {}
    for side in bidmap:
        s2 = side & keepS
        key = s2 if len(s2) <= len(keepS) - len(s2) else (keepS - s2)
        restricted.setdefault(frozenset(key), bidmap[side])
    out = {}
    for side, e in e2side.items():
        key = side if len(side) <= len(foldL) - len(side) else (foldL - side)
        b = restricted.get(frozenset(key))
        if b is not None:
            out[int(e)] = b
    return out


def _run_self_place(refbip, leaves, aln_by_leaf, reftree_path, model, workdir, threads, _log,
                    return_fields=False, mask_mode=None):
    """One EPA-ng run for the whole set, with each query's own branches masked afterwards.

    The accumulated-LWR filter has to be opened up here: a sequence placed on a tree that
    still contains it puts essentially all its weight on its own pendant branch, and the
    default filter would cut away exactly the alternatives the masking then needs.
    """
    mask_mode = (mask_mode or os.environ.get("SATIVA_EPANG_SELF_MASK", "neighbour")).lower()
    keep_max = os.environ.get("SATIVA_EPANG_SELF_MAX", "200")
    os.makedirs(workdir, exist_ok=True)
    ref_path = os.path.join(workdir, "ref.fasta")
    qry_path = os.path.join(workdir, "query.fasta")
    _write_fasta({l: aln_by_leaf[l] for l in leaves}, ref_path)
    _write_fasta({"q_%s" % l: aln_by_leaf[l] for l in leaves}, qry_path)

    flags = ["--filter-acc-lwr", os.environ.get("SATIVA_EPANG_ACC_LWR", "1.0"),
             "--filter-max", keep_max]
    # EPA-ng's default heuristic picks the branches to evaluate thoroughly by accumulated
    # preplacement weight. A sequence still present in the tree takes all of that weight on
    # its own branch, so the default would hand back only the branches the mask is about to
    # remove. The candidate set has to be chosen by rank instead: --fix-heur takes a fixed
    # share of the branches, which is also what SATIVA asks RAxML for on large trees.
    heur = os.environ.get("SATIVA_EPANG_SELF_HEUR", "fix:0.02").lower()
    if heur.startswith("fix"):
        flags += ["--fix-heur", heur.split(":", 1)[1] if ":" in heur else "0.02"]
    elif heur.startswith("baseball"):
        flags.append("--baseball-heur")
    elif heur in ("no", "off", "none"):
        flags.append("--no-heur")
    if os.environ.get("SATIVA_EPANG_BLO", "sliding").lower() == "raxml":
        flags.append("--raxml-blo")
    cmd = [EPANG, "-t", reftree_path, "-s", ref_path, "-q", qry_path, "-m", model,
           "--outdir", workdir, "--redo", "-T", str(threads)] + flags
    r = subprocess.run(cmd, capture_output=True, text=True)
    jpf = os.path.join(workdir, "epa_result.jplace")
    if r.returncode != 0 or not os.path.isfile(jpf):
        raise RuntimeError("EPA-ng self-placement failed: " +
                           (r.stderr[-300:] if r.stderr else "no jplace"))
    with open(jpf) as handle:
        d = json.load(handle)
    ie = d["fields"].index("edge_num")
    ilwr = d["fields"].index("like_weight_ratio")

    restricted = refbip.restrict(leaves)
    epa2b = {}
    for key, e in _epa_edge_map(d["tree"]).items():
        b = restricted.get(key)
        if b is not None:
            epa2b[int(e)] = b
    mask = refbip.leaf_pendant() if mask_mode == "pendant" else refbip.leaf_neighbourhood()

    placements, emptied = [], 0
    for pl in d["placements"]:
        name = (pl.get("n") or pl.get("nm"))[0]
        if isinstance(name, list):
            name = name[0]
        leaf = name[2:] if name.startswith("q_") else name
        drop = mask.get(leaf, ())
        rows, total = [], 0.0
        for row in pl["p"]:
            b = epa2b.get(int(row[ie]))
            if b is None or b in drop:
                continue
            rr = list(row)
            rr[ie] = int(b)
            rows.append(rr)
            total += float(rr[ilwr])
        if not rows:
            emptied += 1
            continue
        if total > 0:
            for rr in rows:
                rr[ilwr] = float(rr[ilwr]) / total
        placements.append({"p": rows, "n": [leaf]})
    _log("self-placement: %d sequences, mask=%s, %d left with no usable branch"
         % (len(placements), mask_mode, emptied))
    return (placements, ie, ilwr) if return_fields else placements


def _subtree_intervals(refbip):
    """DFS numbering of the reference tree, so 'is this branch inside that clade' is O(1)."""
    children, root = refbip.children, refbip.root
    tin = [0] * len(children)
    tout = [0] * len(children)
    parent = [-1] * len(children)
    counter = 0
    stack = [(root, False)]
    while stack:
        i, done = stack.pop()
        if done:
            tout[i] = counter
            continue
        counter += 1
        tin[i] = counter
        stack.append((i, True))
        for c in children[i]:
            parent[c] = i
            stack.append((c, False))
    for i in refbip.post:                      # tout as the max tin below the node
        if children[i]:
            tout[i] = max(tout[c] for c in children[i])
        else:
            tout[i] = tin[i]
    return tin, tout, parent


def _run_screened(refbip, leaves, aln_by_leaf, reftree_path, model, workdir, threads,
                  _log, threshold, height):
    """Two passes: a cheap one over everything, an exact one over what it finds suspicious.

    A sequence that is correctly labelled sits, in a taxonomy-constrained reference tree,
    exactly where free placement puts it back. So a first self-placement pass can tell the
    sequences whose weight stays in their own corner of the tree from the ones whose weight
    goes elsewhere, and only the latter need the expensive leave-one-out. The second pass
    removes ALL of them from the reference at once, which is one EPA-ng run rather than K.

    SATIVA_EPANG_SCREEN=<mass>   below this weight in its own neighbourhood, a sequence is
                                 taken as suspicious (default 0.9)
    SATIVA_EPANG_SCREEN_HEIGHT=<n>  how many ancestors up that neighbourhood reaches (3)
    """
    # Pass 1 masks the pendant branch only. Masking the sister branch as well, which is what
    # the standalone self-placement does, pushes the weight of a perfectly normal sequence
    # away from its own corner too, and that is exactly the signal the screen reads.
    first, ie, ilwr = _run_self_place(refbip, leaves, aln_by_leaf, reftree_path, model,
                                      os.path.join(workdir, "screen_pass1"), threads, _log,
                                      return_fields=True,
                                      mask_mode=os.environ.get("SATIVA_EPANG_SCREEN_MASK",
                                                               "pendant"))

    tin, tout, parent = _subtree_intervals(refbip)
    b2node = {}
    for i, b in refbip.tagged:
        b2node[int(b)] = i
    leaf_node = {refbip.name[i]: i for i in refbip.post
                 if not refbip.children[i] and refbip.name[i] is not None}

    # A sequence's own neighbourhood is the smallest ancestor clade holding at least `height`
    # leaves. Counting ancestors instead would not work: SATIVA resolves the multifurcations
    # of the taxonomy arbitrarily, so a genus of fifty can come out as a ladder in which ten
    # ancestors still only cover eleven of its members.
    nleaves = [0] * len(refbip.children)
    for i in refbip.post:
        kids = refbip.children[i]
        nleaves[i] = sum(nleaves[c] for c in kids) if kids else 1
    local_root = {}
    for name, i in leaf_node.items():
        a = i
        while nleaves[a] < height and parent[a] >= 0:
            a = parent[a]
        local_root[name] = a

    # Two ways to read the first pass. "mass" asks how much of the weight stayed in the
    # sequence's own neighbourhood; "top" only asks whether the single best branch is in it,
    # which does not depend on how flat EPA-ng's weights happen to be for that query.
    rule = os.environ.get("SATIVA_EPANG_SCREEN_RULE", "top").lower()
    suspicious, by_name, masses = [], {}, []
    for pl in first:
        name = pl["n"][0]
        by_name[name] = pl
        a = local_root.get(name)
        mass, best, best_lwr = 0.0, None, -1.0
        if a is not None:
            lo, hi = tin[a], tout[a]
            for row in pl["p"]:
                node = b2node.get(int(row[ie]))
                inside = node is not None and lo <= tin[node] <= hi
                w = float(row[ilwr])
                if not 0.0 <= w <= 1.0:      # a query EPA-ng could not weight; take no view
                    continue
                if inside:
                    mass += w
                if w > best_lwr:
                    best_lwr, best = w, inside
        masses.append(mass)
        if (not best) if rule == "top" else (mass < threshold):
            suspicious.append(name)
    if masses:
        masses.sort()
        _log("screen: weight kept in the sequence's own neighbourhood, deciles %s"
             % " ".join("%.2f" % masses[int(q * (len(masses) - 1) / 10)] for q in range(11)))
    suspicious.extend(l for l in leaves if l not in by_name)   # nothing placed at all

    _log("screen: %d of %d sequences go to the exact pass (%.1f%%)"
         % (len(suspicious), len(leaves), 100.0 * len(suspicious) / max(1, len(leaves))))
    susp = set(suspicious)
    if not suspicious or len(leaves) - len(susp) < 4:
        return first

    # Pass 2 is the real leave-one-out, but only over the suspicious sequences. It is split
    # so that no run takes more than `frac` of the taxa out of the reference: removing them
    # all at once would leave a reference too thin for the placements to mean anything.
    frac = float(os.environ.get("SATIVA_EPANG_SCREEN_FOLD_FRAC", "0.04"))
    per_fold = max(1, int(frac * len(leaves)))
    m = max(1, (len(suspicious) + per_fold - 1) // per_fold)
    suspicious.sort()
    sub_folds = [suspicious[i::m] for i in range(m)]
    _log("screen: exact pass over %d sequences in %d run(s), each taking %.1f%% of the "
         "taxa out of the reference" % (len(suspicious), m, 100.0 * len(sub_folds[0]) / len(leaves)))

    flat = _flatten_tree(Tree(reftree_path, format=1))
    out = [by_name[l] for l in leaves if l in by_name and l not in susp]
    exact = 0
    fold_jobs = max(1, int(os.environ.get("SATIVA_EPANG_FOLD_JOBS", "1")))
    fold_threads = max(1, threads // fold_jobs)

    def exact_run(item):
        fi, fold = item
        if not fold:
            return None
        wd = os.path.join(workdir, "screen_pass2_%03d" % fi)
        os.makedirs(wd, exist_ok=True)
        fold_set = set(fold)
        ref_leaves = [l for l in leaves if l not in fold_set]
        with open(os.path.join(wd, "ref.nwk"), "w") as handle:
            handle.write(_prune_to_newick(flat, ref_leaves) + "\n")
        _write_fasta({l: aln_by_leaf[l] for l in ref_leaves}, os.path.join(wd, "ref.fasta"))
        _write_fasta({l: aln_by_leaf[l] for l in fold}, os.path.join(wd, "query.fasta"))
        cmd = [EPANG, "-t", os.path.join(wd, "ref.nwk"), "-s", os.path.join(wd, "ref.fasta"),
               "-q", os.path.join(wd, "query.fasta"), "-m", model, "--outdir", wd, "--redo",
               "-T", str(fold_threads)] + epang_placement_flags()
        return fi, wd, ref_leaves, subprocess.run(cmd, capture_output=True, text=True)

    if fold_jobs > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=fold_jobs) as pool:
            runs = list(pool.map(exact_run, enumerate(sub_folds)))
    else:
        runs = [exact_run(item) for item in enumerate(sub_folds)]

    for outcome in runs:
        if outcome is None:
            continue
        fi, wd, ref_leaves, r = outcome
        jpf = os.path.join(wd, "epa_result.jplace")
        if r.returncode != 0 or not os.path.isfile(jpf):
            _log("EPA-ng FAIL on exact run %d: %s" % (fi, (r.stderr or "")[-300:]))
            continue
        with open(jpf) as handle:
            d = json.load(handle)
        je = d["fields"].index("edge_num")
        restricted = refbip.restrict(ref_leaves)
        epa2b = {}
        for key, e in _epa_edge_map(d["tree"]).items():
            b = restricted.get(key)
            if b is not None:
                epa2b[int(e)] = b
        for pl in d["placements"]:
            name = (pl.get("n") or pl.get("nm"))[0]
            if isinstance(name, list):
                name = name[0]
            rows = []
            for row in pl["p"]:
                b = epa2b.get(int(row[je]))
                if b is None:
                    continue
                rr = list(row)
                rr[je] = int(b)
                rows.append(rr)
            if rows:
                out.append({"p": rows, "n": [name]})
                exact += 1
    _log("screen: %d sequences kept from the cheap pass, %d re-placed exactly"
         % (len(out) - exact, exact))
    return out


# --- the leave-one-out in three steps -------------------------------------------------
# Dealing the sequences into K folds and writing what each fold needs, placing every fold
# with EPA-ng, and mapping the placements back onto the B= edge numbering. Each step reads
# its input from a directory rather than from the previous call, so the middle one can be
# run elsewhere:
#
#   emit_l1o_tasks()     taskdir/manifest.json, taskdir/fold_XXX/{ref.nwk,ref.fasta,query.fasta}
#   place_l1o_tasks()    taskdir/fold_XXX/epa_result.jplace, or run manifest["command"]
#   collect_l1o_tasks()  the placement list SATIVA classifies
#
# sativa.py exposes them as -stage {loo-tasks,loo-place,loo-score}. run_epang_l1o() below
# calls the three in a row, so there is no second code path.

MANIFEST_VERSION = 1


def _l1o_setup(refjson_tree_str, refaln_path, raxml_outdir, _log):
    """Model, reference bipartitions and the alignment, keyed by leaf name."""
    # EPA-ng model: RAxML_info.mfresolv if present, otherwise GTR+G (EPA-ng re-evaluates)
    # SATIVA_EPANG_MODEL points EPA-ng at a model explicitly. It is needed when SATIVA is
    # started from a ready-made reference (-r): there is no RAxML_info in the working
    # directory then, and EPA-ng would silently fall back to fitting GTR+G itself, once per
    # fold, which is both slower and a different model from the one the tree was built under.
    info = glob.glob(os.path.join(raxml_outdir, "RAxML_info.mfresolv*"))
    model = os.environ.get("SATIVA_EPANG_MODEL") or (info[0] if info else "GTR+G")
    _log("EPA-ng model: %s" % model)

    # SATIVA_EPANG_FAST_MAP=0 goes back to the ete3 bipartition map, which is what the
    # linear one is checked against.
    fast_map = os.environ.get("SATIVA_EPANG_FAST_MAP", "1").lower() not in ("0", "off", "false")
    if fast_map:
        refbip = _RefBipartitions(refjson_tree_str)
        bidmap, leaves = None, refbip.leaves
    else:
        refbip = None
        bidmap, allL = _bip_map(refjson_tree_str, "B")    # bipartition -> B-id
        leaves = list(allL)
    aln = _read_fasta(refaln_path)
    # align the alignment keys onto the leaf names
    aln_by_leaf = {}
    for lf in leaves:
        for cand in (lf, lf[2:] if lf.startswith("r_") else "r_"+lf):
            if cand in aln:
                aln_by_leaf[lf] = aln[cand]; break
    missing = [l for l in leaves if l not in aln_by_leaf]
    if missing:
        raise RuntimeError("epang-l1o: %d leaves without a sequence (e.g. %s)" % (len(missing), missing[:3]))

    # sorted, so the folds and the per-fold alignments do not depend on the order the
    # leaves came out of the tree
    return {"model": model, "fast_map": fast_map, "refbip": refbip, "bidmap": bidmap,
            "leaves": sorted(leaves), "aln_by_leaf": aln_by_leaf}


def _l1o_folds(leaves, reftree_path, folds, _log):
    """Deal the leaves into K folds. Returns the list of folds, in fold order."""
    # How the folds are made up. By sequence name, the folds are random with respect to the
    # tree, so two sequences of the same species regularly leave the reference together and
    # neither can find the other on the way back: that is the k-fold approximation showing.
    # Dealing the leaves round robin in tree order instead guarantees that neighbours end up
    # in different folds, at no cost. SATIVA_EPANG_FOLD_ORDER=tree turns it on.
    ordered = leaves
    if os.environ.get("SATIVA_EPANG_FOLD_ORDER", "name").lower() == "tree":
        with open(reftree_path) as handle:
            ch, nm, _tg, rt = _parse_newick_arrays(handle.read(), None)
        post, _lvl = _tree_orders(ch, rt)
        known = set(leaves)
        in_tree_order = [nm[i] for i in post if not ch[i] and nm[i] in known]
        if len(in_tree_order) == len(leaves):
            ordered = in_tree_order
        else:
            _log("fold order: the tree gave %d of %d leaf names, keeping name order"
                 % (len(in_tree_order), len(leaves)))

    K = min(folds, len(ordered))
    return [ordered[i::K] for i in range(K)]


def emit_l1o_tasks(refjson_tree_str, refaln_path, reftree_path, raxml_outdir, taskdir,
                   folds=5, log=None):
    """Step 1. Write one self-contained directory per fold, plus manifest.json.

    A fold directory holds the reference tree with that fold's leaves pruned away, the
    matching reference alignment, the fold's queries and the model: everything one EPA-ng
    call needs and nothing outside the directory. Returns the manifest.
    """
    def _log(m):
        if log: log.info("[epang-l1o] " + m)
        else: sys.stderr.write("[epang-l1o] " + m + "\n")

    t_start = time.time()
    setup = _l1o_setup(refjson_tree_str, refaln_path, raxml_outdir, _log)
    leaves, aln_by_leaf = setup["leaves"], setup["aln_by_leaf"]
    folds_list = _l1o_folds(leaves, reftree_path, folds, _log)

    os.makedirs(taskdir, exist_ok=True)

    # The model file has to travel with the tasks, and into every fold directory rather
    # than once at the top, since a scheduler stages one directory at a time. It is a 3 kB
    # RAxML_info. A model name (GTR+G) needs no file. The extra copy at the top is for
    # -stage loo-score: the confirmation pass wants the same model, and by then the temp
    # directory the reference was built in is gone.
    model = setup["model"]
    model_is_file = os.path.isfile(model)
    model_arg = "model" if model_is_file else model
    if model_is_file:
        shutil.copyfile(model, os.path.join(taskdir, "model"))

    full_tree = Tree(reftree_path, format=1)
    # copy(method="newick") serialises and reparses the tree on every fold. The string is
    # the same every time, so build it once; the flattened form feeds the fast prune.
    full_newick = full_tree.write(format=1)
    fast_prune = os.environ.get("SATIVA_EPANG_FAST_PRUNE", "1").lower() not in ("0", "off", "false")
    flat = _flatten_tree(full_tree) if fast_prune else None

    # The K folds write K copies of the reference alignment, 31 MB at 5402 sequences in 3K
    # small files. On a network filesystem that is latency rather than throughput, so the
    # folds are written a few at a time; fold_records is rebuilt in fold order below.
    def write_fold(item):
        fi, fold = item
        fold_set = set(fold)
        ref_leaves = [l for l in leaves if l not in fold_set]
        # fewer than four reference leaves leaves no tree to place into: skip, do not fail
        if len(ref_leaves) < 4 or not fold:
            return None
        name = "fold_%03d" % fi
        wd = os.path.join(taskdir, name)
        os.makedirs(wd, exist_ok=True)
        if fast_prune:
            with open(os.path.join(wd, "ref.nwk"), "w") as handle:
                handle.write(_prune_to_newick(flat, ref_leaves) + "\n")
        else:
            tpr = Tree(full_newick, format=1)
            tpr.prune(ref_leaves, preserve_branch_length=True)
            tpr.write(outfile=os.path.join(wd, "ref.nwk"), format=5)
        _write_fasta({l: aln_by_leaf[l] for l in ref_leaves}, os.path.join(wd, "ref.fasta"))
        _write_fasta({l: aln_by_leaf[l] for l in fold},       os.path.join(wd, "query.fasta"))
        if model_is_file:
            shutil.copyfile(model, os.path.join(wd, "model"))
        return {"id": fi, "dir": name, "n_ref": len(ref_leaves), "queries": list(fold)}

    write_jobs = max(1, int(os.environ.get("SATIVA_EPANG_EMIT_JOBS", "4")))
    if write_jobs > 1 and len(folds_list) > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=write_jobs) as pool:
            written = list(pool.map(write_fold, enumerate(folds_list)))
    else:
        written = [write_fold(item) for item in enumerate(folds_list)]
    fold_records = [record for record in written if record is not None]

    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "sativa_epang": "leave-one-out tasks",
        "n_folds": len(fold_records),
        "n_folds_planned": len(folds_list),
        "n_leaves": len(leaves),
        "model": model_arg,
        "epang_args": epang_placement_flags(),
        # every path relative to the fold directory, every file in it, so this is the
        # command wherever the directory is staged. --redo for retries in place; add -T
        # according to how many threads the caller wants to give one placement.
        "command": ["epa-ng", "-t", "ref.nwk", "-s", "ref.fasta", "-q", "query.fasta",
                    "-m", model_arg, "--outdir", ".", "--redo"] + epang_placement_flags(),
        "output": "epa_result.jplace",
        "leaves": leaves,
        "folds": fold_records,
    }
    with open(os.path.join(taskdir, "manifest.json"), "w") as handle:
        json.dump(manifest, handle, indent=1)
        handle.write("\n")

    _log("wrote %d fold tasks to %s (%.1fs)"
         % (len(fold_records), taskdir, time.time() - t_start))
    return manifest


def read_l1o_manifest(taskdir):
    with open(os.path.join(taskdir, "manifest.json")) as handle:
        manifest = json.load(handle)
    got = manifest.get("manifest_version")
    if got != MANIFEST_VERSION:
        raise RuntimeError("epang-l1o: %s was written by manifest version %s, this SATIVA "
                           "reads version %d" % (taskdir, got, MANIFEST_VERSION))
    return manifest


def place_l1o_tasks(taskdir, threads=1, jobs=None, log=None):
    """Step 2. Run EPA-ng once per fold directory, and return how many were placed.

    Convenience: this is the step a workflow manager takes over, by running
    manifest["command"] in each fold directory itself.
    """
    def _log(m):
        if log: log.info("[epang-l1o] " + m)
        else: sys.stderr.write("[epang-l1o] " + m + "\n")

    t_start = time.time()
    manifest = read_l1o_manifest(taskdir)
    model = manifest["model"]
    records = manifest["folds"]

    # Folds are independent: each has its own tree, its own alignment and its own EPA-ng.
    # SATIVA_EPANG_FOLD_JOBS runs that many at once and splits the thread budget between
    # them, which pays because EPA-ng's own scaling flattens well before the core count.
    # It does not change the result: same folds, placements collected in fold order.
    #
    # Default: up to four folds at once. The cap is there because each concurrent EPA-ng
    # holds its own copy of the reference (about 3 GB at 5400 taxa), and past four the
    # memory bandwidth costs more than the parallelism returns. Raise it on a large node,
    # set it to 1 to go back to one fold at a time.
    if jobs is None:
        jobs = int(os.environ.get("SATIVA_EPANG_FOLD_JOBS",
                                  str(min(max(1, len(records)), max(1, threads), 4))))
    jobs = max(1, jobs)
    fold_threads = max(1, threads // jobs)

    def place(record):
        wd = os.path.join(taskdir, record["dir"])
        # manifest["command"] with the paths made absolute, so this runs from anywhere.
        model_arg = os.path.join(wd, "model") if model == "model" else model
        cmd = [EPANG, "-t", os.path.join(wd, "ref.nwk"), "-s", os.path.join(wd, "ref.fasta"),
               "-q", os.path.join(wd, "query.fasta"), "-m", model_arg,
               "--outdir", wd, "--redo", "-T", str(fold_threads)] + manifest["epang_args"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        with open(os.path.join(wd, "epang.log"), "w") as handle:
            handle.write(proc.stdout or "")
            handle.write(proc.stderr or "")
        if proc.returncode != 0 or not os.path.isfile(os.path.join(wd, "epa_result.jplace")):
            _log("EPA-ng FAIL fold %d: %s" % (record["id"], (proc.stderr or "")[-300:]))
            return 0
        if os.environ.get("SATIVA_EPANG_DEBUG"):
            # What EPA-ng makes of the model file it was handed. Above 500 taxa SATIVA
            # builds the reference tree under GTRCAT, and a CAT RAxML_info carries
            # "alpha: 1.000000" -- a placeholder, since CAT fits no gamma shape.
            for line in (proc.stdout or "").splitlines():
                if any(k in line.lower() for k in ("model", "alpha", "rate")):
                    _log("epa-ng says: " + line.strip())
        return 1

    if jobs > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            placed = sum(pool.map(place, records))
    else:
        placed = sum(place(record) for record in records)

    _log("placed %d/%d folds (%d at a time, %d threads each, %.1fs)"
         % (placed, len(records), jobs, fold_threads, time.time() - t_start))
    return placed


def collect_l1o_tasks(taskdir, refjson_tree_str, log=None):
    """Step 3. Map every fold's EPA-ng edges back onto SATIVA's B= numbering.

    Reads the jplace files whatever produced them and returns the placement list SATIVA
    classifies. The mapping is per fold, an edge number meaning nothing outside the fold
    that produced it, hence the manifest and not just the jplace files.
    """
    def _log(m):
        if log: log.info("[epang-l1o] " + m)
        else: sys.stderr.write("[epang-l1o] " + m + "\n")

    t_remap = time.time()
    manifest = read_l1o_manifest(taskdir)
    leaves = manifest["leaves"]

    fast_map = os.environ.get("SATIVA_EPANG_FAST_MAP", "1").lower() not in ("0", "off", "false")
    if fast_map:
        refbip = _RefBipartitions(refjson_tree_str)
        bidmap = None
    else:
        refbip = None
        bidmap, _allL = _bip_map(refjson_tree_str, "B")

    placements = []
    ie = ilwr = None
    for record in manifest["folds"]:
        fi = record["id"]
        wd = os.path.join(taskdir, record["dir"])
        jpf = os.path.join(wd, manifest["output"])
        if not os.path.isfile(jpf):
            _log("no placement for fold %d (%s)" % (fi, jpf)); continue
        d = json.load(open(jpf))
        ie = d["fields"].index("edge_num")
        epa2b = {}
        fold_set = set(record["queries"])
        if fast_map:
            restricted = refbip.restrict([l for l in leaves if l not in fold_set])
            for key, e in _epa_edge_map(d["tree"]).items():
                b = restricted.get(key)
                if b is not None:
                    epa2b[int(e)] = b
            if os.environ.get("SATIVA_EPANG_MAP_CHECK"):
                # runs the ete3 map alongside and reports where the two disagree
                slow_map = _slow_edge_map(refjson_tree_str, d["tree"],
                                          [l for l in leaves if l not in fold_set])
                diff = [(e, epa2b.get(e), slow_map.get(e))
                        for e in set(epa2b) | set(slow_map) if epa2b.get(e) != slow_map.get(e)]
                _log("map check fold %d: %d edges fast, %d slow, %d disagree %s"
                     % (fi, len(epa2b), len(slow_map), len(diff), diff[:5]))
        else:
            e2side, foldL = _bip_map(d["tree"], "EDGE")
            # table bipartition-restreinte -> B
            restricted = {}
            for side, b in ((s, bidmap[s]) for s in bidmap):
                s2 = side & foldL
                key = s2 if len(s2) <= len(foldL)-len(s2) else (foldL - s2)
                restricted.setdefault(frozenset(key), b)
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
    # EPA-ng above one thread returns the queries in whatever order its threads finish, and
    # at the accumulated-LWR boundary it occasionally keeps one placement more or one less.
    # Downstream that shows up as a different proposed label whenever two candidates tie.
    # Sorting here costs nothing, makes the run reproducible, and makes the staged path
    # agree with the single-process one whatever order the folds come back in.
    if os.environ.get("SATIVA_EPANG_SORT", "1").lower() not in ("0", "off", "false") \
            and placements:
        for pl in placements:
            pl["p"].sort(key=lambda r: (-float(r[ilwr]), int(r[ie])) if ilwr is not None
                         else int(r[ie]))
        placements.sort(key=lambda p: p["n"][0])

    _log("placements produits: %d (K=%d folds, remap %.1fs)"
         % (len(placements), len(manifest["folds"]), time.time() - t_remap))
    return placements


def run_epang_l1o(refjson_tree_str, refaln_path, reftree_path, raxml_outdir,
                  workdir, folds=5, threads=1, log=None):
    t_start = time.time()

    def _log(m):
        if log: log.info("[epang-l1o] " + m)
        else: sys.stderr.write("[epang-l1o] " + m + "\n")

    # SATIVA_EPANG_SELF_PLACE=1 drops the folds entirely: every sequence is placed once on
    # the complete reference tree, and the branches that would have vanished with its own
    # leaf are then struck from its placement list and the remaining weights renormalised.
    # One reference setup instead of K is where the time goes. It is an approximation of the
    # leave-one-out, not a reformulation of it: what it costs in agreement is in RESULTS.md.
    #   SATIVA_EPANG_SELF_MASK=neighbour  pendant + sister + parent branch (default)
    #   SATIVA_EPANG_SELF_MASK=pendant    the query's own pendant branch only
    #   SATIVA_EPANG_SELF_MAX=<N>         placements kept per query before masking
    #
    # SATIVA_EPANG_SCREEN keeps the exact leave-one-out but only for the sequences a cheap
    # first pass finds suspicious, which turns K placement runs into two.
    #
    # Neither has folds, so neither works through the staged entry points.
    self_place = os.environ.get("SATIVA_EPANG_SELF_PLACE", "0").lower() in ("1", "on", "true", "yes")
    screen = os.environ.get("SATIVA_EPANG_SCREEN")
    if self_place or screen:
        setup = _l1o_setup(refjson_tree_str, refaln_path, raxml_outdir, _log)
        refbip = setup["refbip"] or _RefBipartitions(refjson_tree_str)
        if self_place:
            return _run_self_place(refbip, setup["leaves"], setup["aln_by_leaf"],
                                   reftree_path, setup["model"], workdir, threads, _log)
        return _run_screened(refbip, setup["leaves"], setup["aln_by_leaf"], reftree_path,
                             setup["model"], workdir, threads, _log, float(screen),
                             int(os.environ.get("SATIVA_EPANG_SCREEN_HEIGHT", "3")))

    # the three staged steps, one after the other, in this process
    t_prep = time.time()
    emit_l1o_tasks(refjson_tree_str, refaln_path, reftree_path, raxml_outdir, workdir,
                   folds=folds, log=log)
    t_place = time.time()
    place_l1o_tasks(workdir, threads=threads, log=log)
    t_collect = time.time()
    placements = collect_l1o_tasks(workdir, refjson_tree_str, log=log)

    _log("timing: fold prep %.1fs, epa-ng %.1fs, remap %.1fs, elapsed %.1fs"
         % (t_place - t_prep, t_collect - t_place, time.time() - t_collect,
            time.time() - t_start))
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
    # SATIVA_EPANG_MODEL points EPA-ng at a model explicitly. It is needed when SATIVA is
    # started from a ready-made reference (-r): there is no RAxML_info in the working
    # directory then, and EPA-ng would silently fall back to fitting GTR+G itself, once per
    # fold, which is both slower and a different model from the one the tree was built under.
    info = glob.glob(os.path.join(raxml_outdir, "RAxML_info.mfresolv*"))
    model = os.environ.get("SATIVA_EPANG_MODEL") or (info[0] if info else "GTR+G")
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
