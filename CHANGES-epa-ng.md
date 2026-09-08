# SATIVA with EPA-ng

This branch replaces SATIVA's placement engine with EPA-ng (Barbera et al. 2019) and leaves
everything else alone.

Measured against unmodified SATIVA v0.9.3 on ITS alignments, with the same reference tree in
both: 257 s to 15.4 s at 1600 sequences on 8 threads, 1536 s to 87 s at 5402 on 16. The gain
is parallelism over queries, which RAxML gets little of on a 242 column alignment. Agreement
at 800 sequences with the default 25 folds is 0.98 recall and 0.90 precision, against
0.92 / 0.90 when unmodified SATIVA merely changes its own substitution model. Of 33 mislabels
injected into three clades, unmodified SATIVA finds 28 and this version 33.

## What differs from upstream

Base: `amkozlov/sativa` v0.9.3 (commit `68284c2`), the authors' own python 3 version. The
whole difference is:

| File | Change |
|---|---|
| `sativa.py` | `run_leave_seq_out_test()` places with EPA-ng in K folds instead of RAxML `-f O`. `run_epa_once()` runs the final confirmation on EPA-ng instead of RAxML `-f v`. Both fall back to RAxML with `SATIVA_L1O_ENGINE=raxml`. `-stage` and `-taskdir` run the leave-one-out one step at a time (below). Temporary files go to the output directory when the install directory is read only. |
| `epac/epang_l1o.py` | New. `run_epang_l1o()` (pass 1) and `run_epang_final()` (pass 2), plus the mapping from EPA-ng edge numbers back to SATIVA's `B=` numbering, and the three staged steps `run_epang_l1o()` is built from. |
| `epac/config.py` | `shutil.rmtree(..., ignore_errors=True)` when cleaning the temp directory, which otherwise races on a parallel filesystem. Carries `-stage` and `-taskdir`. |

`epac/classify_util.py`, `epac/taxonomy_util.py`, `epac/json_util.py`,
`epac/raxml_util.py` and `epac/msa.py` are byte identical to upstream. The decision rule,
the branch labelling, the confidence computation and the reference tree step are untouched.

## Running the leave-one-out in steps

`-stage` splits the analysis into steps a workflow manager can schedule separately, each
with its own resources. `-stage all` is the default and calls the same functions in one
process, so there is no second code path.

```bash
# 1. reference tree from the taxonomy. RAxML-bound, and the step worth its own resources.
sativa.py -s aln.fasta -t taxonomy.tsv -x BOT -n run -o out -stage reference

# 2. one directory per fold
sativa.py -r out/run.refjson -n run -o out -stage loo-tasks

# 3. place every fold. Each fold directory holds ref.nwk, ref.fasta, query.fasta and model
#    and needs nothing outside itself; manifest.json carries the command to run in it.
sativa.py -stage loo-place -taskdir out/run.l1o_tasks -T 8

# 4. map the placements back onto the reference and finish the analysis
sativa.py -r out/run.refjson -n run -o out -stage loo-score -taskdir out/run.l1o_tasks
```

Steps 1 and 2 can still be done together: `-stage loo-tasks` without `-r` builds the
reference first, as before.

### A reference tree inferred elsewhere

`-reftree` takes the topology as given instead of running the constrained RAxML search, and
`-refmodel` the model that goes with it, as a RAxML-NG model string or file:

```bash
sativa.py -s aln.fasta -t taxonomy.tsv -x BOT -n run -o out -stage reference \
          -reftree raxml-ng.bestTree -refmodel raxml-ng.bestModel
```

Nothing about the refjson is tied to RAxML. The branch `B=` values are identifiers shared
between the tree and `branch_tax_map`, and that map, the node heights and the speciation
rate are all computed in python from the tree and the taxonomy. The branch numbering itself
comes from placing a dummy query and reading the numbered tree out of the jplace, which
RAxML and EPA-ng write in the same `{n}` convention — so with `-reftree` it is EPA-ng that
does it, and no RAxML runs at all. `binary_model` is then empty; only the RAxML fallback
ever reads it.

Feeding a reference's own tree back through `-reftree` reproduces its `.mis` byte for byte,
which is mode E of `tests/roundtrip.sh`.

`-stage reference` writes `NAME.model` next to `NAME.refjson`. EPA-ng needs the model the
tree was built under, and it otherwise lives in the temp directory the run deletes, so a
reference reused later with `-r` would silently fall back to fitting GTR+G itself: slower,
and not the model the tree was built under. Any run that builds a reference now writes that
file, and any run given `-r` picks it back up.

`out/run.l1o_tasks/manifest.json` describes the whole job: the folds, which sequences are
held out in each, the model, and the EPA-ng command to run in a fold directory. Step 3
needs that manifest as well as the jplace files, because an EPA-ng edge number means
nothing outside the fold that produced it.

The staged run and the one-shot run produce the same `.mis` file, byte for byte, whatever
order the folds are placed in, because the placements are sorted before SATIVA sees them
(`SATIVA_EPANG_SORT`). `tests/roundtrip.sh` checks it four ways: one shot, staged in place,
staged with every fold copied to a directory of its own and placed in a process that has no
access to the reference or the other folds, the four steps run as four separate invocations,
and the reference rebuilt from a supplied tree. Identical `.mis` at 38, 400 and 1600
sequences.

**On batching.** A batch of placements and a fold are the same thing: two held-out
sequences can only share one EPA-ng call if they are held out together, their references
differing by exactly the leaf under test. What that costs is set by the *fraction* of the
reference a fold removes, `1/K`, not by how many sequences are in it: the measured point is
4% removed, giving 0.98 recall and 0.90 precision against unmodified SATIVA. Batches of
10 000 are 1% of a million sequences and 50% of twenty thousand.

Cost of one EPA-ng call on a 5185-taxon reference: 3.8 s to set the reference up plus
59.6 ms per query at `-T 2`, 3.0 s plus 15.7 ms at `-T 8`. Setup barely parallelises and
placement nearly does, so a fold of a few hundred queries spends about a fifth of its time
on setup and one of a few thousand almost none.

Memory is the other bound, and the one that decides the shape at scale: each concurrent
EPA-ng holds its own copy of the reference, about 3 GB at 5400 taxa on a 242-column
alignment. Past a few hundred thousand taxa a single instance is the constraint rather than
the orchestration, and the budget is better spent on threads inside one placement than on
concurrent placements (`SATIVA_EPANG_FOLD_JOBS=1`).

## Version

`sativa.py -version` prints it, and so does the banner every run writes. `SATIVA_EPANG_BUILD`
in `epac/version.py` carries the fork's version and is bumped on every release tag;
`SATIVA_BUILD` stays at the upstream SATIVA release this is based on.

## Environment variables

| Variable | Default | Effect |
|---|---|---|
| `SATIVA_EPANG_FOLDS` | 25 | Number of folds the leave one out is split into. A value at or above the leaf count gives one sequence per fold, the strict leave one out. |
| `SATIVA_L1O_ENGINE` | `epang` | `raxml` reverts both placement passes to RAxML. |
| `SATIVA_EPANG_HEUR` | `on` | `off` passes `--no-heur` to EPA-ng, which then evaluates every branch, as RAxML does below 1000 taxa. |
| `SATIVA_EPANG_BLO` | `sliding` | `raxml` passes `--raxml-blo`, RAxML style branch length optimisation. |
| `SATIVA_EPANG_ACC_LWR` | 0.99999 | Accumulated likelihood weight kept. RAxML uses 0.999. |
| `SATIVA_EPANG_FINAL_MODEL` | from `RAxML_info` | Model for the confirmation pass. |
| `SATIVA_EPANG_BIN` | from `PATH` | EPA-ng binary. |
| `SATIVA_EPANG_DEBUG` | unset | Logs the model EPA-ng reports and how much likelihood weight the edge remapping drops. |
| `SATIVA_EPANG_FOLD_JOBS` | up to 4 | Folds placed at once, sharing the thread budget. |
| `SATIVA_EPANG_EMIT_JOBS` | 4 | Fold directories written at once. Each fold writes its own copy of the reference alignment, which on a network filesystem is latency rather than throughput. |
| `SATIVA_EPANG_SORT` | 1 | Sorts the placements before SATIVA classifies them. This is what makes a run reproducible, and what makes the staged run agree with the one-shot run whatever order the folds come back in. |
| `SATIVA_EPANG_MODEL` | from `RAxML_info` | Model for the leave-one-out. Needed with `-r`, where there is no `RAxML_info` to find. |
| `SATIVA_EPANG_FAST_MAP` | 1 | 0 uses the ete3 bipartition map instead of the linear one. |
| `SATIVA_EPANG_FAST_PRUNE` | 1 | 0 uses ete3's `prune()` instead of the array pass. |
| `SATIVA_EPANG_MAP_CHECK` | unset | Builds both edge maps and logs where they differ. See the note below. |
| `SATIVA_EPANG_FOLD_ORDER` | `name` | `tree` deals the leaves round robin in tree order, so that neighbours land in different folds. |
| `SATIVA_EPANG_DYN_HEUR` | EPA-ng default | `--dyn-heur`. |
| `SATIVA_EPANG_FIX_HEUR` | unset | `--fix-heur`, the same kind of shortcut SATIVA asks RAxML for above 1000 taxa. |
| `SATIVA_EPANG_PRECISION` | 10 | `--precision`, decimals in the jplace. |

Two further approximations of the leave-one-out are implemented, off, and not recommended.
They agree with the exact leave-one-out on 0.54 and 0.73 of the calls respectively, which is
no better than simply lowering the number of folds at the same cost. They are kept because
that measurement is worth being able to repeat.

| Variable | Default | Effect |
|---|---|---|
| `SATIVA_EPANG_SELF_PLACE` | 0 | Places every sequence once on the whole tree and masks its own branches afterwards, instead of running K folds. |
| `SATIVA_EPANG_SELF_MASK` | `neighbour` | `pendant` masks only the query's own branch. |
| `SATIVA_EPANG_SELF_MAX` | 200 | Placements kept per query before masking. |
| `SATIVA_EPANG_SELF_HEUR` | `fix:0.02` | EPA-ng heuristic for that single run. |
| `SATIVA_EPANG_SCREEN` | unset | Keeps the exact leave-one-out but only for sequences a cheap first pass finds suspicious, turning K runs into two. |
| `SATIVA_EPANG_SCREEN_HEIGHT` | 3 | How many ancestors up the neighbourhood reaches. |
| `SATIVA_EPANG_SCREEN_MASK` | `neighbour` | As `SELF_MASK`, for the screening pass. |
| `SATIVA_EPANG_SCREEN_RULE` | `top` | How the first pass decides a sequence is suspicious. |
| `SATIVA_EPANG_SCREEN_FOLD_FRAC` | 0.04 | Fold size for the exact second pass. |

## The two edge maps do not agree everywhere

`SATIVA_EPANG_MAP_CHECK=1` builds the linear map and the ete3 one side by side. They agree
on almost every edge, and where they differ the linear map has one edge more: 2 folds out of
25 at 400 sequences, 1 out of 25 at 1600, one edge each time.

The cause is a bipartition that splits the pruned tree exactly in half. Both maps key a
bipartition by one canonical side, but the ete3 version picks it by set size, and at a
192-192 split that rule does not pick the same side for the table entry as for the lookup,
so the lookup misses an edge whose complement is sitting in the table. The linear map keys
on `min(hash, all ^ hash)`, which is complement-invariant, and finds it.

So the linear map is the more correct of the two, and the ete3 one is kept only as the
reference to check it against. The `.mis` file is the same either way on everything measured
here: a dropped edge loses its likelihood weight, and one edge out of 765 did not move a
confidence far enough to change a call.

## Building RAxML

Build RAxML once before running anything, as upstream does:

```bash
cd raxml && make
```

SATIVA calls them for the reference tree and for the RAxML fallback. In the benchmark the
same binaries are copied into the two upstream checkouts, so a comparison never changes
binary as well as engine.

## Requirements

Python 3, and `epa-ng` 0.3.8 on `PATH` (or `SATIVA_EPANG_BIN`). A gcc able to build RAxML
8.2.3. No third-party python package: the tree parsing goes through the `epac/ete2` that
SATIVA already vendors, as the rest of the code does. SATIVA is GPL 3, and so is this copy;
see `LICENSE`.

There is also a conda package, submitted to bioconda as `sativa-epang`. It installs the
command as **`sativa-epang`** rather than `sativa.py`, because the `sativa` package already
claims that path and the two are meant to coexist:

```bash
conda install -c bioconda -c conda-forge sativa-epang
sativa-epang -s aln.fasta -t taxonomy.tsv -x BOT -n run -o out
```
