# SATIVA with EPA-ng

This branch replaces SATIVA's placement engine with EPA-ng (Barbera et al. 2019) and leaves
everything else alone. It is the code measured in the benchmark at
https://github.com/Aaramis/sativa-epang-benchmark, whose RESULTS.md carries the speed and
agreement numbers, the ground-truth mislabel test, and the method.

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

`-stage` splits the leave-one-out into three, so a workflow manager can place the folds
itself. `-stage all` is the default and calls the same three functions in one process, so
there is no second code path.

```bash
# 1. build the reference and write one directory per fold
sativa.py -s aln.fasta -t taxonomy.tsv -x BOT -n run -o out -stage loo-tasks

# 2. place every fold. Each fold directory holds ref.nwk, ref.fasta, query.fasta and model
#    and needs nothing outside itself; manifest.json carries the command to run in it.
sativa.py -stage loo-place -taskdir out/run.l1o_tasks -T 8

# 3. map the placements back onto the reference and finish the analysis
sativa.py -r out/run.refjson -n run -o out -stage loo-score -taskdir out/run.l1o_tasks
```

`out/run.l1o_tasks/manifest.json` describes the whole job: the folds, which sequences are
held out in each, the model, and the EPA-ng command to run in a fold directory. Step 3
needs that manifest as well as the jplace files, because an EPA-ng edge number means
nothing outside the fold that produced it.

The three-step run and the one-shot run produce the same `.mis` file, byte for byte,
whatever order the folds are placed in, because the placements are sorted before SATIVA
sees them (`SATIVA_EPANG_SORT`). `tests/roundtrip.sh` in the benchmark repository checks it
three ways, one of which copies every fold to a directory of its own and places it in a
process that has no access to the reference or the other folds.

**On batching.** A batch of placements and a fold are the same thing: two held-out
sequences can only share one EPA-ng call if they are held out together, their references
differing by exactly the leaf under test. `SATIVA_EPANG_FOLDS` therefore trades placement
time against agreement with the strict leave-one-out. On a 5185-taxon reference at `-T 2`,
EPA-ng costs 3.8 s to set the reference up and 59.6 ms per query after that, so a fold of a
few hundred queries spends about a fifth of its time on setup.

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
They exist because `09_optimisations.sh` in the benchmark measures them, and what they cost
is in `OPTIMISATIONS.md`.

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

Python 3 with `ete3`, and `epa-ng` 0.3.8 on `PATH` (or `SATIVA_EPANG_BIN`). A gcc able
to build RAxML 8.2.3. SATIVA is GPL 3, and so is this copy; see `LICENSE`.
