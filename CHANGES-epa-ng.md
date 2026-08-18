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
| `sativa.py` | Two hunks. `run_leave_seq_out_test()` places with EPA-ng in K folds instead of RAxML `-f O`. `run_epa_once()` runs the final confirmation on EPA-ng instead of RAxML `-f v`. Both fall back to RAxML with `SATIVA_L1O_ENGINE=raxml`. |
| `epac/epang_l1o.py` | New. `run_epang_l1o()` (pass 1) and `run_epang_final()` (pass 2), plus the mapping from EPA-ng edge numbers back to SATIVA's `B=` numbering. |
| `epac/config.py` | One hunk. `shutil.rmtree(..., ignore_errors=True)` when cleaning the temp directory, which otherwise races on a parallel filesystem. |

`epac/classify_util.py`, `epac/taxonomy_util.py`, `epac/json_util.py`,
`epac/raxml_util.py` and `epac/msa.py` are byte identical to upstream. The decision rule,
the branch labelling, the confidence computation and the reference tree step are untouched.

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

## Checksums of the changed files (benchmark copy)

```
539cc3fdc30f2215f2e267412e076cde  sativa.py
5e403bd9009929e90c38a8b3c28451bb  epac/epang_l1o.py
79d6c16ab7107082d637e51213459cba  epac/config.py
```

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
