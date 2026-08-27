#!/usr/bin/env bash
# Does splitting the leave-one-out into three steps change the answer?
#
# Runs the same analysis three ways and compares the .mis files byte for byte:
#
#   A  one shot            sativa.py -s ... -t ...
#   B  staged, in place    -stage loo-tasks, -stage loo-place, -stage loo-score
#   C  staged, detached    -stage loo-tasks, then every fold copied to a scratch
#                          directory of its own and placed there by manifest["command"]
#                          with nothing else in scope, then -stage loo-score
#   D  four steps          -stage reference, loo-tasks, loo-place, loo-score, each from a
#                          separate invocation, as four workflow processes would be
#
# C is what a workflow manager does: each fold placed by a process with no access to the
# reference, the taxonomy or the other folds, and only the jplace coming back.
#
# Usage: tests/roundtrip.sh ALIGNMENT TAXONOMY TAXCODE [WORKDIR]
set -euo pipefail

ALN="${1:?alignment}"
TAX="${2:?taxonomy}"
CODE="${3:?taxonomic code, e.g. BOT}"
WORK="${4:-$(mktemp -d)}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SATIVA="${SATIVA_PY:-$HERE/../sativa.py}"
PYTHON="${PYTHON:-python3}"
EPANG="${EPANG_BIN:-epa-ng}"
THREADS="${THREADS:-8}"
# identical in all three, and the only thing that has to be
ARGS=(-x "$CODE" -T "$THREADS" -m ultrafast -C 0.4 -p 42)

mkdir -p "$WORK"
ALN="$(cd "$(dirname "$ALN")" && pwd)/$(basename "$ALN")"
TAX="$(cd "$(dirname "$TAX")" && pwd)/$(basename "$TAX")"
echo "work: $WORK"

# --- A: one shot -------------------------------------------------------------------
echo "== A: one shot"
mkdir -p "$WORK/a" && cd "$WORK/a"
"$PYTHON" "$SATIVA" -s "$ALN" -t "$TAX" -n run -o . "${ARGS[@]}" > a.log 2>&1

# --- B: staged, placed in place ----------------------------------------------------
echo "== B: staged, placed in place"
mkdir -p "$WORK/b" && cd "$WORK/b"
"$PYTHON" "$SATIVA" -s "$ALN" -t "$TAX" -n run -o . "${ARGS[@]}" -stage loo-tasks > b1.log 2>&1
"$PYTHON" "$SATIVA" -stage loo-place -taskdir run.l1o_tasks -T "$THREADS"          > b2.log 2>&1
"$PYTHON" "$SATIVA" -r run.refjson -n run -o . "${ARGS[@]}" \
                    -stage loo-score -taskdir run.l1o_tasks                        > b3.log 2>&1

# --- C: staged, every fold placed somewhere else -----------------------------------
echo "== C: staged, folds placed detached"
mkdir -p "$WORK/c" && cd "$WORK/c"
"$PYTHON" "$SATIVA" -s "$ALN" -t "$TAX" -n run -o . "${ARGS[@]}" -stage loo-tasks > c1.log 2>&1

TASKS="$WORK/c/run.l1o_tasks"
STAGE="$WORK/c/staged"          # stands in for a workflow manager's work directory
mkdir -p "$STAGE"
# The command comes out of the manifest, not out of this script: if it is wrong there, the
# test fails here.
readarray -t FOLD_DIRS < <("$PYTHON" -c "
import json,sys
m=json.load(open(sys.argv[1]+'/manifest.json'))
print('\n'.join(f['dir'] for f in m['folds']))" "$TASKS")
CMD=$("$PYTHON" -c "
import json,sys,shlex
m=json.load(open(sys.argv[1]+'/manifest.json'))
print(' '.join(shlex.quote(a) for a in m['command']))" "$TASKS")
echo "   manifest command: $CMD"
echo "   ${#FOLD_DIRS[@]} folds, each placed in its own directory"

for d in "${FOLD_DIRS[@]}"; do
    # what staging does: copy the task directory in, and nothing else
    rm -rf "${STAGE:?}/$d"
    mkdir -p "$STAGE/$d"
    cp "$TASKS/$d"/* "$STAGE/$d/"
    ( cd "$STAGE/$d" && eval "${CMD/#epa-ng/$EPANG} -T 1" ) > "$STAGE/$d.log" 2>&1
    # and what collecting does: bring the placement back, nothing else
    cp "$STAGE/$d/epa_result.jplace" "$TASKS/$d/epa_result.jplace"
done

"$PYTHON" "$SATIVA" -r run.refjson -n run -o . "${ARGS[@]}" \
                    -stage loo-score -taskdir run.l1o_tasks > c3.log 2>&1

# --- D: reference, tasks, placement and scoring as four separate runs ---------------
echo "== D: four steps, four invocations"
mkdir -p "$WORK/d" && cd "$WORK/d"
"$PYTHON" "$SATIVA" -s "$ALN" -t "$TAX" -n run -o . "${ARGS[@]}" -stage reference   > d1.log 2>&1
"$PYTHON" "$SATIVA" -r run.refjson -n run -o . "${ARGS[@]}" -stage loo-tasks        > d2.log 2>&1
"$PYTHON" "$SATIVA" -stage loo-place -taskdir run.l1o_tasks -T "$THREADS"           > d3.log 2>&1
"$PYTHON" "$SATIVA" -r run.refjson -n run -o . "${ARGS[@]}" \
                    -stage loo-score -taskdir run.l1o_tasks                         > d4.log 2>&1

# --- verdict -----------------------------------------------------------------------
echo
cd "$WORK"
md5sum a/run.mis b/run.mis c/run.mis d/run.mis
fail=0
for m in b c d; do
    diff -q a/run.mis $m/run.mis > /dev/null || { echo "FAIL: $m differs from a"; diff a/run.mis $m/run.mis | head -20; fail=1; }
done
if [ $fail -eq 0 ]; then
    echo "PASS: $(wc -l < a/run.mis) mislabels, identical in all four"
    exit 0
fi
exit 1
