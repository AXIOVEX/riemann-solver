# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Tristen Kyle Pierson / BitConcepts Research
"""Run the complete proof verification pipeline.

Steps (all fail-closed; a step counts as passed only on exit code 0):
  1. Arb/FLINT primary IA: Q_Phi < 0 on [0, 1.0] (55,892 subintervals)
  2. Rigorous mpmath.iv IA: Q_Phi < 0 on [0, 1.0] (52,898 subintervals)
  3. Truncation error certification (uniform, [0,1]) + cross-validation
  4. Arb/FLINT primary cert: (log Phi)'' < 0 on [1.0, 3.0]
     (direct DeltaQ enclosure, adaptive ball grid)
  5. mpmath.iv cross-check: (log Phi)'' < 0 on [1.0, 3.0] (101 checkpoints)
  6. Weighted tail comparisons for u >= 3 (far-tail lemma constants)
  7. Galerkin stabilization check

The final verdict additionally requires cross_validation_passed = true
in the truncation_crosscheck.json written by THIS run.  A step that
reports SKIP (missing optional dependency) counts as FAILED: in the
certified pipeline a step that verified nothing is not a pass.

Retired from the pipeline in the Spec 004 repair pass:
verify_algebraic_core.py (legacy; printed a withdrawn proof claim),
verify_ia_1_to_1_5.py (legacy C=204 perturbation script, formerly
mislabeled as the [1,3] cross-check), verify_debruijn_condition.py
(exploratory; not part of the proof).  All three are banner-marked
SUPERSEDED in place.

Usage:
    python verify.py           # Full pipeline
    python verify.py --quick   # Skip the two [0,1] IA steps (partial verdict)
"""
import argparse
import subprocess
import sys
import os
import json
import time

PROOF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "proof")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

STEPS = [
    ("verify_logconcavity_arb.py",       "Arb/FLINT primary IA: 55,892 subintervals on [0, 1.0]"),
    ("verify_logconcavity_rigorous.py",  "Rigorous IA: 52,898 subintervals on [0, 1.0]"),
    ("verify_truncation_and_crosscheck.py", "Truncation error (uniform) + cross-validation"),
    ("verify_ia_1_to_3_arb.py",          "Arb/FLINT primary cert: (log Phi)'' < 0 on [1.0, 3.0] (direct DeltaQ enclosure)"),
    ("verify_ia_1_to_3.py",              "mpmath.iv cross-check: (log Phi)'' < 0 on [1.0, 3.0] (101 checkpoints)"),
    ("verify_weighted_comparisons.py",   "Weighted tail comparisons: S_6 <= 5*S_4, S_8 <= 32*S_4 for u >= 3"),
    ("verify_galerkin_extended.py",      "Galerkin stabilization: CvS form plateau, eigenvector c-invariance (--quick)"),
]

QUICK_SKIP = {"verify_logconcavity_arb.py", "verify_logconcavity_rigorous.py"}

# Scripts that require extra flags when run via the pipeline
# CI-specific: limit Galerkin to 3 cutoffs (c=7,11,13) for speed on the test matrix;
# the full 9-cutoff quick run is preserved for local use via `python proof/verify_galerkin_extended.py --quick`
_GALERKIN_FLAGS = ["--quick", "--cutoff-max", "13"] if os.getenv("CI") else ["--quick"]
_SCRIPT_FLAGS = {
    "verify_galerkin_extended.py": _GALERKIN_FLAGS,
}


def run_step(script, desc, quiet=False):
    path = os.path.join(PROOF_DIR, script)
    if not os.path.exists(path):
        print("  FAIL: %s not found" % script)
        return False

    extra_flags = _SCRIPT_FLAGS.get(script, [])
    # Propagate UTF-8 encoding to child processes (needed on Windows CI)
    child_env = {**os.environ, "PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8"}
    result = subprocess.run(
        [sys.executable, path] + extra_flags,
        capture_output=True, text=True, timeout=1200,
        cwd=os.path.dirname(PROOF_DIR),
        env=child_env,
    )

    if result.returncode != 0:
        print("  FAIL (exit code %d)" % result.returncode)
        if result.stderr:
            print(result.stderr[-500:])
        return False

    # A SKIP report (optional dependency missing) is a FAILED step in the
    # certified pipeline: the step verified nothing.
    if result.stdout and result.stdout.strip().startswith("SKIP:"):
        print("  FAIL (step skipped: %s)" % result.stdout.strip().split("\n")[0])
        return False

    if not quiet:
        # Print last few meaningful lines
        lines = [l for l in result.stdout.strip().split("\n") if l.strip()]
        for line in lines[-6:]:
            print("  " + line)

    return True


def main():
    parser = argparse.ArgumentParser(description="Run RH proof verification pipeline")
    parser.add_argument("--quick", action="store_true",
                        help="Skip the two [0,1] IA steps (partial verdict)")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress step output")
    args = parser.parse_args()

    steps = [(s, d) for (s, d) in STEPS
             if not (args.quick and s in QUICK_SKIP)]

    print("=" * 72)
    print("  RH PROOF VERIFICATION PIPELINE")
    print("  %d steps%s" % (len(steps), " (quick mode — [0,1] IA steps skipped)" if args.quick else ""))
    print("=" * 72)

    run_start = time.time()
    passed = 0
    for i, (script, desc) in enumerate(steps, 1):
        print("\n[%d/%d] %s" % (i, len(steps), desc))
        if run_step(script, desc, quiet=args.quiet):
            passed += 1

    # Cross-validation enters the final verdict, and only from the file
    # written by this run (a stale file from an earlier run does not count).
    xv_path = os.path.join(RESULTS_DIR, "truncation_crosscheck.json")
    xv_ok = False
    if os.path.exists(xv_path) and os.path.getmtime(xv_path) >= run_start - 1:
        with open(xv_path) as f:
            xv = json.load(f)
        xv_ok = bool(xv.get("cross_validation_passed", False)
                     and xv.get("passed", False))

    print("\n" + "=" * 72)
    print("  RESULTS: %d/%d steps passed" % (passed, len(steps)))
    print("  Cross-validation: %s" % ("PASSED" if xv_ok else "NOT RUN or FAILED"))

    if not args.quick:
        ia_path = os.path.join(RESULTS_DIR, "verify_logconcavity_rigorous.json")
        if os.path.exists(ia_path):
            with open(ia_path) as f:
                ia = json.load(f)
            print("  IA certified: %d/%d subintervals" % (ia["certified"], ia["n_subintervals"]))

    print("=" * 72)
    ok = (passed == len(steps)) and xv_ok
    if args.quick:
        print("  NOTE: quick mode — verdict covers only the steps run.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
