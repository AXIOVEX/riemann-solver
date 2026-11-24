# Source comparison — executed original vs round-4 corrected operative script

Written 2026-11-24 in the Spec 013 round-4 adjudication
(items E2-R4-1..3; adjudication:
`../coscientist-adjudication-round4.md`).

The recorded full-grid replay (`results_e2cert3.json`,
`FULL-GRID-REPLAY.md`) was executed by the original operative
script. That script is preserved byte-identical in this package
at `executed-original/e2cert3.py` (23,872 bytes, sha256
1e77032b48e6544b22633ddea293772d65aac0681ed291d040610e280da02db1).
The operative script going forward is the corrected copy at
`e2cert3.py` (25,346 bytes, sha256
5a080410456bb1f11e501cf8b0c739a9accc79ca904b79be9f8ef8a94ce428c6).

**Token-level comparison (comments stripped).** Tokenizing both
files with Python's `tokenize` and dropping COMMENT and NL
tokens gives 4,061 tokens (original) and 4,112 tokens
(corrected). The streams are identical except at exactly three
sites, all of them the E2-R4-3 metadata repair and its import:

1. The import line gains `, socket` (2 tokens).
2. The payload's venue value: the string literal
   `"tower-compute"` is replaced by the call
   `socket.gethostname()` (the actual execution host).
3. The `grid_registered` construction gains a mode conditional:
   the control-block dictionary now applies to smoke mode only,
   and selftest mode records its actual grid, built from the
   run's own grid variable (`selftest_points`, `count`, and
   `levels` keys).

The two header comment corrections (E2-R4-1, the STEP 1
widening placement; E2-R4-2, the noncontainment scope sentence)
change no tokens at all. No certification code, setting,
precision, tolerance, or numeric path differs between the two
files.

**Validation re-run.** The corrected copy was run end to end on
tower-compute (run dir
`~/.syncyard/worktrees/riemann/riemann-runs/e2cert3-r4-validation-20261124T201525Z`;
Python 3.12.3, python-flint 0.9.0). Full mode: payload exit
code 0, both families 231/231 CERTIFIED with `label_outcome`
strings identical to the recorded replay, D_k ball strings
33/33 identical to the certificate of record, W_0 identical,
both cross-check flags true, all 462 records certified, worst
discriminant relative width 1.495401392981139e-11 at
(γ, J^{12,20}) — the recorded replay's value — and the payload's
venue field reading `tower-compute`, now because that is the
host's actual name. Selftest mode: exit code 0, 7/7 certified
per family, and `grid_registered` recording exactly the seven
points run (2,0), (3,0), (3,1), (4,0), (8,10), (12,0), (12,20).
The validation payloads are run records on the venue; the
recorded payloads in this package are byte-untouched.
