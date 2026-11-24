# Schema crosswalk — certificate of record vs replay payload

Written 2026-11-24 in the Spec 013 round-3 repair pass
(adjudication: `../coscientist-adjudication-round3.md`, items
E2-R3-1 and E2-R3-2). The two JSON payloads in this package
have different schemas and different evidentiary roles. This
document is the map between them: what each file actually
contains, which comparisons between them are direct, which
assertions live on the replay side only, and which expectations
come from outside both files. Read it before applying the
README's comparison contract.

## The two payloads, as they stand

**`results_adj3.json` — the aggregate certificate-of-record
payload** (generator: adj3_r3_1.py; see the README's
Provenance). Per family it stores an aggregate outcome, not
individual records. Its complete schema:

- `venue`, `kmax` (= 32), `seconds_total`
- `r_star`: `{ball, sup}` — the r* ball string and a binary64
  rendering of its upper endpoint (see "Bound-like numeric
  fields" below)
- `D_k.<k>` for k = 0..32: `{ball, mid, sup}` — the widening
  bound's ball string and binary64 renderings
- `failure_set.<k>` for k = 0..32: `{delta: {ball, mid, lower,
  upper}, I_k_printed_rad, I_k_internal_rad,
  delta_lower_gt_diameter, delta_mid_gt_printed_rad}` — the
  round-3 noncontainment analysis of the *unwidened* balls
- `proven_outside_k` (= [0, 1, 2]), `mid_gt_rad_k` (= [0, 1, 2])
- `W_0` — the widened ball W_0, as a ball string (the only
  widened coefficient stored)
- `widened_reruns.<family>` for family in {borel, gamma}:
  `{tally: {certified, failed, indeterminate}, non_certified:
  {}}` — the aggregate widened re-run outcome. There are **no**
  per-polynomial records, no enumerated certified indices, no
  `label_outcome`, no control-block flag, and no cross-check
  flags in this file.
- `worst_disc_rel_width`: `[2.612920159759354e-11, "12,20",
  "gamma"]`

The grid the tallies refer to is **not** defined inside this
file. Grid attribution comes from the generating procedure
(adj3_r3_1.py, whose grid is `d = 2..12, n = 0..20` for both
families) and the method records (`CERTIFICATION.md`,
`EXTENSION.md`, `DK-DERIVATION.md`).

**`results_e2cert3.json` — the recorded replay payload**
(generator: the executed-original script — the operative
script's original bytes, sha256 1e77032b…, now packaged at
`executed-original/e2cert3.py` — full mode). This is where
the individual polynomial records live:

- `results_by_family.<family>.tallies`:
  `{certified, failed, indeterminate, total}` (note: this
  schema differs from the certificate's tally — it adds
  `total`, and it sits at a different path)
- `results_by_family.<family>.polynomials`: a dict keyed
  `"d,n"` holding one record per polynomial (231 per family),
  each with `status` ("certified" / "failed" / "indeterminate")
  and its Sturm/discriminant evidence
- `results_by_family.<family>.control_block_77` (the 77-point
  control-block tally), `.control_block_all_certified`
  (boolean), `.label_outcome` (the family's verdict string)
- `levels[].cross_check_vs_E2_all_intersect` and
  `levels[].cross_check_vs_E1_all_intersect`: the moment
  cross-check flags (booleans, with per-k detail alongside)
- `widening.D_k.<k>.ball`: the widening bounds as ball strings
  (string-identical to the certificate's, 33/33)
- `moments_widened_final_level."0".ball`: the W_0 ball string
- `grid_registered`: `{d: "2..12", n: "0..20", count: 231}`

## Direct comparisons (fields present on both sides)

| Comparison | Certificate path | Replay path |
|---|---|---|
| Per-family tally | `widened_reruns.<family>.tally` {certified, failed, indeterminate} | `results_by_family.<family>.tallies` (same three counts; replay adds `total`) |
| Non-certified set | `widened_reruns.<family>.non_certified` (empty) | the set of `"d,n"` keys in `results_by_family.<family>.polynomials` whose `status` is not "certified" (empty) |
| Widening vector | `D_k.<k>.ball`, k = 0..32 | `widening.D_k.<k>.ball` |
| W_0 | `W_0` | `moments_widened_final_level."0".ball` |
| Recorded internal radius floats | `failure_set.<k>.I_k_internal_rad` | `moments_unwidened_final_level.<k>.rad` (the recorded float representation of the recomputed internal radius) |
| Recorded printed radius floats | `failure_set.<k>.I_k_printed_rad` | parse the radius from `moments_unwidened_final_level.<k>.ball` (the replay's printed ball string), then apply the diagnostic conversion — convert that printed decimal radius to a float — if comparing floats |
| Worst width triple | `worst_disc_rel_width` | the maximum `rel_width` over `results_by_family.<family>.polynomials.<key>.discriminant`, with its key and family (a regression check, per the README's contract — the values legitimately differ) |

Agreement on the tally and the empty non-certified set is
agreement with the certificate's **aggregate** all-certified
result. It is not, by itself, an inspection of individual
records; the individual records exist only on the replay side.

The two radius rows are separate comparisons with different
replay-side extractions (the internal radius is a stored float
field in the replay payload; the printed radius is parsed from
the replay's printed ball string). Both were verified in the
round-3 adjudication and re-verified in the round-4
adjudication: 33/33 agree in each case.

## Replay-only assertions (no counterpart in the certificate)

These are properties of the replay payload itself. A reader
verifies them against the replay payload and the externally
specified expectations below — never against
`results_adj3.json`, which does not contain the fields:

- **Coverage check (required):** for each family, the keys of
  `results_by_family.<family>.polynomials` are exactly the 231
  pairs `"d,n"` with d = 2..12 and n = 0..20, each exactly
  once. (Verified for the recorded replay in the round-3
  adjudication: exact, both families, all 462 records with
  `status` "certified".)
- Every individual record's `status`.
- `control_block_all_certified` true, and the control-block
  tally `control_block_77` = 77/0/0, for each family.
- `label_outcome` equal to the expected verdict string quoted
  in the README's "Commands and expected outcomes".
- Both cross-check flags (`cross_check_vs_E2_all_intersect`,
  `cross_check_vs_E1_all_intersect`) true at the deciding
  level.

## Externally specified expectations (in neither payload)

- The family definitions (standard ξ vs Borel) and the grid
  (d = 2..12, n = 0..20; control block d = 2..8, n = 0..10) —
  stated in the README's claim section.
- The expected `label_outcome` strings — quoted in the README.
- The grid attribution of the certificate's tallies — from the
  generating procedure and method records, as above.

## Bound-like numeric fields: diagnostics, not bounds

The certificate's numeric `mid`, `sup`, `lower`, and `upper`
fields (in `r_star`, `D_k.<k>`, and `failure_set.<k>.delta`)
are binary64 renderings of Arb ball midpoints and endpoints —
`float(ball.mid())`, `float(ball.upper())`, and so on, written
by the generating script for convenience. **They are not
outward bounds.** A binary64 rendering of an endpoint can fall
on either side of the true endpoint: for example, the exact
binary64 value of `D_k["6"].sup` lies *below* the lower
endpoint of its own ball, and the exact binary64 value of
`failure_set["0"].delta.lower` lies *above* the upper endpoint
of its ball (both reproduced exactly in the round-3
adjudication). Do not use these fields as directed bounds.

- **Authoritative representation:** the `ball` string fields.
  A ball string `"[mid +/- rad]"` denotes the enclosure
  [mid − rad, mid + rad] in exact decimal arithmetic on the
  printed strings.
- **Safe endpoint extraction:** parse the ball string and
  compute mid − rad and mid + rad in exact decimal (or
  rational) arithmetic. `CERTIFICATE-ENDPOINTS.json` in this
  package does exactly that for `r_star`, every `D_k`, every
  `failure_set` delta, and `W_0` (all 68 balls in the
  certificate), and is hash-pinned in `INVENTORY.md`. The
  certificate itself is byte-untouched.
- **Script usage (verified in the round-3 adjudication):** in
  `e2cert3.py`, every `float(...)` value is a written
  diagnostic; all certification decisions are Arb ball
  comparisons. In adj3_r3_1.py, the widened re-runs likewise
  decide entirely in Arb (the widened balls are built from the
  Arb endpoint `Dk[k].upper()`, never from the float). The one
  float-mediated test in adj3_r3_1.py is the Part 1
  noncontainment flag `delta_lower_gt_diameter`, computed as
  `float(delta.lower()) > 2 * I_k_internal_rad * (1 + 1e-9)`
  (the diameter deliberately inflated for float rounding). It
  is margin-dominated at every index — the smallest ratio on
  the true side is 5.58 (k = 2), the largest on the false side
  is 4.0e-4 (k = 3), against a binary64 relative error of
  about 1.1e-16 — and it scopes only the *unwidened* balls'
  noncontainment (`proven_outside_k`), not the widened
  certification.
