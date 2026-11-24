# E2 Certification Extension — Jensen grid extended to d=2..12, n=0..20, all CERTIFIED (2026-11-22)

**Correction (2026-11-23, Co-scientist round-2 adjudication R2-1, SUSTAINED; Spec 005 R5-01).** The family certified in this run is misidentified in the original text. Definitions: H(z) := ∫₀^∞ Φ(u)cosh(u√z)du = Σ_{k≥0} a_k z^k = Σ_{k≥0} γ_k z^k/k!, so γ_k = k!·a_k. The standard Jensen polynomials of Ξ are J_H^{d,n}(x) = Σ C(d,k) γ_{n+k} x^k. The polynomials certified here, Σ C(d,k) a_{n+k} x^k, are the Jensen polynomials of B(z) = Σ a_k z^k/k!, the Borel transform of H. The original Reading below said: "The certified Jensen range for Ξ now covers d ≤ 12, n ≤ 20: for every one of the 231 grid points, J^{d,n} is proved … to have d distinct real roots." Corrected label: **CERTIFIED (Borel family of the moment sequence)**. The root counts stand verbatim as theorems about the Borel family. The standard grid's CERTIFIED label rests on the round-2 adjudication's consequence run, reported in the addendum at the end of this file. The falsifier direction for the Borel family holds via the 1/k! multiplier argument (Pólya–Schur; Σ z^k/(k!)² = I₀(2√z) ∈ LP), not via the identification originally printed. The identification phrases below are corrected in place under this note, and the definitions are also inserted after the moment definition in the Method section.

**Venue:** tower-compute, run dir ~/riemann-runs/e2cert2-20261122T154437, TC-03 venv, python-flint 0.9.0 (FLINT ball arithmetic). **Payload:** results_e2cert2.json (305,347 bytes) read back from the run directory; every number below is quoted from that file, not from the launch wrapper (Spec 002 TC-03). Script: e2cert2.py, preserved in the run dir and alongside this report (md5 b1a6a914908bfbfcc5a7c7a5aecc7af1). *Documentation note (2026-11-23):* the alongside copies of e2cert.py and e2cert2.py now carry corrected header documentation (R2-1/R2-1b — identification, truncation scoping, summation tail; comments only, no executable line changed), so their hashes no longer match the run-dir originals; the md5 quoted here is the run-dir original's, which is unchanged and remains the artifact of record for this run. Total payload wall time **0.77 s**; **all 231 polynomials certified at the first precision level (256 bits), with zero retries, zero indeterminate signs, and zero exceptions** — the escalation ladder (512/1024/2048 bits) was never entered (231 attempts logged, one per polynomial).

**Classification (Spec 002 R1): CERTIFIED (Borel family of the moment sequence)** — for this grid only. Hyperbolicity on any finite grid, certified or not, proves nothing beyond the grid (spec §1).

## R5 label outcome, stated exactly

> **CERTIFIED (Borel family of the moment sequence) for the extended grid d=2..12, n=0..20 (231/231) — the Spec 002 R5 label is earned for this grid only.**

This is the expected outcome under RH, it extends the rigorously verified range, and it says nothing about (d, n) off the grid.

## Registration (fixed in the script header before running; not extended after seeing results)

- Grid: d = 2..12, n = 0..20, J^{d,n}(x) = Σ_{k=0}^{d} C(d,k) a_{n+k} x^k — 11 × 21 = **231 points**.
- **Control block:** the 77 points (d = 2..8, n = 0..10) CERTIFIED by run e2cert-20261122T145019, re-run here under this run's independently computed interval moments; the batch is NO VERDICT unless all 77 re-certify.
- New territory: the remaining 154 points. Interval moments a_0..a_32.

## Tallies (from the payload)

| Block | CERTIFIED | FAILED | INDETERMINATE | Total |
|---|---|---|---|---|
| Control block (d=2..8, n=0..10) | 77 | 0 | 0 | 77 |
| New territory | 154 | 0 | 0 | 154 |
| **Full extended grid** | **231** | **0** | **0** | **231** |

**Control-block outcome: all 77 re-certified** under this run's independent moments, so the batch counts (not a NO VERDICT). The four TC-02 controls J^{2,0}, J^{3,0}, J^{3,1}, J^{4,0} are inside the control block and all CERTIFIED. There are no INDETERMINATE polynomials, so there are no obstructions to name, and no FAILED (falsifier) candidates.

## Method (identical to the cert run, restated where the extension changes it)

**Interval moments.** Each a_k = ∫₀^∞ Φ(u) u^{2k} du / (2k)!, k = 0..32, enclosed as (∫₀⁸ + tail ball) / (2k)!:

*Definitions (inserted 2026-11-23, R2-1).* H(z) := ∫₀^∞ Φ(u)cosh(u√z)du = Σ_{k≥0} a_k z^k = Σ_{k≥0} γ_k z^k/k!, with γ_k = k!·a_k. The standard Jensen polynomials of Ξ are J_H^{d,n}(x) = Σ C(d,k) γ_{n+k} x^k. The polynomials certified in this run, Σ C(d,k) a_{n+k} x^k, are the Jensen polynomials of B(z) = Σ a_k z^k/k!, the Borel transform of H.

- ∫₀⁸ by FLINT's `acb_integral` in ball arithmetic (guaranteed enclosure; the integrand as evaluated is the finite six-term truncation sum (n = 1..6), entire in u — the full Φ series is analytic only in the strip |Im u| < π/4; scoping corrected 2026-11-23, R2-1b; true integral real and in the ball's real part). Requested rel_tol = 2⁻¹⁹² at 256-bit precision.
- **Summation tail (n ≥ 7), documented 2026-11-23 (R2-1b).** As in the cert run: the omitted n ≥ 7 terms contribute at most 5.28×10⁻⁶² relative to each moment (at k = 0, decreasing in k), about 3.5 orders below the enclosures' maximum relative width; see CERTIFICATION.md for the full statement and its source.
- **Tail bound (restated for k ≤ 32).** Φ(u) ≤ 4550π² e^{9u/2 − πe^{2u}} for u ≥ 0 (same per-term majorant as the cert run). For u ≥ 8 and k ≤ 32, u^{2k} ≤ u^{64}; with g(u) = (9/2)u + 64 ln u − πe^{2u}, g′(u) ≤ 9/2 + 8 − 2πe^{16} =: −λ < 0, giving ∫₈^∞ Φ(u)u^{2k}du ≤ 4550π² e^{g(8)}/λ, g(8) = 36 + 64 ln 8 − πe^{16}. Evaluated in ball arithmetic the bound is **≈ 10^{−12,123,929}** (payload log10 ≈ −12123928.73) — zero at any working precision, but bounded and added as a ball containing [0, bound] to every enclosure before the per-k division by (2k)!.

**Interval widths and cross-checks.** Max relative width over a_0..a_32 at 256 bits: **1.90×10⁻⁵⁸** (the same floor as the cert run). Cross-checks: every enclosure intersects the corresponding stored numerical value — against the E2 payload (40 significant digits) for k = 0..18 **and** against the E1 payload moments (80 dps, moments_partial.json) for all k = 0..32. All intersect, at the single level computed.

**Sturm certification with ball coefficients.** Identical straight-line Sturm construction from the coefficient balls C(d,k)·a_{n+k}: every divisor and Sturm leading-coefficient ball excluded 0 on the first attempt for all 231 polynomials, so each computed ball sequence encloses the exact Sturm sequence of every member polynomial; sign variations at ±∞ give the exact distinct-real-root count, and the final Sturm constant ball excluded 0 in all 231 cases (separability proved). Every Sturm sequence had the full length d+1 (lengths 3..13 across the grid) and root count exactly d. Secondary check: Sylvester-determinant discriminant enclosures **positive for all 231**, as distinct real roots require. No internal inconsistencies.

## Discriminant margins — the carry-forward dataset

The payload records the full discriminant ball for every polynomial (mid, radius, relative width, sign margin (|mid| − rad)/|mid|, log10|mid|). At float resolution every sign margin reads 1.000000 — no discriminant enclosure comes anywhere near straddling 0 — so the informative margin variable is the enclosure **relative width**, which grows toward the high-(d, n) corner. The tightest enclosures on the grid:

| Polynomial | disc rel. width | log10\|disc\| |
|---|---|---|
| J^{12,20} | 2.12×10⁻¹⁸ | −1677.17 |
| J^{12,19} | 1.37×10⁻¹⁸ | −1602.06 |
| J^{12,18} | 1.65×10⁻¹⁹ | −1527.44 |
| J^{12,13} | 1.42×10⁻¹⁹ | −1162.48 |
| J^{12,12} | 1.21×10⁻¹⁹ | −1091.28 |

Even the worst case (J^{12,20}) pins the discriminant's sign with ~18 digits of headroom at 256 bits, which is why no precision escalation was needed anywhere. The loosest scale on the grid is at the opposite corner (J^{2,0}: rel. width 5.9×10⁻⁷¹, log10|disc| = −4.75). If a future extension pushes (d, n) higher, this relative-width curve — worst at maximal d and n — is the quantity that decides when the 512-bit rung of the ladder first gets used; the per-polynomial balls in the payload are the dataset to extrapolate from.

## Reading

The certified Jensen range for B, the Borel transform of H, now covers d ≤ 12, n ≤ 20: for every one of the 231 grid points, J^{d,n} is proved — against interval enclosures of the coefficients, so the conclusion covers the true polynomials — to have d distinct real roots. The 77-point predecessor grid re-certified in full under independent moments, so the extension rests on a reproduced base, not on the earlier run's say-so. This remains a finite grid: it is the expected outcome under RH, it extends the rigorously verified range, and it says nothing about (d, n) off the grid. A single FAILED verdict anywhere on any future grid would be a falsifier (one non-hyperbolic Jensen polynomial refutes RH itself — for this Borel family via the 1/k! multiplier argument in the correction note; for the standard family via the Griffin–Ono–Rolen–Zagier equivalence); none appeared.

## Addendum (2026-11-23) — the standard ξ Jensen grid, certified by the γ run

The round-2 adjudication's consequence run re-ran the polynomial-level certification with the standard coefficients γ_{n+k} = (n+k)!·a_{n+k} from the same interval-moment machinery (no new integrations; only the coefficient construction changed). Script: e2gamma.py (sha256 09850a704a0447e159b775c0567b7cc24145ddc092939c1e48b0868f558aed61), e2cert2.py with the single coefficient line changed to `comb(d,k) * factorial(n+k) * a[n+k]`. Venue: tower-compute, run dir ~/riemann-runs/e2gamma-adj2-20261123T114853, TC-03 venv; payload exit code 0 (mechanically captured). Payload: results_e2gamma.json (sha256 438d4e7f4b0f151aa99ba803099421e6be860c5e7ed98afff8a3906f0fa55acc), copy alongside this report.

**Result: the standard ξ Jensen grid is CERTIFIED hyperbolic 231/231** (d = 2..12, n = 0..20) at 256 bits — control block 77/77, new territory 154/154, zero failed, zero indeterminate. Moment enclosures match the runs above in character (maximum relative width 1.90×10⁻⁵⁸; cross-checks against the E1 and E2 stored moments all intersect). The standard family sits closer to the boundary than the Borel family: the worst discriminant relative width is 1.50×10⁻¹¹ at J^{12,20} (against 2.12×10⁻¹⁸ for the Borel family), still ~11 digits of headroom at the first precision level. The standard grid and the Borel grid agree on the full range; any standard-grid CERTIFIED claim cites this payload, and the 2026-11-22 runs are cited for the Borel family only.

## Round-3 completion note (2026-11-23, Co-scientist round-3 adjudication R3-1, SUSTAINED; Spec 006) — containment failure and widened-ball repair

The summation-tail bullet in the Method section above states that the omitted n ≥ 7 terms contribute "at most 5.28×10⁻⁶² relative to each moment (at k = 0, decreasing in k), about 3.5 orders below the enclosures' maximum relative width", and the γ-run addendum above rests the standard grid's CERTIFIED label on the same moment balls ("Moment enclosures match the runs above in character (maximum relative width 1.90×10⁻⁵⁸…)"). **The containment conclusion those statements support is false as stated, and both are superseded by this note.** Comparing the omitted tail against the *maximum* relative width across all 33 moments is invalid per ball: a_0's own relative width is 4.36e-75. The round-3 adjudication enclosed the omitted contributions δa_k directly in Arb at 512 bits (tower run dir ~/riemann-runs/adj3-20261123T0725; position-independent criterion: δa_k's lower bound vs the accepted ball's full diameter):

- **Containment failure set: exactly k = 0, 1, 2** — δa_0 = 2.148808e-65, δa_1 = 2.259718e-70, δa_2 = 2.316534e-75, each exceeding its accepted ball's diameter, so the true moments lie outside the balls I_k recorded in `results_e2cert2.json` and `results_e2gamma.json` (the moments are identical across all three original payloads). For k ≥ 3, δa_k is at least three orders of magnitude below the ball radius; no failure is demonstrable there.
- **Repair (executed in the round-3 adjudication):** widened balls W_k = I_k + [0, D_k] with the §B bounds D_k (D_0 = 2.184234e-65; independently re-derived, and δa_k ≤ D_k confirmed for every k), then polynomial-level Sturm re-runs only, no new integrations. Verdicts: **Borel 231/231 CERTIFIED and standard γ 231/231 CERTIFIED** (0 failed, 0 indeterminate), worst discriminant relative width 2.61e-11 (γ family at J^{12,20}; the unwidened γ run's was 1.50×10⁻¹¹ — the honest balls cost almost nothing).
- **Payload supersession:** `results_e2cert.json`, `results_e2cert2.json`, and `results_e2gamma.json` are **SUPERSEDED as true-kernel certificates**; they stand as records of the originally enclosed coefficient sets. The operative certificate of record for both families is now `results_adj3.json` (sha256 981078c5f7f25177b9791f717e97af4772d7babbc770ce866aa5cd7a1d399769), copy alongside this report. Any future standard-grid or Borel-grid CERTIFIED claim cites `results_adj3.json`, not the payloads above.
- **The root counts themselves never changed**: the recorded Sturm counts were correct for the coefficients actually enclosed, and the widened re-runs return identical counts for the true kernel's coefficients. Any future E2 run must include the [0, D_k] widening (or an equivalent n ≥ 7 tail ball) in its moment enclosure.

R2-1b is accordingly DONE (repair complete) as of this note; see the ledger and `specs/006-round3-repair-pass/REPORT.md`.


## R4-1 correction note (2026-11-23, round-4 adjudication, SUSTAINED) — scope of the failure-set statement

The round-3 completion note above calls the containment failure set "exactly k = 0, 1, 2". That phrasing overstates the test, as the round-4 audit showed and the R4-1 adjudication (appended to `coscientist-adjudication-round3.md`, 2026-11-23) sustained: the position-independent criterion (δa_k's lower bound exceeding the accepted ball's diameter) proves the true moments lie outside the old balls for k = 0, 1, 2, but for k ≥ 3 — where δa_k is at least three orders of magnitude below the ball radius — the test establishes neither containment nor noncontainment in the old balls, since the truncated moment could sit anywhere in its ball (near an endpoint, the small omitted contribution would still push it outside). The corrected sentence, adopted verbatim from the audit, is:

> Position-independent noncontainment is proved for k=0,1,2. This test does not establish containment in the old balls for the remaining moments. After widening, containment is established for all k=0,…,32.

This corrects the scope of the retrospective diagnosis only. The widening repair (W_k = I_k + [0, D_k] for every k; `results_adj3.json` as the certificate of record, both families 231/231) is unaffected.
