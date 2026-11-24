# Proof Map — Certified Log-Concavity of the Riemann–Jacobi Kernel Φ

One page. Every claim in the paper, what proves it, and what a human can
check without a computer. Φ(u) = Σ_{n≥1} (2π²n⁴e^{9u/2} − 3πn²e^{5u/2})
exp(−πn²e^{2u}); φ₁ is its n = 1 term; Q(f) := f f″ − (f′)², so
(log f)″ < 0 is Q(f) < 0.

## The main theorem: (log Φ)″(u) < 0 for all u ≥ 0

This is Csordas's theorem (2015, Thm 4.5), resting on Csordas–Varga
(1988) Thm 2.1 via Thm 4.2(b). Our contribution is an independent
interval-arithmetic certification, in three ranges:

| Range | What proves it | Human-checkable unaided |
|---|---|---|
| Algebraic core (all u ≥ 0) | (log φ₁)″ < 0 by one explicit computation: the numerator is a manifestly negative combination. Hand argument (Lemma, §3). | Yes — differentiate and inspect; a page of algebra. |
| [0, 1] | Arb ball-arithmetic certificate, 55,892 intervals (`proof/verify_logconcavity_arb.py`), cross-checked by an independent mpmath.iv run, 52,898 subintervals (`verify_logconcavity_rigorous.py`). The n > 5 truncation is bounded by a uniform hand estimate: total perturbation ≤ 5.34×10⁻³⁶ against a certified margin 3.36×10⁻¹² (`verify_truncation_and_crosscheck.py`). | The truncation estimate: yes (majorants + one propagation line). The interval runs: no — machine certificates; the scripts are 225 + 207 + 224 lines. |
| [1, 3] | Direct ΔQ enclosure in Arb ball arithmetic over 1,813 adaptive balls (`proof/verify_ia_1_to_3_arb.py`): at each ball, φ₁ and the finite tail n = 2..10 are enclosed from their closed forms, ΔQ is enclosed from its definition, the n > 10 tail is bounded by explicit majorants (≤ 6.56×10⁻¹¹⁷⁴), and the sign test W₁ + B < 0 is an Arb comparison. Cross-checked by mpmath.iv at 101 checkpoints (`verify_ia_1_to_3.py`). Max B ≤ 8.54×10⁻²⁴, min margin 91.29. | The ΔQ formula and the tail majorants: yes. The 1,813 ball evaluations: no — machine certificate; the script is 421 lines. |
| [3, ∞) | Hand argument: weighted comparisons S₆ ≤ 5S₄, S₈ ≤ 32S₄ (n = 2 domination; the two majorant sums are certified once, at u = 3, by `verify_weighted_comparisons.py`, 164 lines), then termwise monotonicity gives |W_tail| ≤ 2600·e^{4u}ε*(u), and W₁(u) ≤ −2500 dominates. C_tail = 2600. | Yes, except the two certified sums (each < 10⁻²⁷⁴⁹, enormous slack). |

## The counterexample (the criterion is false)

K(t) = (5 + t⁴)·exp(−t² − t⁴/1000) satisfies all five hypotheses of the
Hypothetical Criterion, and its Fourier transform has a non-real zero.

| Step | What proves it | Human-checkable |
|---|---|---|
| K is strictly log-concave | Displayed identity for (log K)″: numerator −5 − 20(t²−3/2)² − 4t⁶ − 2t⁸ < 0. | Yes — one differentiation. |
| The other hypotheses | Decay e^{−t²/1000}-type bound, evenness, entireness of the transform: direct estimates. | Yes. |
| A non-real zero exists | Hadamard factorization (order ≤ 4/3) forces the necessary condition m₀m₄ ≤ 3m₂² on the moments; exact Gaussian-moment bounds give m₀m₄ − 3m₂² ≥ (26544867/40960000)·π > 0. Contradiction. | Yes — the moments are elementary Gaussian integrals; the arithmetic is exact rationals. |
| Corroboration | A zero at 2.8086 + 1.3359i, winding-number certified (numerical, in the program's adjudication record). | No — corroboration only; the proof above stands alone. |

Precedent: Csordas (2015), Example 3.12, is the same phenomenon.

## The source audit (hand scholarship, no computation)

- No audited primary source states the criterion; Pólya's Satz I/II are
  multiplier theorems, not the criterion. Hilfssatz II is Cardon (2005),
  Lemma 2, after Pólya (1926) — a shift-sum preservation lemma.
- Csordas–Varga (1988) Thm 2.1 / Csordas (2015) Thm 4.2(b) prove the
  stronger statement that log Φ(√t) is concave, under a Mellin condition
  not established for Φ here.
- Consequence: the conditional theorem "criterion ⇒ RH" is vacuous.
  This paper claims nothing about RH.

## What a human can check in one sitting

The algebraic core, the whole counterexample, the far-tail argument, and
every lemma statement. What requires trusting a machine: the two
interval-arithmetic ranges, [0,1] and [1,3] — that is what the
certificates are, and the certificate path is 8 scripts, 2,001 lines,
run end-to-end by `verify.py` (fail-closed; exit 0 only if every step
passes).
