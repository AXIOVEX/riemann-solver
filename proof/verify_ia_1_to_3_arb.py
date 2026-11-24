# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Tristen Kyle Pierson / BitConcepts Research
"""Arb/FLINT ball-arithmetic certificate: (log Phi)''(u) < 0 on [1.0, 3.0].

This is the PRIMARY rigorous interval certificate for [1,3].
The mpmath.iv verifier (verify_ia_1_to_3.py) is retained as an independent cross-check.

Method (Appendix E.5, as executed): certify W_1(u) + B(u) < 0 on adaptive
overlapping ball-arithmetic intervals covering [0.99, 3.01] ⊃ [1, 3]
(half-width min(0.01, 0.35/(2*pi*e^{2u})), so the exponent -pi*e^{2u} of
the dominant term varies by at most ~0.35 across any one ball), where
  W_1(u) = (log phi_1)''(u)  (closed form, both terms negative),
  B(u)   = |DeltaQ(u)| / phi_1(u)^2 + 3 * eps*(u) * |W_1(u)|
is a rigorous upper bound on |W_tail(u)| (Lemma Wtail_bound), and DeltaQ is
enclosed DIRECTLY from its definition
  DeltaQ = R''*phi_1 + phi_1''*R - 2*phi_1'*R' + R''*R - (R')^2,
with R = sum_{n>=2} phi_n split into the finite part n = 2..10 (evaluated in
ball arithmetic from the exact closed forms) and the tail n > 10, whose
contribution is bounded by explicit derivative-aware majorants
(tail_allowances_arb).  eps*(u) is likewise a finite ball evaluation
(n = 2..15) plus an explicit majorant tail bound (n > 15).

All proof-critical comparisons are performed in Arb: the sign test is the
ball comparison (W_1 + B) < 0, never a float endpoint sum.  Floats are used
only to format display values and JSON summaries.

History: the previous version of this script certified a surrogate bound
(2500*e^{4u}*eps* + 15*pi*e^{2u}*eps*) imported from a lemma stated only
for u >= 3, instead of enclosing DeltaQ.  The bound was numerically true
with wide margin but was not the method Appendix E.5 describes.  This
version implements the direct enclosure.  See Spec 004 / the round-1
adjudication, Finding C-3.3.
"""

import json
import time
import sys
import os

try:
    from flint import arb, ctx
except ImportError:
    print("SKIP: python-flint not installed. Install with: pip install python-flint")
    print("The Arb/FLINT [1,3] certificate requires python-flint.")
    sys.exit(0)

PREC = 256          # bits of precision
HW_MAX = 0.01       # maximum checkpoint half-width
EXP_HW = 0.35       # max half-width of the exponent -pi*e^{2u} per ball
N_FINITE = 10       # finite terms in the DeltaQ sums (n = 1..10)
N_EPS_TERMS = 15    # finite terms for epsilon* (n = 2..15), tail bounded
N_TAIL_SUM = 80     # majorant sums for tails run n = 11..80 (resp. 16..80)
U_LO = arb("0.99")  # global cover endpoints, used by the tail majorants
U_HI = arb("3.01")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# Arb helper functions
# ---------------------------------------------------------------------------

def upper(x):
    """Rigorous upper endpoint of an arb ball, as a float (display only)."""
    return float(x.mid()) + float(x.rad())


def fmt(x):
    """Short decimal string for an arb's upper endpoint (safe for values
    that underflow float)."""
    return x.upper().str(6)


def make_ball(center, radius):
    """Create an arb ball [center ± radius]."""
    return arb(str(center), str(radius))


# ---------------------------------------------------------------------------
# Mathematical functions using Arb ball arithmetic
# ---------------------------------------------------------------------------

def log_phi1_d2_arb(u_ball):
    """Compute W_1(u) = (log phi_1)''(u) = -24*pi*e^{2u}/h(u)^2 - 4*pi*e^{2u}.

    h(u) = 2*pi*e^{2u} - 3.
    """
    pi = arb.pi()
    two = arb(2)
    e2u = (two * u_ball).exp()
    h = two * pi * e2u - arb(3)
    return -arb(24) * pi * e2u / (h * h) - arb(4) * pi * e2u


def phi_n_all_arb(u_ball, n):
    """Compute (phi_n, phi_n', phi_n'') using Arb ball arithmetic.

    Exact closed-form expressions from Appendix E:
      g_n(u) = 2*pi^2*n^4*e^{9u/2} - 3*pi*n^2*e^{5u/2}
      E_n(u) = exp(-pi*n^2*e^{2u})
      phi_n = g_n * E_n, etc.
    """
    pi = arb.pi()
    n_arb = arb(n)
    n2 = n_arb * n_arb
    n4 = n2 * n2

    # Shared exponentials
    two = arb(2)
    e5h = (arb(5) / two * u_ball).exp()     # e^{5u/2}
    e9h = (arb(9) / two * u_ball).exp()     # e^{9u/2}
    e2u = (two * u_ball).exp()               # e^{2u}
    e4u = e2u * e2u                          # e^{4u}

    pi2 = pi * pi

    # g_n and derivatives
    g   = two * pi2 * n4 * e9h - arb(3) * pi * n2 * e5h
    gd1 = arb(9) * pi2 * n4 * e9h - arb(15) / two * pi * n2 * e5h
    gd2 = arb(81) / two * pi2 * n4 * e9h - arb(75) / arb(4) * pi * n2 * e5h

    # E_n and derivatives
    E   = (-pi * n2 * e2u).exp()
    Ed1 = -two * pi * n2 * e2u * E
    Ed2 = (-arb(4) * pi * n2 * e2u + arb(4) * pi2 * n4 * e4u) * E

    # Product rule
    phi    = g * E
    phi_d1 = gd1 * E + g * Ed1
    phi_d2 = gd2 * E + two * gd1 * Ed1 + g * Ed2

    return phi, phi_d1, phi_d2


# ---------------------------------------------------------------------------
# Tail majorants (explicit, derivative-aware)
#
# For n >= 11 and u in [0.99, 3.01], with A = e^{9*U_HI/2} (the largest
# growth factor e^{9u/2} on the cover) and D_n = exp(-pi*n^2*e^{2*U_LO})
# (the weakest decay on the cover):
#   |g_n|   <= (2 pi^2 n^4 + 3 pi n^2) A            <= 30 n^4 A
#   |g_n'|  <= (9 pi^2 n^4 + (15/2) pi n^2) A       <= 100 n^4 A
#   |g_n''| <= ((81/2) pi^2 n^4 + (75/4) pi n^2) A  <= 460 n^4 A
#   |E_n|   <= D_n
#   |E_n'|  <= 2 pi n^2 e^{2 U_HI} D_n
#   |E_n''| <= (4 pi n^2 e^{2 U_HI} + 4 pi^2 n^4 e^{4 U_HI}) D_n
# Product rule then gives |phi_n^{(k)}| <= C_k * n^{4+2k} * A * D_n with
#   C_0 = 30,
#   C_1 = 100 + 30 * 2 pi e^{2 U_HI}                       (< 7.8e4),
#   C_2 = 460 + 230 * 2 pi e^{2 U_HI}
#         + 30 * (4 pi e^{2 U_HI} + 4 pi^2 e^{4 U_HI})     (< 2.1e8),
# evaluated below in Arb (the code uses the Arb-evaluated constants, not
# the rounded values in this comment).  The majorant ratio satisfies
#   M_k(n+1)/M_k(n) <= ((n+1)/n)^{4+2k} * exp(-pi*(2n+1)*e^{2 U_LO}),
# which the code certifies to be < 1/2 at n = N_TAIL_SUM; the tail beyond
# N_TAIL_SUM is therefore at most twice the first neglected majorant.
# ---------------------------------------------------------------------------

def _majorant_constants():
    pi = arb.pi()
    e2uhi = (2 * U_HI).exp()
    c0 = arb(30)
    c1 = arb(100) + arb(30) * 2 * pi * e2uhi
    c2 = (arb(460) + arb(230) * 2 * pi * e2uhi
          + arb(30) * (4 * pi * e2uhi + 4 * pi * pi * e2uhi * e2uhi))
    A = (arb(9) * U_HI / 2).exp()
    return (c0, c1, c2), A


def _majorant(k, n, consts, A):
    """M_k(n) = C_k * n^{4+2k} * A * exp(-pi*n^2*e^{2 U_LO}), an arb >= 0."""
    pi = arb.pi()
    na = arb(n)
    decay = (-pi * na * na * (2 * U_LO).exp()).exp()
    return consts[k] * na ** (4 + 2 * k) * A * decay


def tail_allowances_arb():
    """Arb upper bounds (tau_0, tau_1, tau_2) on the n > 10 tail sums
    |sum_{n>10} phi_n^{(k)}(u)|, valid for all u in [0.99, 3.01],
    and the eps* tail bound for n > 15.  All returned values are arbs
    whose (point) value is a rigorous upper bound."""
    consts, A = _majorant_constants()
    pi = arb.pi()
    taus = []
    for k in range(3):
        total = arb(0)
        for n in range(N_FINITE + 1, N_TAIL_SUM + 1):
            total += _majorant(k, n, consts, A)
        # Certified ratio bound at n = N_TAIL_SUM.
        n = N_TAIL_SUM
        r = ((arb(n + 1) / arb(n)) ** (4 + 2 * k)
             * (-pi * (2 * n + 1) * (2 * U_LO).exp()).exp())
        assert r < arb(1) / 2, "majorant ratio bound not < 1/2"
        total += _majorant(k, n + 1, consts, A) / (1 - r)
        taus.append(total.upper())
    # eps* tail: 2 * sum_{n>15} n^4 e^{-pi*(n^2-1)*e^{2u}}
    #          <= 2 * sum_{n>15} n^4 e^{-pi*(n^2-1)*e^{2 U_LO}}.
    e2ulo = (2 * U_LO).exp()
    etotal = arb(0)
    for n in range(N_EPS_TERMS + 1, N_TAIL_SUM + 1):
        na = arb(n)
        etotal += 2 * na ** 4 * (-pi * (na * na - 1) * e2ulo).exp()
    n = N_TAIL_SUM
    r = ((arb(n + 1) / arb(n)) ** 4
         * (-pi * (2 * n + 1) * e2ulo).exp())
    assert r < arb(1) / 2, "eps* majorant ratio bound not < 1/2"
    na = arb(n + 1)
    etotal += (2 * na ** 4 * (-pi * (na * na - 1) * e2ulo).exp()) / (1 - r)
    return taus[0], taus[1], taus[2], etotal.upper()


# Tail allowances are u-independent (global majorants on the whole cover),
# so they are computed once.
TAU0, TAU1, TAU2, EPS_TAIL = (None, None, None, None)


def epsilon_upper_arb(u_ball, n_terms=N_EPS_TERMS):
    """Arb upper bound on eps*(u) = 2*sum_{n>=2} n^4 * exp(-pi*(n^2-1)*e^{2u}).

    Finite ball evaluation for n = 2..n_terms plus the explicit majorant
    tail bound EPS_TAIL for n > n_terms (computed in tail_allowances_arb).
    """
    pi = arb.pi()
    two = arb(2)
    e2u = (two * u_ball).exp()
    total = arb(0)
    for n in range(2, n_terms + 1):
        n_arb = arb(n)
        n2 = n_arb * n_arb
        n4 = n2 * n2
        term = two * n4 * (-pi * (n2 - arb(1)) * e2u).exp()
        total = total + term
    return total + EPS_TAIL


def certify_checkpoint_arb(center, radius):
    """Certify (log Phi)''(u) < 0 on [center-radius, center+radius].

    Direct enclosure of DeltaQ (Appendix E.5 as executed):
      B(u) = |DeltaQ(u)| / phi_1(u)^2 + 3 * eps*(u) * |W_1(u)|
    bounds |W_tail(u)| (Lemma Wtail_bound).  DeltaQ is evaluated from its
    definition with R split at n = 10; the n > 10 tail enters through the
    explicit allowances TAU0/TAU1/TAU2.  The sign test is the Arb ball
    comparison (W_1 + B) < 0.

    Returns (ok, margin, W1_upper, B_upper, B_str): floats for display,
    plus B as a decimal string (B underflows float at large u).
    """
    u_ball = make_ball(center, radius)

    # W_1(u) — both terms negative, no cancellation
    W1 = log_phi1_d2_arb(u_ball)

    # phi_1 and derivatives (balls)
    phi1, phi1d1, phi1d2 = phi_n_all_arb(u_ball, 1)
    assert phi1 > 0, "phi_1 not certifiably positive on this interval"

    # Finite tail sums A = sum_{n=2}^{10} phi_n (balls)
    A = Ad1 = Ad2 = arb(0)
    for n in range(2, N_FINITE + 1):
        p, pd1, pd2 = phi_n_all_arb(u_ball, n)
        A += p
        Ad1 += pd1
        Ad2 += pd2

    # DeltaQ from the finite part (ball enclosure)
    dQ = Ad2 * phi1 + phi1d2 * A - 2 * phi1d1 * Ad1 + Ad2 * A - Ad1 ** 2

    # |DeltaQ| <= |dQ| + tail cross terms.  Writing the full tail sums as
    # T, T', T'' with |T^{(k)}| <= TAU_k, the difference DeltaQ - dQ is
    #   T''*phi_1 + phi_1''*T - 2*phi_1'*T'
    #   + (T''*A + A''*T - 2*A'*T') + (T''*T - (T')^2),
    # bounded term by term with abs_upper (rigorous point bounds); the
    # final .upper() makes the result a single certified point bound.
    abs_dQ = (dQ.abs_upper()
              + TAU2 * phi1.abs_upper() + phi1d2.abs_upper() * TAU0
              + 2 * phi1d1.abs_upper() * TAU1
              + TAU2 * A.abs_upper() + Ad2.abs_upper() * TAU0
              + 2 * Ad1.abs_upper() * TAU1
              + TAU2 * TAU0 + TAU1 * TAU1).upper()

    # eps*(u) upper bound (ball + explicit tail)
    eps_up = epsilon_upper_arb(u_ball).upper()

    # B(u) as a rigorous point bound (arb), then the Arb sign test.
    phi1_sq_lower = (phi1.lower() * phi1.lower()).lower()
    B = (abs_dQ / phi1_sq_lower + 3 * eps_up * W1.abs_upper()).upper()
    total = W1 + B
    ok = bool(total < 0)
    margin = -upper(total)
    return ok, margin, upper(W1), float(B.upper()), fmt(B)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    global TAU0, TAU1, TAU2, EPS_TAIL
    ctx.prec = PREC
    TAU0, TAU1, TAU2, EPS_TAIL = tail_allowances_arb()

    # Adaptive grid covering [0.99, 3.01]: half-width at left endpoint u is
    # min(HW_MAX, EXP_HW / (2*pi*e^{2u})), so the dominant term's exponent
    # -pi*e^{2u} varies by at most ~EXP_HW across each ball and the ball
    # evaluation of phi_1 stays positive and tight on the whole range.
    import math
    grid = []
    left = 0.99
    while left < 3.01 - 1e-12:
        hw = min(HW_MAX, EXP_HW / (2 * math.pi * math.exp(2 * left)))
        if left + 2 * hw > 3.01:
            hw = (3.01 - left) / 2
        grid.append((left + hw, hw))
        left += 2 * hw
    n_grid = len(grid)

    print("=" * 72)
    print("  ARB/FLINT CERTIFICATE: (log Phi)''(u) < 0 on [1.0, 3.0]")
    print(f"  Method: Arb ball arithmetic, {PREC}-bit precision")
    print("  Direct DeltaQ enclosure (finite n<=10) + explicit tail majorants")
    print(f"  {n_grid} adaptive checkpoints covering [0.99, 3.01], "
          f"half-width <= {HW_MAX}")
    print(f"  Tail allowances: tau_0={fmt(TAU0)}, tau_1={fmt(TAU1)}, "
          f"tau_2={fmt(TAU2)}, eps_tail={fmt(EPS_TAIL)}")
    print("=" * 72)
    print()
    print(f"  {'u':>6}  {'W1_upper':>14}  {'B_upper':>14}  {'margin':>14}  status")
    print("  " + "-" * 62)

    cert_count = failed_count = 0
    min_margin = float("inf")
    worst_u = 1.0
    max_B_upper = 0.0
    max_sign_upper = -float("inf")
    interval_results = []
    t0 = time.time()

    for i, (u, hw) in enumerate(grid):
        ok, margin, W1_up, B_up, B_str = certify_checkpoint_arb(u, hw)

        sign_up = W1_up + B_up

        if ok:
            cert_count += 1
            if margin < min_margin:
                min_margin = margin
                worst_u = u
        else:
            failed_count += 1

        max_B_upper = max(max_B_upper, B_up)
        max_sign_upper = max(max_sign_upper, sign_up)

        status = "OK" if ok else "*** FAIL ***"
        print(f"  {u:6.3f}  {W1_up:>14.4e}  {B_up:>14.4e}  {margin:>14.4e}  {status}")

        interval_results.append({
            "i": i,
            "center": f"{u:.6f}",
            "radius": f"{hw:.6f}",
            "B_upper": B_str,
            "sign_upper": f"{sign_up:.6e}",
            "margin": f"{margin:.6e}",
        })

    elapsed = time.time() - t0
    all_ok = failed_count == 0

    print()
    print("=" * 72)
    print(f"  Certified: {cert_count} / {n_grid}")
    print(f"  Failed:    {failed_count} / {n_grid}")
    print(f"  Min margin: {min_margin:.4e} at u={worst_u:.4f}")
    print(f"  Max B_upper: {max_B_upper:.4e}")
    print(f"  Time: {elapsed:.1f}s")

    if all_ok:
        print()
        print(f"  *** ALL {n_grid} CHECKPOINTS CERTIFIED ***")
        print("  (log Phi)''(u) < 0 on [1.0, 3.0] — Arb/FLINT primary certificate")
    else:
        print("  *** FAILURES — REVIEW ABOVE ***")
    print("=" * 72)

    result = {
        "claim": "(log Phi)''(u) < 0 on [1,3]",
        "method": "Arb/FLINT ball arithmetic via python-flint; direct DeltaQ "
                  "enclosure (finite sum n=1..10) plus explicit derivative-aware "
                  "tail majorants for n>10; sign test (W_1 + B) < 0 in Arb",
        "coverage": "[0.99,3.01]",
        "intervals_checked": n_grid,
        "grid": "adaptive: half-width min(0.01, 0.35/(2*pi*e^{2u})) at the "
                "left endpoint, contiguous cover of [0.99, 3.01]",
        "precision_bits": PREC,
        "finite_terms": N_FINITE,
        "finite_terms_eps_star": N_EPS_TERMS,
        "tail_allowance_B": fmt(TAU2),
        "tail_bound_n_gt_10": f"|R''| tail < {fmt(TAU2)}; "
                              f"eps* tail < {fmt(EPS_TAIL)}",
        "max_B_upper": f"{max_B_upper:.6e}",
        "max_sign_upper": f"{max_sign_upper:.6e}",
        "minimum_margin": f"{min_margin:.6e}",
        "worst_u": worst_u,
        "passed": all_ok,
        "certified": cert_count,
        "failed": failed_count,
        "time_s": round(elapsed, 2),
        "interval_results": interval_results,
    }

    out_path = os.path.join(ROOT, "results", "verify_ia_1_to_3_arb.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  -> results/verify_ia_1_to_3_arb.json")
    return result


if __name__ == "__main__":
    result = main()
    sys.exit(0 if result["passed"] else 1)
