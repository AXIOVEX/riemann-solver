# ROUND-2 ADJUDICATION CONSEQUENCE RUN: identical machinery to e2cert2.py
# except the Sturm coefficients use gamma_{n+k} = (n+k)! * a_{n+k} (the
# standard xi Jensen family), from the same interval moments.
# E2 CERTIFICATION EXTENSION: extend the CERTIFIED Jensen grid for the
# Phi kernel to higher degrees and shifts, reusing the E2 certification
# machinery (runs/e2/e2cert.py) UNCHANGED IN METHOD.
#
# REGISTERED EXTENSION (fixed in this header before running; the grid is
# NOT extended further after seeing results):
#   d = 2..12, n = 0..20  --  11 x 21 = 231 grid points total.
#   CONTROL BLOCK: the 77 points (d = 2..8, n = 0..10) already CERTIFIED
#   by run e2cert-20261122T145019 are re-run here under THIS run's
#   independently computed interval moments; they must ALL certify again
#   or the batch is NO VERDICT. New territory: the remaining 154 points.
#   This run needs interval moments a_0..a_32.
#
# Definitions (paper's normalization, identical to E1/E2/E2-cert):
#   Phi(u) = sum_{n=1}^{6} (2 pi^2 n^4 e^{4.5u} - 3 pi n^2 e^{2.5u})
#                          * e^{-pi n^2 e^{2u}}
#   a_k = int_0^inf Phi(u) u^{2k} du / (2k)!
#   J^{d,n}(x) = sum_{k=0}^{d} C(d,k) a_{n+k} x^k
#
# STEP 1 -- INTERVAL MOMENTS. a_0..a_32 are enclosed rigorously:
#   (a) int_0^8 Phi(u) u^{2k} du by FLINT's acb_integral (python-flint
#       acb.integral): adaptive Gauss-Legendre in ball arithmetic; the
#       returned acb ball is a guaranteed enclosure of the integral of
#       the integrand as evaluated in ball arithmetic (the integrand is
#       entire in u, so no branch-cut handling arises). The true
#       integral is real and lies in the ball's real part.
#   (b) Tail int_8^inf, bounded analytically (derivation below) and
#       added as a ball containing [0, tail_bound] before dividing by
#       (2k)! (the bound's derivation does not depend on k except
#       through the u^{2k} majorant; the (2k)! division is applied per
#       k, as in the cert run).
#
#   TAIL BOUND DERIVATION (restated for KMAX = 32). For u >= 0 and each
#   n in 1..6,
#     term_n(u) = (2 pi^2 n^4 e^{9u/2} - 3 pi n^2 e^{5u/2})
#                 * e^{-pi n^2 e^{2u}}
#               <= 2 pi^2 n^4 e^{9u/2} e^{-pi e^{2u}},
#   because the subtracted term is >= 0 and n^2 >= 1 gives
#   e^{-pi n^2 e^{2u}} <= e^{-pi e^{2u}}. Summing (sum_{n=1}^6 n^4 =
#   2275):
#     Phi(u) <= 4550 pi^2 e^{9u/2 - pi e^{2u}}.
#   For u >= 8 and k <= 32, u^{2k} <= u^{64}; set
#     g(u) = (9/2) u + 64 ln u - pi e^{2u}.
#   g'(u) = 9/2 + 64/u - 2 pi e^{2u} <= 9/2 + 8 - 2 pi e^{16} =: -lambda
#   for u >= 8, with lambda = 2 pi e^{16} - 12.5 > 5.5e7 > 0. Hence
#   g(u) <= g(8) - lambda (u - 8) and
#     int_8^inf Phi(u) u^{2k} du <= 4550 pi^2 e^{g(8)} / lambda,
#     g(8) = 36 + 64 ln 8 - pi e^{16} ~ -2.79e7,
#   so the bound is ~10^{-1.21e7}: the tail is zero at any working
#   precision, but it is BOUNDED, not asserted away; the script computes
#   the bound in ball arithmetic and adds the ball [0, bound] (as a
#   containing ball) to every moment enclosure before dividing by (2k)!.
#
# STEP 2 -- STURM CERTIFICATION WITH BALL COEFFICIENTS (identical to
# the cert run). For each J^{d,n}, the Sturm sequence is computed by a
# FIXED straight-line program (polynomial long division) in ball
# arithmetic from the coefficient balls C(d,k) * a_{n+k}. By inclusion
# isotonicity of ball arithmetic, if every divisor leading-coefficient
# ball and every Sturm leading-coefficient ball excludes 0, the
# computed ball sequence encloses the exact Sturm sequence of EVERY
# polynomial whose coefficients lie in the moment balls, all members
# share the same degree pattern, and the signs of all leading
# coefficients are the same for every member. Sign variations at
# +/-infinity are read off the leading-coefficient signs (times
# (-1)^degree at -infinity). V(-inf) - V(+inf) is then the exact
# number of distinct real roots of every member. The final Sturm
# element is a constant ball; if it excludes 0, every member is
# separable (gcd(p, p') constant). So:
#   count == d and final constant excludes 0  =>  CERTIFIED
#   count <  d and all signs determinate      =>  FAILED (falsifier:
#       the true J^{d,n} is not hyperbolic)
#   any sign indeterminate at the precision cap => INDETERMINATE.
# Secondary consistency check: the discriminant, via the Sylvester
# determinant in ball arithmetic; if its sign is determinate it must
# be > 0 for a CERTIFIED polynomial (distinct real roots); a
# determinate disc <= 0 alongside a Sturm count of d would be an
# internal inconsistency and demotes the polynomial to INDETERMINATE
# with a flag. The discriminant ball itself is recorded per polynomial
# (mid, radius, relative width, sign margin) as the carry-forward
# margin dataset.
#
# PRECISION ESCALATION: levels 256, 512, 1024, 2048 bits (cap 2048),
# identical ladder to the cert run. Moments are recomputed at each
# level (rel_tol = 2^{-(prec-64)}). Every attempt is logged.
#
# CROSS-CHECKS vs numerical records: E2 payload stored moments a_0..18
# rounded to 40 significant digits (100 dps run); E1 payload stored
# moments a_0..41 at 80 dps in moments_partial.json. Each stored value
# is treated as a ball of radius half a unit in its last stored digit,
# and consistency = that ball intersects our enclosure. k = 0..18 are
# checked against BOTH records; k = 19..32 against E1. Max enclosure
# relative width is reported per level.
#
# Payload: results_e2gamma.json in the run directory. The verdict comes
# from this file, never from the wrapper (Spec 002 TC-03).
import json, math, os, sys, time, traceback
from flint import arb, acb, arb_mat, ctx

RUN_DIR = os.path.dirname(os.path.abspath(__file__))
E2_PAYLOAD = os.path.expanduser(
    "~/riemann-runs/e2-20261122T143947/results_e2.json")
E1_MOMENTS = os.path.expanduser(
    "~/riemann-runs/e1-20261122/moments_partial.json")
KMAX = 32
GRID = [(d, n) for d in range(2, 13) for n in range(0, 21)]
CONTROL_BLOCK = [(d, n) for d in range(2, 9) for n in range(0, 11)]
CONTROLS_TC02 = [(2, 0), (3, 0), (3, 1), (4, 0)]
LEVELS = [256, 512, 1024, 2048]


def Phi_acb(u):
    s = acb(0)
    pi = acb.pi()
    for n in range(1, 7):
        s += (2 * pi**2 * n**4 * (acb("4.5") * u).exp()
              - 3 * pi * n**2 * (acb("2.5") * u).exp()) \
            * (-pi * n**2 * (2 * u).exp()).exp()
    return s


def tail_ball():
    """Ball T with [0, 4550 pi^2 e^{g(8)} / lambda] subset of T, plus
    the bound ball itself (see header derivation, KMAX = 32:
    g(u) = (9/2) u + 64 ln u - pi e^{2u}, lambda = 2 pi e^{16} - 12.5,
    g(8) = 36 + 64 ln 8 - pi e^{16})."""
    pi = arb.pi()
    e16 = arb(16).exp()
    g8 = arb(36) + arb(2 * KMAX) * arb(8).log() - pi * e16
    lam = 2 * pi * e16 - (arb(9) / 2 + arb(2 * KMAX) / 8)
    assert lam > 0
    bound = arb(4550) * pi**2 * g8.exp() / lam
    assert bound > 0
    T = arb(0).union(bound)  # contains 0 and all of bound => contains [0, sup]
    return T, bound


def compute_moments(prec):
    ctx.prec = prec
    rel = arb(2) ** (-(prec - 64))
    abst = arb(2) ** (-(prec - 64)) * arb("1e-12")
    T, bound = tail_ball()
    a = []
    for k in range(KMAX + 1):
        def f(u, analytic, k=k):
            return Phi_acb(u) * u ** (2 * k)
        r = acb.integral(f, acb(0), acb(8), rel_tol=rel, abs_tol=abst,
                         eval_limit=10**7, depth_limit=1000)
        assert r.imag.contains(arb(0)) or r.imag == 0
        a.append((r.real + T) / arb(math.factorial(2 * k)))
        print(f"  moment a_{k} done (prec {prec})", flush=True)
    return a, bound


def sign_of(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def poly_rem(A, B):
    """Enclosure of A mod B (lists of arb, ascending). Requires the
    leading coefficient ball of B to exclude 0 (checked by caller)."""
    dB = len(B) - 1
    R = list(A)
    if len(R) < len(B):
        R = R + [arb(0)] * (len(B) - len(R))
    leadB = B[-1]
    for k in range(len(A) - 1 - dB, -1, -1):
        factor = R[dB + k] / leadB
        for j in range(dB + 1):
            R[k + j] = R[k + j] - factor * B[j]
    return R[:dB]


def discriminant_ball(coeffs):
    """disc(p) as an arb ball via the Sylvester determinant in ball
    arithmetic. Returns the ball, or None if unavailable."""
    d = len(coeffs) - 1
    if d < 2:
        return None
    try:
        P = [coeffs[d - i] for i in range(d + 1)]            # descending
        Q = [(d - j) * coeffs[d - j] for j in range(d)]      # p' descending
        size = 2 * d - 1
        rows = []
        for i in range(d - 1):
            rows.append([arb(0)] * i + P + [arb(0)] * (d - 2 - i))
        for i in range(d):
            rows.append([arb(0)] * i + Q + [arb(0)] * (d - 1 - i))
        res = arb_mat(rows).det()
        disc = res / coeffs[-1]
        if (d * (d - 1) // 2) % 2 == 1:
            disc = -disc
        return disc
    except Exception:
        return None


def disc_record(disc):
    """Serializable record of a discriminant ball: sign, and the
    margin dataset (mid, radius, relative width, sign margin =
    (|mid| - rad)/|mid|, log10 |mid|)."""
    if disc is None:
        return {"sign": 0, "available": False}
    rec = {"sign": sign_of(disc), "available": True, "ball": str(disc)}
    try:
        mid = disc.mid()
        rad = disc.rad()
        rec["mid"] = float(mid)
        rec["rad"] = float(rad)
        if mid != 0:
            rec["rel_width"] = float(rad / abs(mid))
            rec["sign_margin"] = float((abs(mid) - rad) / abs(mid))
            rec["log10_abs_mid"] = float(
                abs(mid).log() / arb(10).log())
    except Exception:
        pass
    return rec


def sturm_attempt(coeffs):
    """One Sturm certification attempt at the current precision.
    coeffs: ascending arb list for J^{d,n}. Identical logic to the
    cert run, plus the discriminant ball record."""
    d = len(coeffs) - 1
    P0 = list(coeffs)
    P1 = [(j + 1) * coeffs[j + 1] for j in range(d)]
    for P in (P0, P1):
        if sign_of(P[-1]) == 0:
            return {"status": "indeterminate",
                    "reason": "leading coefficient sign of J or J'"}
    seq = [P0, P1]
    while True:
        A, B = seq[-2], seq[-1]
        if len(B) == 1:
            break
        if sign_of(B[-1]) == 0:
            return {"status": "indeterminate",
                    "reason": "divisor leading coefficient sign"}
        R = poly_rem(A, B)
        if len(R) == 0:
            return {"status": "indeterminate",
                    "reason": "empty remainder"}
        if sign_of(R[-1]) == 0:
            return {"status": "indeterminate",
                    "reason": "remainder leading coefficient sign",
                    "sturm_length_so_far": len(seq)}
        seq.append([-c for c in R])
        if len(R) == 1:
            break
    final = seq[-1]
    if len(final) != 1 or sign_of(final[0]) == 0:
        return {"status": "indeterminate",
                "reason": "final Sturm constant sign (separability)"}

    def variations(at_minus):
        signs = []
        for P in seq:
            s = sign_of(P[-1])
            if s == 0:
                return None
            if at_minus and ((len(P) - 1) % 2 == 1):
                s = -s
            signs.append(s)
        return sum(1 for i in range(len(signs) - 1)
                   if signs[i] != signs[i + 1])

    vm, vp = variations(True), variations(False)
    if vm is None or vp is None:
        return {"status": "indeterminate",
                "reason": "sign variation evaluation"}
    count = vm - vp
    drec = disc_record(discriminant_ball(coeffs))
    rec = {"sturm_length": len(seq), "real_root_count": count,
           "final_sturm_constant": str(final[0]),
           "discriminant_sign": drec["sign"],
           "discriminant": drec}
    if count == d:
        if drec["sign"] < 0:
            rec.update({"status": "indeterminate",
                        "reason": "INTERNAL INCONSISTENCY: Sturm count "
                                  "= d but discriminant <= 0"})
            return rec
        rec["status"] = "certified"
        return rec
    if count < d:
        rec["status"] = "failed"
        return rec
    rec.update({"status": "indeterminate",
                "reason": f"root count {count} > degree {d} (impossible; "
                          "enclosure failure)"})
    return rec


def stored_ball(s):
    mant = s.lower().split('e')[0].replace('-', '').replace('+', '') \
        .replace('.', '').lstrip('0')
    sig = len(mant) if mant else 1
    v = arb(s)
    rho = v * arb(10) ** (-(sig - 1)) / 2
    return (v - rho).union(v + rho)


def ballrec(x):
    return {"ball": str(x), "mid": float(x.mid()), "rad": float(x.rad()),
            "rel_width": float(x.rad() / abs(x.mid())) if x.mid() != 0
            else None}


def tally(recs):
    return {
        "certified": sum(1 for r in recs if r["status"] == "certified"),
        "failed": sum(1 for r in recs if r["status"] == "failed"),
        "indeterminate": sum(1 for r in recs
                             if r["status"] == "indeterminate"),
        "total": len(recs),
    }


def main():
    selftest = len(sys.argv) > 1 and sys.argv[1] == "selftest"
    t_start = time.time()
    with open(E2_PAYLOAD) as fh:
        e2 = json.load(fh)
    e2_moments = e2["moments_a"]
    with open(E1_MOMENTS) as fh:
        e1_moments = json.load(fh)["a"]

    if selftest:
        grid = [(2, 0), (3, 0), (3, 1), (4, 0), (8, 10), (12, 0),
                (12, 20)]
        levels = LEVELS[:1]
    else:
        grid = GRID
        levels = LEVELS

    decided = {}
    attempts = []
    level_info = []
    moments_final = None
    tail_bound_rec = None

    for level, prec in enumerate(levels):
        if all(k in decided for k in grid):
            break
        t0 = time.time()
        print(f"LEVEL {level}: prec {prec} bits", flush=True)
        a, bound = compute_moments(prec)
        moments_final = a
        tail_bound_rec = {"ball": str(bound),
                          "log10_approx": float(bound.log()
                                                / arb(10).log())}
        # cross-checks: k=0..18 vs E2 (40 sig digits) AND E1 (80 dps);
        # k=19..32 vs E1 only.
        cross_e2, cross_e1 = {}, {}
        max_rel_width = 0.0
        for k in range(KMAX + 1):
            if k <= 18:
                cross_e2[str(k)] = bool(
                    a[k].overlaps(stored_ball(e2_moments[str(k)])))
            cross_e1[str(k)] = bool(
                a[k].overlaps(stored_ball(e1_moments[k])))
            rw = float(a[k].rad() / abs(a[k].mid()))
            max_rel_width = max(max_rel_width, rw)
        level_info.append({
            "level": level, "prec_bits": prec,
            "seconds": time.time() - t0,
            "max_moment_rel_width": max_rel_width,
            "cross_check_vs_E2_all_intersect": all(cross_e2.values()),
            "cross_check_vs_E2": cross_e2,
            "cross_check_vs_E1_all_intersect": all(cross_e1.values()),
            "cross_check_vs_E1": cross_e1,
        })
        print(f"  cross-checks: E2 {all(cross_e2.values())}, "
              f"E1 {all(cross_e1.values())}, "
              f"max rel width {max_rel_width:.3e}", flush=True)
        for (d, n) in grid:
            if (d, n) in decided:
                continue
            coeffs = [arb(math.comb(d, k)) * arb(math.factorial(n + k))
                      * a[n + k] for k in range(d + 1)]
            try:
                rec = sturm_attempt(coeffs)
            except Exception as e:
                rec = {"status": "indeterminate",
                       "reason": f"exception: {e}"}
            rec.update({"d": d, "n": n, "level": level,
                        "prec_bits": prec})
            attempts.append(rec)
            print(f"  J^{d},{n}: {rec['status']}"
                  + (f" ({rec.get('reason')})"
                     if rec["status"] != "certified" else ""),
                  flush=True)
            if rec["status"] in ("certified", "failed"):
                decided[(d, n)] = rec

    records = {}
    for (d, n) in grid:
        if (d, n) in decided:
            records[f"{d},{n}"] = decided[(d, n)]
        else:
            last = [r for r in attempts
                    if r["d"] == d and r["n"] == n][-1]
            last = dict(last)
            last["status"] = "indeterminate"
            last["cap_bits"] = levels[-1]
            records[f"{d},{n}"] = last

    if selftest:
        label = "SELFTEST ONLY"
        ctrl_ok = None
        tallies_split = None
    else:
        ctrl_recs = [records[f"{d},{n}"] for (d, n) in CONTROL_BLOCK]
        new_recs = [records[f"{d},{n}"] for (d, n) in GRID
                    if (d, n) not in CONTROL_BLOCK]
        ctrl_ok = all(r["status"] == "certified" for r in ctrl_recs)
        tallies_split = {"control_block_77": tally(ctrl_recs),
                         "new_territory_154": tally(new_recs),
                         "full_grid_231": tally(list(records.values()))}
        tt = tallies_split["full_grid_231"]
        if not ctrl_ok:
            label = ("NO VERDICT: the 77-point control block did not "
                     "re-certify in full under this run's independent "
                     "interval moments "
                     f"(control tallies {tallies_split['control_block_77']})")
        elif tt["failed"] > 0:
            failed_named = sorted(k for k, r in records.items()
                                  if r["status"] == "failed")
            label = (f"FALSIFIER: {tt['failed']} polynomial(s) proved "
                     f"NOT hyperbolic ({', '.join(failed_named)}); "
                     f"{tt['certified']} CERTIFIED, "
                     f"{tt['indeterminate']} INDETERMINATE")
        elif tt["certified"] == len(GRID):
            label = ("CERTIFIED for the STANDARD (gamma) grid d=2..12, "
                     "n=0..20 (231/231) -- standard xi Jensen family, "
                     "for this grid only")
        else:
            ind_named = sorted(k for k, r in records.items()
                               if r["status"] == "indeterminate")
            label = (f"CERTIFIED for a subset ({tt['certified']}/231); "
                     f"INDETERMINATE for {tt['indeterminate']} "
                     f"({', '.join(ind_named)}); the pre-existing E2 "
                     "NUMERICAL result (100 dps) stands for the "
                     "indeterminate points")

    out = {
        "experiment": "E2 STANDARD-GRID (gamma) certification (Sturm/interval), round-2 adjudication consequence run, "
                      "Spec 002 R5",
        "classification": "CERTIFIED (ball arithmetic) where stated; "
                          "see label_outcome",
        "venue": "tower-compute",
        "grid_registered": {"d": "2..12", "n": "0..20", "count": 231},
        "control_block": {"d": "2..8", "n": "0..10", "count": 77,
                          "all_certified": ctrl_ok},
        "controls_tc02": [f"{d},{n}" for (d, n) in CONTROLS_TC02],
        "label_outcome": label,
        "tallies": tallies_split if not selftest else
                   tally(list(records.values())),
        "tail_bound": tail_bound_rec,
        "levels": level_info,
        "moments_final_level": {
            str(k): ballrec(moments_final[k]) for k in range(KMAX + 1)},
        "polynomials": records,
        "attempts": attempts,
        "seconds_total": time.time() - t_start,
    }
    name = "results_e2gamma_selftest.json" if selftest \
        else "results_e2gamma.json"
    with open(os.path.join(RUN_DIR, name), "w") as fh:
        json.dump(out, fh, indent=2)
    print("LABEL:", label)
    print("TALLIES:", out["tallies"])
    print("CONTROL BLOCK ALL CERTIFIED:", ctrl_ok)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        err = traceback.format_exc()
        print(err, flush=True)
        try:
            with open(os.path.join(RUN_DIR,
                           "results_e2gamma_ERROR.json"), "w") as fh:
                json.dump({"error": err}, fh, indent=2)
        except Exception:
            pass
        sys.exit(1)
