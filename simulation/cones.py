"""The two support cones and their comparison.

For a dominant real cocharacter chi of G put x_i = <chi, alpha_i>, z = <chi, theta_G>, R = n.x and
W_i = (x . M[:, i]) / m_i = <chi, varpi_i> / m_i. With kappa = <rho, rho^vee> > 1:

  automorphic cone (A):  x >= 0,  R/(kappa-1) >= z >= max( (R - x_j)/(kappa-1), W_i )      [Chen-Venkatesh (66)]
  spectral cone (S_d):   R >= z/d >= R - min_j x_j                                      [Chen-Venkatesh (70)]

where d = <theta_G, nu_G> is the primitive central pairing (theta_G <-> d theta_vee of the dual group).
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product

import numpy as np
from scipy.optimize import linprog

import rootdata as rd


def _R(D, x):
    return sum(a * b for a, b in zip(D["n"], x))


def _W(D, x, i):
    r = D["rank"]
    return sum(x[j] * D["M"][j][i] for j in range(r)) / D["m"][i]


def in_A(D, x, z):
    k = D["kappa"]
    if any(v < 0 for v in x):
        return False
    R = _R(D, x)
    lower = max([(R - xj) / (k - 1) for xj in x] + [_W(D, x, i) for i in range(D["rank"])])
    return R / (k - 1) >= z >= lower


def in_S(D, x, z, d):
    R = _R(D, x)
    return R >= F(z) / d >= R - min(x)


def in_58(D, x, z):
    """The simplified cone (58): R/(kappa-1) >= z >= max_j (R - x_j)/(kappa-1), x >= 0."""
    k = D["kappa"]
    if any(v < 0 for v in x):
        return False
    R = _R(D, x)
    return R / (k - 1) >= z >= max((R - xj) / (k - 1) for xj in x)


def slice_intervals(D, d):
    """z-intervals of the two cones on the slice x = (1, ..., 1)."""
    k = D["kappa"]
    x = [F(1)] * D["rank"]
    R = _R(D, x)
    lowA = max([(R - xj) / (k - 1) for xj in x] + [_W(D, x, i) for i in range(D["rank"])])
    return {"A": [lowA, R / (k - 1)], "S": [d * (R - 1), d * R]}


def barycentric(D):
    """Matrix B with columns B[:, i] = n - (kappa-1)/m_i M[:, i]: varpi_i/m_i = sum_j B_ji (rho - alpha_j)/(kappa-1)."""
    r, k = D["rank"], D["kappa"]
    return [[D["n"][j] - (k - 1) / D["m"][i] * D["M"][j][i] for i in range(r)] for j in range(r)]


def _rat(v, den=10 ** 6):
    return F(v).limit_denominator(den)


def redundancy(D):
    """Decide exactly whether every W_i <= max_j f_j on the orthant, f_j = (R - x_j)/(kappa-1).

    For each i solve  max t  s.t.  W_i(x) - f_j(x) >= t for all j, x >= 0, sum x = 1.
    If the optimum is <= 0 the inequality z >= W_i is redundant; an exact certificate is a convex
    combination lambda with sum_j lambda_j f_j - W_i having nonnegative coefficients. If > 0 a
    rational witness x with W_i(x) > f_j(x) for all j is returned.
    """
    r, k = D["rank"], D["kappa"]
    n, M, m = D["n"], D["M"], D["m"]
    fcoef = [[(n[l] - (1 if l == j else 0)) / (k - 1) for l in range(r)] for j in range(r)]
    out = []
    for i in range(r):
        wcoef = [M[l][i] / m[i] for l in range(r)]
        # variables (x_1..x_r, t); maximize t
        A_ub, b_ub = [], []
        for j in range(r):
            A_ub.append([-(float(wcoef[l] - fcoef[j][l])) for l in range(r)] + [1.0])
            b_ub.append(0.0)
        res = linprog(c=[0.0] * r + [-1.0], A_ub=A_ub, b_ub=b_ub, A_eq=[[1.0] * r + [0.0]], b_eq=[1.0],
                      bounds=[(0, None)] * r + [(None, None)], method="highs")
        t = -res.fun
        if t > 1e-9:
            x = [_rat(v) for v in res.x[:r]]
            ok = all(sum(wcoef[l] * x[l] for l in range(r)) > sum(fcoef[j][l] * x[l] for l in range(r)) for j in range(r))
            out.append({"i": i, "redundant": False, "witness": [str(v) for v in x], "exact": ok})
        else:
            lam = [max(0.0, -v) for v in res.ineqlin.marginals]
            s = sum(lam)
            lam = [_rat(v / s) for v in lam] if s > 0 else [F(1, r)] * r
            lam[-1] = 1 - sum(lam[:-1])
            coef = [sum(lam[j] * fcoef[j][l] for j in range(r)) - wcoef[l] for l in range(r)]
            ok = all(v >= 0 for v in lam) and all(v >= 0 for v in coef)
            if not ok:  # fall back to the barycentric certificate when it exists
                B = barycentric(D)
                lam = [B[j][i] for j in range(r)]
                coef = [sum(lam[j] * fcoef[j][l] for j in range(r)) - wcoef[l] for l in range(r)]
                ok = all(v >= 0 for v in lam) and all(v >= 0 for v in coef)
            out.append({"i": i, "redundant": True, "certificate": [str(v) for v in lam], "exact": ok})
    return out


def classify(max_rank: int = 8):
    """Every semisimple type of total rank <= max_rank with kappa > 1."""
    rows = []
    for comps in rd.products(max_rank):
        D = rd.data(rd.cartan_of(comps))
        if D["kappa"] <= 1:
            continue
        d_req = 1 / (D["kappa"] - 1)
        B = barycentric(D)
        convex = all(B[j][i] >= 0 for i in range(D["rank"]) for j in range(D["rank"]))
        red = redundancy(D)
        rows.append({
            "type": rd.name(comps), "rank": D["rank"], "kappa": str(D["kappa"]), "H": D["H"],
            "required_d": str(d_req), "required_d_integral": d_req.denominator == 1,
            "convex_hull": convex, "redundant": all(x["redundant"] for x in red),
            "certificates_exact": all(x["exact"] for x in red),
            "data_checks": D["check_n"] and D["check_m"] and D["check_kappa"],
        })
    return rows


def lattice_compare(D, d, bound, zmax):
    """Compare (A) and (S_d) on integer x in [0, bound]^r and integer 0 <= z <= zmax."""
    agree = total = 0
    sep = None
    for x in product(range(bound + 1), repeat=D["rank"]):
        x = [F(v) for v in x]
        for z in range(zmax + 1):
            a, s = in_A(D, x, z), in_S(D, x, z, d)
            total += 1
            agree += a == s
            if a != s and sep is None:
                sep = {"x": [int(v) for v in x], "z": z, "in_A": a, "in_S": s}
    return {"points": total, "agree": agree, "first_disagreement": sep}


def node_bound(D):
    """q(l) = n_l m_l / M_ll. Testing x = e_l shows that redundancy of z >= W_l requires kappa - 1 <= q(l)."""
    q = [D["n"][l] * D["m"][l] / D["M"][l][l] for l in range(D["rank"])]
    return {"q": [str(v) for v in q], "q_min": min(q), "kappa_minus_1": D["kappa"] - 1, "passes": D["kappa"] - 1 <= min(q)}


def node_bound_irreducible():
    """Node bound for every irreducible type of rank <= 8, and the closed forms used in the proof for node 1."""
    rows = {}
    for kind, r in rd.IRREDUCIBLE:
        D = rd.data(rd.cartan(kind, r))
        nb = node_bound(D)
        q1 = D["n"][0] * D["m"][0] / D["M"][0][0]
        closed = {"A": F(r * (r + 1), 4), "B": F(r * (2 * r - 1), 2), "C": F(r * (2 * r - 1), 2), "D": F((r - 1) ** 2)}.get(kind)
        kappa_closed = {"A": F(r * (r + 1) * (r + 2), 12), "B": F(r * (r + 1) * (4 * r - 1), 12),
                        "C": F(r * (r + 1) * (4 * r - 1), 12), "D": F(r * (r - 1) * (2 * r - 1), 6)}.get(kind)
        rows[f"{kind}{r}"] = {"kappa": str(D["kappa"]), "q_min": str(nb["q_min"]), "q_node1": str(q1),
                              "node_bound_passes": nb["passes"],
                              "q_node1_closed_form_ok": None if closed is None else q1 == closed,
                              "kappa_closed_form_ok": None if kappa_closed is None else D["kappa"] == kappa_closed}
    return rows
