"""Explicit separating labels and support comparisons for specific global forms."""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product
from math import prod

import cones
import rootdata as rd

A2 = rd.data(rd.cartan("A", 2))
A1_3 = rd.data(rd.cartan_of([("A", 1)] * 3))
A1_4 = rd.data(rd.cartan_of([("A", 1)] * 4))


def gl3():
    """GL3: chi = (a, b, c) gives x = (a - b, b - c), z = a + b + c; d = 3 (det(tI) = t^3)."""
    out = {}
    for chi in ((1, 0, 0), (2, 1, 0)):
        a, b, c = chi
        x, z = [F(a - b), F(b - c)], a + b + c
        out[str(chi)] = {"x": [int(v) for v in x], "z": z, "in_A": cones.in_A(A2, x, z),
                         "in_S_d3": cones.in_S(A2, x, z, 3), "in_S_d1": cones.in_S(A2, x, z, 1)}
    return out


def a1_4_degree_one():
    """(PGL2)^4 x Gm, d = 1: labels at z = 1 in each cone and the summed dimensions of the SL2^4 representations."""
    A, S = [], []
    for x in product(range(6), repeat=4):
        xf = [F(v) for v in x]
        if cones.in_A(A1_4, xf, 1):
            A.append(x)
        if cones.in_S(A1_4, xf, 1, 1):
            S.append(x)
    dimsum = lambda L: sum(prod(v + 1 for v in x) for x in L)
    return {"labels_A": len(A), "labels_S": len(S), "dim_A": dimsum(A), "dim_S": dimsum(S),
            "S_minus_A": sorted(set(S) - set(A)), "A_subset_S": set(A) <= set(S)}


def a1_witness():
    x = [F(2), F(0), F(0), F(0)]
    return {"in_A": cones.in_A(A1_4, x, 1), "in_S_d1": cones.in_S(A1_4, x, 1, 1)}


def gl2_pgl2_pgl2(zmax=8):
    """Identity (1) for G = GL2 x PGL2 x PGL2, dual GL2 x SL2 x SL2.

    Labels chi in X_*(T_G) = X^*(T_dual): coordinates (p1, p2, p3, z), p1 = z mod 2 (GL2 factor), p2, p3 free.
    Automorphic support: chi in cone (A). Spectral support (corrected Proposition 5.1 for the dual group):
    chi in cone (S_2) and <chi, rho_dual - theta_vee_dual> = (p1 + p2 + p3 - z)/2 integral.
    Dimension of the dual-group representation: (p1 + 1)(p2 + 1)(p3 + 1).
    """
    rows = []
    for z in range(zmax + 1):
        A, S = [], []
        for p in product(range(z + 1), repeat=3):
            if (p[0] - z) % 2:
                continue
            pf = [F(v) for v in p]
            inA = cones.in_A(A1_3, pf, z)
            inS = cones.in_S(A1_3, pf, z, 2) and (sum(p) - z) % 2 == 0
            if inA:
                A.append(p)
            if inS:
                S.append(p)
        dimsum = lambda L: sum(prod(v + 1 for v in x) for x in L)
        rows.append({"z": z, "labels_A": len(A), "labels_S": len(S), "dim_A": dimsum(A), "dim_S": dimsum(S),
                     "S_subset_A": set(S) <= set(A), "A_minus_S": [list(x) for x in sorted(set(A) - set(S))]})
    return rows


def calibration_table():
    rows = {}
    for comps in ([("A", 2)], [("A", 1)] * 3, [("A", 1)] * 4, [("B", 2)], [("G", 2)]):
        D = rd.data(rd.cartan_of(comps))
        rows[rd.name(comps)] = {"kappa": str(D["kappa"]), "H": D["H"], "required_d": str(1 / (D["kappa"] - 1)),
                                "B": [[str(v) for v in row] for row in cones.barycentric(D)]}
    return rows


def slices():
    out = {}
    for nm, D, d in (("A2,d=1", A2, 1), ("A2,d=3", A2, 3), ("A1^3,d=2", A1_3, 2), ("A1^4,d=1", A1_4, 1)):
        s = cones.slice_intervals(D, d)
        out[nm] = {k: [str(v) for v in val] for k, val in s.items()}
    return out


def gl2_pgl2_pgl2_reverse(zmax=8):
    """Identity (2) for the same pair: automorphic side on the dual group GL2 x SL2 x SL2, spectral side k[X] for
    G = GL2 x PGL2 x PGL2. Labels mu in X_*(T_dual) = X^*(T_G): coordinates (p1, p2, p3, r) with p_i = <mu, alpha_i^vee>
    and r = <mu, nu_G>; the lattice X^*(T_G) has p1 = r mod 2 and p2, p3 even. Spectral integrality:
    <mu, rho_G - theta_vee_G> = (p1 + p2 + p3 - r)/2."""
    rows = []
    for z in range(zmax + 1):
        A, S = [], []
        for p in product(range(z + 1), repeat=3):
            if (p[0] - z) % 2 or p[1] % 2 or p[2] % 2:
                continue
            pf = [F(v) for v in p]
            if cones.in_A(A1_3, pf, z):
                A.append(p)
            if cones.in_S(A1_3, pf, z, 2) and (sum(p) - z) % 2 == 0:
                S.append(p)
        rows.append({"z": z, "labels_A": len(A), "labels_S": len(S), "equal": set(A) == set(S)})
    return rows


def _series(num, den_power, nmax):
    """Coefficients of num(t) / (1 - t)^den_power up to t^nmax, num given as {power: coefficient}."""
    from math import comb
    return [sum(c * comb(n - k + den_power - 1, den_power - 1) for k, c in num.items() if k <= n) for n in range(nmax + 1)]


def hilbert_checks(nmax=10):
    """Graded dimensions of k[X] predicted by the corrected Proposition 5.1, against known Hilbert series.

    PGL3 x Gm (X = nilpotent cone of sl3, a complete intersection of the invariants of degrees 2 and 3):
      relevant (a, b) with a = b mod 3 contribute V_(a,b) in degrees n = a + b - s, s = 0..min(a, b).
    Garrett's G of (51) (X = quartic hypersurface in the 8-dimensional tensor space):
      relevant (p1, p2, p3) of equal parity contribute V_p in degrees n = p1 + p2 + p3 - 2s, s = 0..min(p).
    """
    nil = [0] * (nmax + 1)
    for a in range(nmax + 1):
        for b in range(nmax + 1):
            if (a - b) % 3:
                continue
            dim = (a + 1) * (b + 1) * (a + b + 2) // 2
            for s in range(min(a, b) + 1):
                n = a + b - s
                if n <= nmax:
                    nil[n] += dim
    hyp = [0] * (nmax + 1)
    for p in product(range(nmax + 1), repeat=3):
        if len({v % 2 for v in p}) > 1:
            continue
        dim = prod(v + 1 for v in p)
        for s in range(min(p) + 1):
            n = sum(p) - 2 * s
            if n <= nmax:
                hyp[n] += dim
    nil_ref = _series({0: 1, 2: -1, 3: -1, 5: 1}, 8, nmax)
    hyp_ref = _series({0: 1, 4: -1}, 8, nmax)
    return {"nilcone_sl3": nil, "nilcone_sl3_reference": nil_ref, "nilcone_match": nil == nil_ref,
            "hyperdeterminant": hyp, "hyperdeterminant_reference": hyp_ref, "hyperdeterminant_match": hyp == hyp_ref}
