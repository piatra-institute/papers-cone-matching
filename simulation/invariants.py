"""Direct computation of dim V_lambda^S, S = T'U', by exact linear algebra.

Two families are realized explicitly:
  * root type A1^k: V_lambda = Sym^{p_1} C^2 (x) ... (x) Sym^{p_k} C^2 twisted by a central character c = <lambda, theta_vee>;
    raising operators e_i = x_i d/dy_i; Lie(U') = span(e_1 - e_i); the T'-weight of a weight vector is
    sum_i mu_i / 2 + c, where mu_i is its SL2-weight.
  * root type A2: V_(a,b) = harmonic polynomials of bidegree (a, b) in x_1..x_3 and y_1..y_3 (kernel of
    sum_i d/dx_i d/dy_i), with E_ij = x_i d/dx_j - y_j d/dy_i; Lie(U') = span(E_12 - E_23, E_13); the T'-weight is
    mu_1 - mu_3 + c.
The computed dimension is compared with the criterion of Chen-Venkatesh Proposition 5.1,
0 <= <lambda, rho'> <= min_i <lambda, alpha_i^vee>, with and without the condition <lambda, rho'> in Z.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product

from functools import lru_cache

import sympy as sp


def _nullity(rows, ncols):
    if ncols == 0:
        return 0
    if not rows:
        return ncols
    return ncols - sp.Matrix(rows).rank()


def dim_invariants_A1(p, c):
    """dim of S-invariants in (x)_i Sym^{p_i} C^2 with central character c (T'-weight shift)."""
    k = len(p)
    # basis: a = (a_1..a_k), monomial prod x_i^{a_i} y_i^{p_i - a_i}; SL2 weight mu_i = 2 a_i - p_i
    zero = [a for a in product(*[range(pi + 1) for pi in p]) if sum(F(2 * ai - pi, 2) for ai, pi in zip(a, p)) + c == 0]
    if not zero:
        return 0
    if k == 1:
        return len(zero)  # U' is trivial
    # e_i maps monomial a to (p_i - a_i) * monomial a + e_i ; (e_1 - e_i) maps the zero-weight space into weight +1
    targets = {}
    rows = []
    for i in range(1, k):
        img = {}
        for col, a in enumerate(zero):
            for sgn, j in ((1, 0), (-1, i)):
                if a[j] < p[j]:
                    b = list(a)
                    b[j] += 1
                    b = tuple(b)
                    img.setdefault(b, {})
                    img[b][col] = img[b].get(col, 0) + sgn * (p[j] - a[j])
        for b, coeffs in img.items():
            rows.append([coeffs.get(col, 0) for col in range(len(zero))])
    return _nullity(rows, len(zero))


def criterion_A1(p, c, integral=True):
    s = sum(F(x, 2) for x in p) + c  # <lambda, rho'> with rho^vee = sum alpha_i^vee / 2
    ok = 0 <= s <= min(p)
    if integral:
        ok = ok and s.denominator == 1
    return ok


@lru_cache(maxsize=None)
def _harmonic_basis(a, b):
    X = sp.symbols("x1:4")
    Y = sp.symbols("y1:4")
    monos_x = [m for m in product(range(a + 1), repeat=3) if sum(m) == a]
    monos_y = [m for m in product(range(b + 1), repeat=3) if sum(m) == b]
    monos = [(mx, my) for mx in monos_x for my in monos_y]
    index = {m: i for i, m in enumerate(monos)}
    # Laplacian sum_i d/dx_i d/dy_i : bidegree (a, b) -> (a-1, b-1)
    if a == 0 or b == 0:
        return monos, index, sp.eye(len(monos))
    lower = {}
    for col, (mx, my) in enumerate(monos):
        for i in range(3):
            if mx[i] and my[i]:
                nx, ny = list(mx), list(my)
                nx[i] -= 1
                ny[i] -= 1
                key = (tuple(nx), tuple(ny))
                lower.setdefault(key, {})
                lower[key][col] = lower[key].get(col, 0) + mx[i] * my[i]
    rows = [[v.get(col, 0) for col in range(len(monos))] for v in lower.values()]
    ns = sp.Matrix(rows).nullspace()
    return monos, index, sp.Matrix.hstack(*ns)


def dim_invariants_A2(a, b, c):
    """dim of S-invariants in V_(a,b) for SL3 x Gm (theta_vee the central Gm), central character c."""
    monos, index, basis = _harmonic_basis(a, b)
    # T'-weight of a monomial: gl3 weight mu_k = deg_x_k - deg_y_k; weight = mu_1 - mu_3 + c
    wt = [mx[0] - my[0] - (mx[2] - my[2]) + c for mx, my in monos]

    def act(i, j, vec):
        """E_ij = x_i d/dx_j - y_j d/dy_i on a coefficient vector over the monomial basis."""
        out = {}
        for col, coef in enumerate(vec):
            if coef == 0:
                continue
            mx, my = monos[col]
            if mx[j]:
                nx = list(mx)
                nx[j] -= 1
                nx[i] += 1
                key = (tuple(nx), my)
                out[key] = out.get(key, 0) + coef * mx[j]
            if my[i]:
                ny = list(my)
                ny[i] -= 1
                ny[j] += 1
                key = (mx, tuple(ny))
                out[key] = out.get(key, 0) - coef * my[i]
        v = [0] * len(monos)
        for key, val in out.items():
            v[index[key]] += val
        return v

    # restrict to the T'-weight-zero part of V_(a,b): basis columns are weight vectors combinations, so project
    zero_cols = [k for k in range(len(monos)) if wt[k] == 0]
    if not zero_cols:
        return 0
    # vectors of V with weight zero: V is spanned by weight vectors; take V cap span(zero monomials)
    Bz = basis.extract(zero_cols, list(range(basis.shape[1])))
    others = [k for k in range(len(monos)) if wt[k] != 0]
    Bo = basis.extract(others, list(range(basis.shape[1]))) if others else sp.zeros(0, basis.shape[1])
    coeffs = Bo.nullspace() if others else [sp.eye(basis.shape[1])[:, i] for i in range(basis.shape[1])]
    if not coeffs:
        return 0
    Z = sp.Matrix.hstack(*[basis * v for v in coeffs])  # columns: weight-zero vectors of V
    rows = []
    for col in range(Z.shape[1]):
        v = list(Z[:, col])
        u1 = [p - q for p, q in zip(act(0, 1, v), act(1, 2, v))]
        u2 = act(0, 2, v)
        rows.append(u1 + u2)
    M = sp.Matrix(rows).T  # each column = image of a basis vector
    return Z.shape[1] - M.rank()


def criterion_A2(a, b, c):
    s = a + b + c  # <lambda, rho'> with rho^vee = diag(1, 0, -1) and theta_vee the central Gm
    return 0 <= s <= min(a, b)


def run_checks():
    out = {"A1": [], "A2": []}
    for k in (1, 2, 3):
        for p in product(range(4), repeat=k):
            if k == 3 and max(p) > 2:
                continue
            for c2 in range(-8, 9):  # c in (1/2) Z
                c = F(c2, 2)
                dim = dim_invariants_A1(list(p), c)
                out["A1"].append({"p": list(p), "c": str(c), "dim": dim,
                                  "criterion": criterion_A1(p, c), "criterion_without_integrality": criterion_A1(p, c, False)})
    for a in range(4):
        for b in range(4):
            for c in range(-7, 2):
                dim = dim_invariants_A2(a, b, c)
                out["A2"].append({"a": a, "b": b, "c": c, "dim": dim, "criterion": criterion_A2(a, b, c)})
    a1 = out["A1"]
    out["summary"] = {
        "A1_cases": len(a1),
        "A1_dims_at_most_one": all(r["dim"] <= 1 for r in a1),
        "A1_criterion_with_integrality_exact": all((r["dim"] == 1) == r["criterion"] for r in a1),
        "A1_cases_inequalities_hold_but_no_invariant": sum(r["criterion_without_integrality"] and r["dim"] == 0 for r in a1),
        "A1_cases_inequalities_hold": sum(r["criterion_without_integrality"] for r in a1),
        "A2_cases": len(out["A2"]),
        "A2_criterion_exact": all((r["dim"] == 1) == r["criterion"] and r["dim"] <= 1 for r in out["A2"]),
        "A2_cases_with_invariant": sum(r["dim"] == 1 for r in out["A2"]),
    }
    return out
