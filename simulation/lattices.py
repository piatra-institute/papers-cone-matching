"""Global forms of root type A1^3 with primitive central pairing d = 2, encoded by binary codes.

A character lambda of such a group G has integral coordinates (p_1, p_2, p_3, r) with p_i = <lambda, alpha_i^vee> and
r = <lambda, nu_G>. The roots are 2e_i and the primitive character trivial on the derived group is theta_G = 2e_4, so
2Z^4 <= X^*(T) <= Z^4 and the lattice is determined by the code C = X^*(T) / 2Z^4 in F_2^4. The primitivity of
theta_G and of nu_G requires e_4 not in C and some element of C with fourth coordinate 1.

<lambda, rho^vee + theta^vee> = (p_1 + p_2 + p_3 + r)/2, so rho'_G is a cocharacter iff C <= u^perp, u = (1,1,1,1).
The dual group has code C^perp (mod-2 dot product), so rho'_dual is a cocharacter iff u in C.
"""
from __future__ import annotations

from itertools import product

VEC = [v for v in product((0, 1), repeat=4)]
U = (1, 1, 1, 1)
E4 = (0, 0, 0, 1)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b)) % 2


def add(a, b):
    return tuple((x + y) % 2 for x, y in zip(a, b))


def span(gens):
    s = {(0, 0, 0, 0)}
    for g in gens:
        s |= {add(x, g) for x in s}
    return frozenset(s)


def all_subspaces():
    subs = set()
    for k in range(5):
        for gens in product(VEC, repeat=k):
            subs.add(span(gens))
    return sorted(subs, key=lambda s: (len(s), sorted(s)))


def perp(C):
    return frozenset(v for v in VEC if all(dot(v, c) == 0 for c in C))


def dim(C):
    return {1: 0, 2: 1, 4: 2, 8: 3, 16: 4}[len(C)]


def admissible(C):
    return E4 not in C and any(c[3] == 1 for c in C)


def integral_G(C):
    return all(dot(c, U) == 0 for c in C)


def integral_dual(C):
    return U in C


# Named groups, with their codes derived from explicit characters (comments give the derivation).
NAMED = {
    # Garrett's G (Chen-Venkatesh (51)): the tensor representation z^{-1} C^2 (x) C^2 (x) C^2 has p = (1,1,1), r = 1.
    "Garrett G (CV 51)": span([U]),
    # its dual (CV 53): characters (a_i, b_i) with p_i = a_i - b_i and r = sum(a_i + b_i)/3, so r = p_1 + p_2 + p_3 mod 2.
    "Garrett dual (CV 53)": span([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]),
    # SO4 x GL2 with SO4 = (SL2 x SL2)/mu2: p_1 = p_2 mod 2 and p_3 = r mod 2.
    "SO4 x GL2": span([(1, 1, 0, 0), (0, 0, 1, 1)]),
    # GL2 x PGL2 x PGL2: p_1 = r mod 2 on the GL2 factor, p_2 and p_3 even.
    "GL2 x PGL2 x PGL2": span([(1, 0, 0, 1)]),
    # GL2 x SL2 x SL2: p_1 = r mod 2, p_2 and p_3 arbitrary.
    "GL2 x SL2 x SL2": span([(1, 0, 0, 1), (0, 1, 0, 0), (0, 0, 1, 0)]),
}


def run():
    subs = all_subspaces()
    adm = [C for C in subs if admissible(C)]
    passing = [C for C in adm if integral_G(C) and integral_dual(C)]
    one_sided = [C for C in adm if integral_G(C) != integral_dual(C)]
    named = {}
    for nm, C in NAMED.items():
        named[nm] = {"dim": dim(C), "admissible": admissible(C), "rho_prime_G_integral": integral_G(C),
                     "rho_prime_dual_integral": integral_dual(C), "dual_code_matches": None}
    named["Garrett G (CV 51)"]["dual_code_matches"] = perp(NAMED["Garrett G (CV 51)"]) == NAMED["Garrett dual (CV 53)"]
    named["SO4 x GL2"]["dual_code_matches"] = perp(NAMED["SO4 x GL2"]) == NAMED["SO4 x GL2"]
    named["GL2 x PGL2 x PGL2"]["dual_code_matches"] = perp(NAMED["GL2 x PGL2 x PGL2"]) == NAMED["GL2 x SL2 x SL2"]
    return {
        "subspaces_F2_4": len(subs),
        "admissible": len(adm),
        "admissible_dual_closed": all(admissible(perp(C)) for C in adm),
        "passing_both": len(passing),
        "passing_dims": sorted(dim(C) for C in passing),
        "passing_self_dual": sum(perp(C) == C for C in passing),
        "passing_equals_interval": set(passing) == {C for C in subs if U in C and C <= perp(span([U]))},
        "one_sided": len(one_sided),
        "neither": len(adm) - len(passing) - len(one_sided),
        "named": named,
    }
