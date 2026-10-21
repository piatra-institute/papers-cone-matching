"""Reproduce every computation and figure: root data, cone comparison, redundancy classification,
invariant-vector checks, binary-code lattices and support comparisons.

    cd simulation && uv run python run_all.py
"""
from __future__ import annotations

import json
import platform
import time
from fractions import Fraction as F
from pathlib import Path

import matplotlib
import numpy
import scipy
import sympy

import cones
import examples as ex
import figures
import invariants as iv
import lattices
import rootdata as rd

OUT = Path(__file__).parent / "output"
MAX_RANK = 8


def box_compare_58(D, bound, zmax):
    agree = total = 0
    from itertools import product
    for x in product(range(bound + 1), repeat=D["rank"]):
        xf = [F(v) for v in x]
        for z in range(zmax + 1):
            total += 1
            agree += cones.in_A(D, xf, z) == cones.in_58(D, xf, z)
    return {"points": total, "agree": agree}


def main() -> None:
    t0 = time.time()
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    r = {"config": {"max_rank": MAX_RANK, "python": platform.python_version(), "numpy": numpy.__version__,
                    "scipy": scipy.__version__, "sympy": sympy.__version__, "matplotlib": matplotlib.__version__}}

    rows = cones.classify(MAX_RANK)
    red = [x["type"] for x in rows if x["redundant"]]
    match = [x["type"] for x in rows if x["redundant"] and x["required_d_integral"]]
    r["classification"] = {
        "types_with_kappa_gt_1": len(rows),
        "redundant_types": red,
        "convex_hull_types": [x["type"] for x in rows if x["convex_hull"]],
        "integral_calibration_types": [x["type"] for x in rows if x["required_d_integral"]],
        "cone_matching_types": match,
        "all_certificates_exact": all(x["certificates_exact"] for x in rows),
        "all_root_data_consistent": all(x["data_checks"] for x in rows),
        "rows": rows,
    }
    r["node_bound"] = cones.node_bound_irreducible()
    r["calibration"] = ex.calibration_table()
    r["slices"] = ex.slices()

    A2, A13, A14 = ex.A2, ex.A1_3, ex.A1_4
    B2 = rd.data(rd.cartan("B", 2))
    G2 = rd.data(rd.cartan("G", 2))
    r["box"] = {
        "A2_d1": cones.lattice_compare(A2, 1, 8, 16),
        "A1^3_d2": cones.lattice_compare(A13, 2, 5, 15),
        "A2_d3": cones.lattice_compare(A2, 3, 8, 16),
        "A1^4_d1": cones.lattice_compare(A14, 1, 3, 8),
        "B2_A_vs_58": box_compare_58(B2, 10, 30),
        "G2_A_vs_58": box_compare_58(G2, 10, 60),
    }
    r["examples"] = {"gl3": ex.gl3(), "a1_4_degree_one": ex.a1_4_degree_one(), "a1_4_witness": ex.a1_witness(),
                     "gl2_pgl2_pgl2": ex.gl2_pgl2_pgl2(8), "gl2_pgl2_pgl2_reverse": ex.gl2_pgl2_pgl2_reverse(8),
                     "hilbert": ex.hilbert_checks(12)}
    inv = iv.run_checks()
    r["invariants"] = {"summary": inv["summary"],
                       "sl2_gm_standard": {"dim": iv.dim_invariants_A1([1], F(0)), "inequalities": iv.criterion_A1([1], F(0), False),
                                           "rho_prime_pairing": "1/2"},
                       "gl2_sl2_sl2_label": {"dim": iv.dim_invariants_A1([1, 1, 2], F(-3, 2)),
                                             "inequalities": iv.criterion_A1([1, 1, 2], F(-3, 2), False)}}
    r["lattices"] = lattices.run()

    gl = r["examples"]["gl2_pgl2_pgl2"]
    first_gap = next(row["z"] for row in gl if row["labels_A"] != row["labels_S"])
    gap3 = next(row for row in gl if row["z"] == 3)
    r["gl2_pgl2_pgl2_summary"] = {"first_degree_with_gap": first_gap, "z8_labels_A": gl[8]["labels_A"], "z8_labels_S": gl[8]["labels_S"], "z3_labels_A": gap3["labels_A"], "z3_labels_S": gap3["labels_S"],
                                  "z3_dim_A": gap3["dim_A"], "z3_dim_S": gap3["dim_S"], "z3_missing": gap3["A_minus_S"]}

    checks = {
        "root_data_consistent": r["classification"]["all_root_data_consistent"],
        "redundancy_exactly_A2_B2_G2_A1^3": sorted(red) == sorted(["A2", "B2", "G2", "A1^3"]),
        "convex_hull_equals_redundancy": sorted(r["classification"]["convex_hull_types"]) == sorted(red),
        "certificates_exact": r["classification"]["all_certificates_exact"],
        "cone_matching_exactly_A2_A1^3": sorted(match) == sorted(["A2", "A1^3"]),
        "integral_calibration_exactly_A2_A1^3_A1^4": sorted(r["classification"]["integral_calibration_types"]) == sorted(["A2", "A1^3", "A1^4"]),
        "node_bound_excludes_all_other_irreducible": all((v["node_bound_passes"]) == (k in ("A1", "A2", "B2", "G2")) for k, v in r["node_bound"].items()),
        "closed_forms_node1": all(v["q_node1_closed_form_ok"] in (True, None) and v["kappa_closed_form_ok"] in (True, None) for v in r["node_bound"].values()),
        "box_A2_d1_agree": r["box"]["A2_d1"]["agree"] == r["box"]["A2_d1"]["points"],
        "box_A1^3_d2_agree": r["box"]["A1^3_d2"]["agree"] == r["box"]["A1^3_d2"]["points"],
        "box_A2_d3_disagree": r["box"]["A2_d3"]["first_disagreement"] is not None,
        "box_A1^4_d1_disagree": r["box"]["A1^4_d1"]["first_disagreement"] is not None,
        "box_B2_G2_reduce_to_58": all(r["box"][k]["agree"] == r["box"][k]["points"] for k in ("B2_A_vs_58", "G2_A_vs_58")),
        "gl3_separation_both_directions": r["examples"]["gl3"]["(1, 0, 0)"]["in_A"] and not r["examples"]["gl3"]["(1, 0, 0)"]["in_S_d3"]
                                           and r["examples"]["gl3"]["(2, 1, 0)"]["in_S_d3"] and not r["examples"]["gl3"]["(2, 1, 0)"]["in_A"],
        "a1_4_52_vs_40": (r["examples"]["a1_4_degree_one"]["dim_S"], r["examples"]["a1_4_degree_one"]["dim_A"]) == (52, 40),
        "invariants_A1_match_corrected_criterion": inv["summary"]["A1_criterion_with_integrality_exact"],
        "invariants_A1_uncorrected_criterion_fails": inv["summary"]["A1_cases_inequalities_hold_but_no_invariant"] > 0,
        "invariants_A2_match_criterion": inv["summary"]["A2_criterion_exact"],
        "f2_subspaces_67": r["lattices"]["subspaces_F2_4"] == 67,
        "five_codes_pass_both": r["lattices"]["passing_both"] == 5 and r["lattices"]["passing_equals_interval"],
        "garrett_pair_passes": all(r["lattices"]["named"][k]["rho_prime_G_integral"] and r["lattices"]["named"][k]["rho_prime_dual_integral"]
                                   for k in ("Garrett G (CV 51)", "Garrett dual (CV 53)")),
        "gl2_pgl2_pgl2_one_sided": r["lattices"]["named"]["GL2 x PGL2 x PGL2"]["rho_prime_G_integral"]
                                   and not r["lattices"]["named"]["GL2 x PGL2 x PGL2"]["rho_prime_dual_integral"],
        "gl2_pgl2_pgl2_spectral_subset": all(row["S_subset_A"] for row in gl),
        "hilbert_nilcone_sl3": r["examples"]["hilbert"]["nilcone_match"],
        "hilbert_hyperdeterminant": r["examples"]["hilbert"]["hyperdeterminant_match"],
        "gl2_pgl2_pgl2_reverse_identity_supports_equal": all(row["equal"] for row in r["examples"]["gl2_pgl2_pgl2_reverse"]),
    }
    r["checks"] = checks
    r["n_checks"] = len(checks)
    r["execution"] = {"status": "passed" if all(checks.values()) else "failed", "runtime_s": round(time.time() - t0, 2)}

    figures.make_all(r, OUT / "figures")
    (OUT / "results.json").write_text(json.dumps(r, indent=1, default=str))
    print(json.dumps({"checks": checks, "status": r["execution"]["status"]}, indent=1))


if __name__ == "__main__":
    main()
