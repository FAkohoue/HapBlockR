#!/usr/bin/env python3
"""Independent Clarabel and OSQP oracle for DesiredGainR validation.

Input and output are JSON files.  This script deliberately contains no
DesiredGainR transformations and receives only H, f, lower and upper.
"""

import argparse
import contextlib
import importlib.metadata
import io
import json
import math
from pathlib import Path

import clarabel
import numpy as np
import osqp
import scipy.sparse as sparse


def decode_bound(values):
    if isinstance(values, (str, int, float)):
        values = [values]
    decoded = []
    for value in values:
        if value == "Inf":
            decoded.append(np.inf)
        elif value == "-Inf":
            decoded.append(-np.inf)
        else:
            decoded.append(float(value))
    return np.asarray(decoded, dtype=float)


def objective(x, hessian, linear):
    return float(0.5 * x @ hessian @ x + linear @ x)


def solve_clarabel(
    hessian, linear, lower, upper, tolerance, reverse_constraints=False
):
    dimension = len(linear)
    identity = np.eye(dimension)
    rows = []
    bounds = []
    for index in range(dimension):
        if math.isfinite(upper[index]):
            rows.append(identity[index])
            bounds.append(upper[index])
        if math.isfinite(lower[index]):
            rows.append(-identity[index])
            bounds.append(-lower[index])
    constraint = sparse.csc_matrix(
        np.asarray(rows, dtype=float).reshape((-1, dimension))
    )
    rhs = np.asarray(bounds, dtype=float)
    if reverse_constraints:
        constraint = constraint[::-1, :]
        rhs = rhs[::-1]
    hessian_sparse = sparse.csc_matrix(np.triu(hessian))
    cones = [clarabel.NonnegativeConeT(len(bounds))]
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.tol_gap_abs = tolerance
    settings.tol_gap_rel = tolerance
    settings.tol_feas = tolerance
    settings.max_iter = 500
    solver = clarabel.DefaultSolver(
        hessian_sparse, linear, constraint, rhs, cones, settings
    )
    solution = solver.solve()
    x = np.asarray(solution.x, dtype=float)
    return {
        "status": str(solution.status),
        "solution": x.tolist(),
        "objective": objective(x, hessian, linear),
        "dual": np.asarray(solution.z, dtype=float).tolist(),
        "slack": np.asarray(solution.s, dtype=float).tolist(),
        "iterations": int(solution.iterations),
    }


def solve_osqp(hessian, linear, lower, upper, tolerance, polishing):
    dimension = len(linear)
    solver = osqp.OSQP()
    solver.setup(
        P=sparse.csc_matrix(np.triu(hessian)),
        q=linear,
        A=sparse.eye(dimension, format="csc"),
        l=lower,
        u=upper,
        eps_abs=tolerance,
        eps_rel=tolerance,
        max_iter=200000,
        polishing=polishing,
        verbose=False,
    )
    # OSQP 1.1 prints a polishing note even with verbose=False. Keep the JSON
    # oracle machine-readable without weakening or skipping the polished run.
    with contextlib.redirect_stdout(io.StringIO()):
        solution = solver.solve(raise_error=False)
    x = np.asarray(solution.x, dtype=float)
    return {
        "status": str(solution.info.status),
        "solution": x.tolist(),
        "objective": objective(x, hessian, linear),
        "dual": np.asarray(solution.y, dtype=float).tolist(),
        "iterations": int(solution.info.iter),
    }


def solve_problem(problem, tolerance):
    hessian = np.atleast_2d(np.asarray(problem["H"], dtype=float))
    linear = np.atleast_1d(np.asarray(problem["f"], dtype=float))
    try:
        lower = decode_bound(problem["lower"])
        upper = decode_bound(problem["upper"])
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"Invalid bound encoding in {problem['problem_id']}: "
            f"lower={problem['lower']}, upper={problem['upper']}"
        ) from error
    if hessian.shape != (len(linear), len(linear)):
        raise ValueError("H has incompatible dimensions")
    symmetry_scale = max(1.0, float(np.max(np.abs(hessian))))
    if np.max(np.abs(hessian - hessian.T)) > 1e-12 * symmetry_scale:
        raise ValueError("H is not symmetric")
    hessian = 0.5 * (hessian + hessian.T)
    if np.min(np.linalg.eigvalsh(hessian)) <= 0:
        raise ValueError("H is not strictly positive definite")
    return {
        "problem_id": problem["problem_id"],
        "clarabel": solve_clarabel(
            hessian, linear, lower, upper, tolerance
        ),
        "clarabel_reordered": solve_clarabel(
            hessian, linear, lower, upper, tolerance, True
        ),
        "osqp_unpolished": solve_osqp(
            hessian, linear, lower, upper, tolerance, False
        ),
        "osqp_polished": solve_osqp(
            hessian, linear, lower, upper, tolerance, True
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    results = [
        solve_problem(problem, args.tolerance)
        for problem in payload["problems"]
    ]
    args.output.write_text(
        json.dumps(
            {
                "versions": {
                    "clarabel": importlib.metadata.version("clarabel"),
                    "osqp": importlib.metadata.version("osqp"),
                    "numpy": importlib.metadata.version("numpy"),
                    "scipy": importlib.metadata.version("scipy"),
                },
                "results": results,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
