#!/usr/bin/env python3
"""Small, dependency-light 2-D SIMP topology-optimization benchmark.

This is a clean-room reconstruction created during the recovery audit.  It is
not code recovered from the Google Drive archives.  The implementation uses a
Q4 plane-stress finite element, SIMP interpolation, a sensitivity filter, and
an optimality-criteria update for a cantilever compliance problem.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class Config:
    nelx: int = 18
    nely: int = 6
    volfrac: float = 0.40
    penal: float = 3.0
    rmin: float = 1.5
    max_iter: int = 80
    tolerance: float = 0.01
    young: float = 1.0
    young_min: float = 1.0e-9
    poisson: float = 0.30

    def validate(self) -> None:
        if self.nelx < 2 or self.nely < 2:
            raise ValueError("nelx and nely must both be at least 2")
        if not 0.0 < self.volfrac < 1.0:
            raise ValueError("volfrac must lie strictly between 0 and 1")
        if self.penal <= 1.0:
            raise ValueError("penal must be greater than 1")
        if self.rmin <= 0.0:
            raise ValueError("rmin must be positive")
        if self.max_iter < 1:
            raise ValueError("max_iter must be positive")
        if self.tolerance <= 0.0:
            raise ValueError("tolerance must be positive")
        if not 0.0 <= self.poisson < 0.5:
            raise ValueError("poisson must be in [0, 0.5)")


@dataclass
class Iteration:
    iteration: int
    compliance: float
    volume_fraction: float
    maximum_change: float


@dataclass
class Result:
    density: np.ndarray
    displacement: np.ndarray
    history: list[Iteration]
    converged: bool
    initial_compliance: float
    final_compliance: float


def element_stiffness(poisson: float = 0.30) -> np.ndarray:
    """Return the standard unit-modulus Q4 plane-stress stiffness matrix."""
    a11 = np.array(
        [[12, 3, -6, -3], [3, 12, 3, 0], [-6, 3, 12, -3], [-3, 0, -3, 12]],
        dtype=float,
    )
    a12 = np.array(
        [[-6, -3, 0, 3], [-3, -6, -3, -6], [0, -3, -6, 3], [3, -6, 3, -6]],
        dtype=float,
    )
    b11 = np.array(
        [[-4, 3, -2, 9], [3, -4, -9, 4], [-2, -9, -4, -3], [9, 4, -3, -4]],
        dtype=float,
    )
    b12 = np.array(
        [[2, -3, 4, -9], [-3, 2, 9, -2], [4, 9, 2, 3], [-9, -2, 3, 2]],
        dtype=float,
    )
    return (
        np.block([[a11, a12], [a12.T, a11]])
        + poisson * np.block([[b11, b12], [b12.T, b11]])
    ) / (24.0 * (1.0 - poisson**2))


def element_dofs(nelx: int, nely: int) -> np.ndarray:
    """Map each rectangular element to its eight global displacement DOFs."""
    rows: list[list[int]] = []
    for ex in range(nelx):
        for ey in range(nely):
            n0 = ex * (nely + 1) + ey
            n1 = (ex + 1) * (nely + 1) + ey
            n2 = (ex + 1) * (nely + 1) + ey + 1
            n3 = ex * (nely + 1) + ey + 1
            rows.append(
                [
                    2 * n0,
                    2 * n0 + 1,
                    2 * n1,
                    2 * n1 + 1,
                    2 * n2,
                    2 * n2 + 1,
                    2 * n3,
                    2 * n3 + 1,
                ]
            )
    return np.asarray(rows, dtype=int)


def cantilever_load_and_supports(nelx: int, nely: int) -> tuple[np.ndarray, np.ndarray]:
    """Clamp the left edge and apply a unit downward load at right mid-height."""
    ndof = 2 * (nelx + 1) * (nely + 1)
    force = np.zeros(ndof, dtype=float)
    load_node = nelx * (nely + 1) + nely // 2
    force[2 * load_node + 1] = -1.0
    fixed: list[int] = []
    for ey in range(nely + 1):
        fixed.extend((2 * ey, 2 * ey + 1))
    return force, np.asarray(fixed, dtype=int)


def filter_neighbours(nelx: int, nely: int, radius: float) -> list[tuple[np.ndarray, np.ndarray]]:
    """Precompute element indices and cone-filter weights."""
    neighbours: list[tuple[np.ndarray, np.ndarray]] = []
    reach = int(np.ceil(radius)) - 1
    for ex in range(nelx):
        for ey in range(nely):
            indices: list[int] = []
            weights: list[float] = []
            for nx in range(max(0, ex - reach), min(nelx, ex + reach + 1)):
                for ny in range(max(0, ey - reach), min(nely, ey + reach + 1)):
                    weight = radius - float(np.hypot(ex - nx, ey - ny))
                    if weight > 0.0:
                        indices.append(nx * nely + ny)
                        weights.append(weight)
            neighbours.append((np.asarray(indices, dtype=int), np.asarray(weights, dtype=float)))
    return neighbours


def finite_element_analysis(
    density: np.ndarray,
    config: Config,
    edofs: np.ndarray,
    ke: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Solve displacements and return per-element energy plus compliance."""
    ndof = 2 * (config.nelx + 1) * (config.nely + 1)
    stiffness = np.zeros((ndof, ndof), dtype=float)
    flat = density.T.reshape(-1)
    moduli = config.young_min + flat**config.penal * (config.young - config.young_min)
    for dofs, modulus in zip(edofs, moduli):
        stiffness[np.ix_(dofs, dofs)] += modulus * ke

    force, fixed = cantilever_load_and_supports(config.nelx, config.nely)
    free = np.setdiff1d(np.arange(ndof), fixed, assume_unique=True)
    displacement = np.zeros(ndof, dtype=float)
    displacement[free] = np.linalg.solve(stiffness[np.ix_(free, free)], force[free])

    energy = np.empty(len(edofs), dtype=float)
    for index, dofs in enumerate(edofs):
        local = displacement[dofs]
        energy[index] = float(local @ ke @ local)
    compliance = float(moduli @ energy)
    return displacement, energy, compliance


def filtered_sensitivity(
    density: np.ndarray,
    raw: np.ndarray,
    neighbours: list[tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    flat_x = density.T.reshape(-1)
    flat_raw = raw.T.reshape(-1)
    result = np.empty_like(flat_raw)
    for index, (indices, weights) in enumerate(neighbours):
        result[index] = float(weights @ (flat_x[indices] * flat_raw[indices])) / (
            max(1.0e-3, flat_x[index]) * float(weights.sum())
        )
    return result.reshape((density.shape[1], density.shape[0])).T


def optimality_criteria(density: np.ndarray, sensitivity: np.ndarray, target: float) -> np.ndarray:
    lower, upper = 0.0, 1.0e9
    move = 0.20
    safe = np.maximum(1.0e-30, -sensitivity)
    updated = density.copy()
    while (upper - lower) / (upper + lower + 1.0e-30) > 1.0e-4:
        midpoint = 0.5 * (lower + upper)
        candidate = density * np.sqrt(safe / midpoint)
        updated = np.maximum(
            1.0e-3,
            np.maximum(density - move, np.minimum(1.0, np.minimum(density + move, candidate))),
        )
        if float(updated.mean()) > target:
            lower = midpoint
        else:
            upper = midpoint
    return updated


def optimize(config: Config, *, verbose: bool = True) -> Result:
    config.validate()
    density = np.full((config.nely, config.nelx), config.volfrac, dtype=float)
    edofs = element_dofs(config.nelx, config.nely)
    ke = element_stiffness(config.poisson)
    neighbours = filter_neighbours(config.nelx, config.nely, config.rmin)
    history: list[Iteration] = []
    converged = False
    displacement, energy, compliance = finite_element_analysis(density, config, edofs, ke)
    initial_compliance = compliance

    for iteration in range(1, config.max_iter + 1):
        raw_flat = (
            -config.penal
            * (config.young - config.young_min)
            * density.T.reshape(-1) ** (config.penal - 1.0)
            * energy
        )
        raw = raw_flat.reshape((config.nelx, config.nely)).T
        sensitivity = filtered_sensitivity(density, raw, neighbours)
        new_density = optimality_criteria(density, sensitivity, config.volfrac)
        change = float(np.max(np.abs(new_density - density)))
        density = new_density
        displacement, energy, compliance = finite_element_analysis(density, config, edofs, ke)
        entry = Iteration(iteration, compliance, float(density.mean()), change)
        history.append(entry)
        if verbose:
            print(
                f"iter={iteration:03d} compliance={compliance:.6f} "
                f"volume={density.mean():.4f} change={change:.4f}"
            )
        if iteration >= 5 and change <= config.tolerance:
            converged = True
            break

    return Result(
        density=density,
        displacement=displacement,
        history=history,
        converged=converged,
        initial_compliance=initial_compliance,
        final_compliance=compliance,
    )


def _write_history(path: Path, history: Iterable[Iteration]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["iteration", "compliance", "volume_fraction", "maximum_change"],
            lineterminator="\n",
        )
        writer.writeheader()
        for entry in history:
            writer.writerow(asdict(entry))


def _write_density(path: Path, density: np.ndarray) -> None:
    np.savetxt(path, density, delimiter=",", fmt="%.12f")


def _write_svg(path: Path, density: np.ndarray, cell: int = 18) -> None:
    nely, nelx = density.shape
    grid_width, grid_height = nelx * cell, nely * cell
    left, right, top, bottom = 48, 50, 44, 34
    width, height = left + grid_width + right, top + grid_height + bottom
    load_x = left + grid_width
    load_y = top + grid_height / 2.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-label="Optimized cantilever material density with fixed support and downward unit load">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="4" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 Z" fill="#b42318"/></marker></defs>',
        f'<text x="{left}" y="24" font-family="system-ui,sans-serif" font-size="13" fill="#52606d">Material density — black = solid</text>',
    ]
    for ey in range(nely):
        for ex in range(nelx):
            value = float(np.clip(density[ey, ex], 0.0, 1.0))
            gray = int(round(255.0 * (1.0 - value)))
            parts.append(
                f'<rect x="{left + ex * cell}" y="{top + (nely - ey - 1) * cell}" '
                f'width="{cell}" height="{cell}" fill="rgb({gray},{gray},{gray})"/>'
            )
    parts.append(f'<line x1="{left - 4}" y1="{top}" x2="{left - 4}" y2="{top + grid_height}" stroke="#315b7d" stroke-width="4"/>')
    for offset in range(0, grid_height + 1, 12):
        parts.append(
            f'<line x1="{left - 18}" y1="{top + offset + 8}" x2="{left - 4}" y2="{top + offset}" stroke="#315b7d" stroke-width="1.5"/>'
        )
    parts.extend(
        [
            f'<line x1="{load_x}" y1="{load_y - 34}" x2="{load_x}" y2="{load_y + 24}" stroke="#b42318" stroke-width="3" marker-end="url(#arrow)"/>',
            f'<text x="{left - 2}" y="{height - 10}" text-anchor="middle" font-family="system-ui,sans-serif" font-size="12" fill="#315b7d">Fixed</text>',
            f'<text x="{load_x + 8}" y="{load_y - 23}" font-family="system-ui,sans-serif" font-size="12" fill="#b42318">F = 1</text>',
        ]
    )
    parts.append('</svg>')
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def save_result(output_dir: Path, config: Config, result: Result) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_history(output_dir / "history.csv", result.history)
    _write_density(output_dir / "final_density.csv", result.density)
    _write_svg(output_dir / "final_topology.svg", result.density)
    summary = {
        "provenance": "clean-room reconstruction; no recovered user source code was available",
        "benchmark": "2-D cantilever compliance minimization",
        "method": "Q4 plane-stress FEM + SIMP + sensitivity filter + OC update",
        "config": asdict(config),
        "iterations": len(result.history),
        "converged": result.converged,
        "initial_compliance": result.initial_compliance,
        "final_compliance": result.final_compliance,
        "compliance_reduction_percent": 100.0
        * (result.initial_compliance - result.final_compliance)
        / result.initial_compliance,
        "final_volume_fraction": float(result.density.mean()),
        "final_maximum_change": result.history[-1].maximum_change,
        "maximum_displacement": float(np.max(np.abs(result.displacement))),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nelx", type=int, default=18)
    parser.add_argument("--nely", type=int, default=6)
    parser.add_argument("--volfrac", type=float, default=0.40)
    parser.add_argument("--penal", type=float, default=3.0)
    parser.add_argument("--rmin", type=float, default=1.5)
    parser.add_argument("--max-iter", type=int, default=80)
    parser.add_argument("--tolerance", type=float, default=0.01)
    parser.add_argument("--output", type=Path, default=Path("outputs"))
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = Config(
        nelx=args.nelx,
        nely=args.nely,
        volfrac=args.volfrac,
        penal=args.penal,
        rmin=args.rmin,
        max_iter=args.max_iter,
        tolerance=args.tolerance,
    )
    result = optimize(config, verbose=not args.quiet)
    save_result(args.output, config, result)
    print(f"saved reproducible outputs to {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
