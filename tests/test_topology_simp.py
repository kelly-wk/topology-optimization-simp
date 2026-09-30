import tempfile
import unittest
from pathlib import Path
import sys
import json

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topology_simp import (
    Config,
    element_dofs,
    element_stiffness,
    finite_element_analysis,
    optimize,
    save_result,
)


class TopologySimpTests(unittest.TestCase):
    def test_element_stiffness_is_symmetric_semidefinite(self) -> None:
        matrix = element_stiffness()
        self.assertEqual(matrix.shape, (8, 8))
        self.assertTrue(np.allclose(matrix, matrix.T, atol=1.0e-12))
        self.assertGreaterEqual(float(np.linalg.eigvalsh(matrix).min()), -1.0e-10)

    def test_small_problem_respects_volume_and_improves_compliance(self) -> None:
        config = Config(nelx=8, nely=4, max_iter=15, tolerance=0.03)
        result = optimize(config, verbose=False)
        self.assertAlmostEqual(float(result.density.mean()), config.volfrac, delta=2.0e-3)
        self.assertLess(result.final_compliance, result.initial_compliance)
        self.assertAlmostEqual(result.history[-1].compliance, result.final_compliance, places=10)
        self.assertTrue(np.all(result.density >= 1.0e-3))
        self.assertTrue(np.all(result.density <= 1.0))

    def test_final_compliance_matches_direct_analysis_and_summary(self) -> None:
        config = Config(nelx=7, nely=3, max_iter=8, tolerance=0.05)
        result = optimize(config, verbose=False)
        _, _, direct = finite_element_analysis(
            result.density,
            config,
            element_dofs(config.nelx, config.nely),
            element_stiffness(config.poisson),
        )
        self.assertAlmostEqual(result.final_compliance, direct, places=10)
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)
            save_result(target, config, result)
            saved_density = np.loadtxt(target / "final_density.csv", delimiter=",")
            _, _, saved_direct = finite_element_analysis(
                saved_density,
                config,
                element_dofs(config.nelx, config.nely),
                element_stiffness(config.poisson),
            )
            summary = json.loads((target / "summary.json").read_text(encoding="utf-8"))
            self.assertAlmostEqual(summary["final_compliance"], saved_direct, places=9)

    def test_outputs_are_written(self) -> None:
        config = Config(nelx=6, nely=3, max_iter=5, tolerance=0.10)
        result = optimize(config, verbose=False)
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)
            save_result(target, config, result)
            for name in ("history.csv", "final_density.csv", "final_topology.svg", "summary.json"):
                self.assertTrue((target / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
