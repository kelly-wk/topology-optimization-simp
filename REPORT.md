# 再現レポート / Reproduction Report

## 目的 / Objective

設計領域を有限要素に分割し、材料体積率を一定に保ちながら外力仕事（コンプライアンス）を小さくする密度分布を求めます。対象は左端固定、右端中央に鉛直集中荷重を受ける2次元片持ち梁です。

The benchmark minimizes compliance for a two-dimensional cantilever with a fixed left edge and a unit downward load at the midpoint of the right edge, subject to a material-volume constraint.

## 方法 / Method

1. 各要素を4節点平面応力要素として組み立てます。
2. 要素剛性は `E(x)=E_min+x^p(E_0-E_min)` で補間します。
3. 線形方程式 `K(x)u=f` を解き、要素ひずみエネルギーから目的関数と感度を求めます。
4. メッシュ依存性を抑えるため円錐重み付き感度フィルタを適用します。
5. 最適性基準法で密度を更新し、平均密度を指定体積率に合わせます。

The implementation assembles Q4 elements, applies SIMP interpolation, solves the linear elastic equilibrium equation, filters sensitivities, and performs an optimality-criteria update that enforces the prescribed mean density.

## 再現条件 / Reproduction configuration

The checked example uses Python 3.9.6 with NumPy 2.0.2, an 18 × 6 mesh, volume fraction 0.40, penalization 3.0, filter radius 1.5 element widths, and at most 80 iterations. Exact generated metrics are recorded in `outputs/summary.json`; convergence values are in `outputs/history.csv`.

The canonical run converged after 64 updates: compliance decreased from 1802.5563 to 563.5914 (68.7338%), final mean density was 0.3999934, and final maximum density change was 0.009516 below the 0.01 tolerance. Each history row now records compliance, volume, and change for the same post-update density state; the summary final compliance is also verified by a direct finite-element re-analysis of `final_density.csv`.

## 検証 / Verification

- The element stiffness matrix is symmetric and positive semidefinite within numerical tolerance.
- The design remains in `[0.001, 1.0]`.
- Mean density remains within `2×10⁻³` of the requested volume fraction in the test problem.
- Compliance decreases between the initial and final tested designs.
- All four expected output artifacts are created.

## 解釈上の限界 / Limitations

- This is a compact educational benchmark, not a validated production CAE solver.
- It uses dense linear algebra, so large meshes are intentionally out of scope.
- Point loading and a structured rectangular mesh are simplified benchmark conditions.
- No claim is made that this exact solver, parameter set, or result existed in the user's original Drive work.

## 結論 / Conclusion

監査対象から原作コードは回収できませんでしたが、トポロジー最適化の主要な数理要素を明示した、テスト可能かつ再実行可能な基準実装を新規に作成しました。成果を公開する場合も「復元した原作」ではなく「監査後に作成した再現実装」と表示する必要があります。

No original source was recovered. The result is therefore a tested, reproducible clean-room benchmark and must be presented publicly as a new reconstruction rather than recovered prior work.
