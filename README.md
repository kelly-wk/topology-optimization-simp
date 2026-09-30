# 2次元SIMPトポロジー最適化 / 2-D SIMP Topology Optimization

> **復旧監査で新規作成した再現実装です。Google Drive から回収した本人原作コードではありません。**
>
> **This is a clean-room reproduction created during a recovery audit, not recovered original source code.**

片持ち梁のコンプライアンス最小化を、依存を抑えた NumPy 実装で再現します。解析は4節点平面応力要素、SIMP材料補間、感度フィルタ、体積制約付き最適性基準法から構成されます。

This NumPy implementation reproduces a cantilever compliance-minimization benchmark using four-node plane-stress elements, SIMP material interpolation, a sensitivity filter, and a volume-constrained optimality-criteria update.

## 主な成果 / Deliverables

- `topology_simp.py` — deterministic solver and command-line interface
- `outputs/final_topology.svg` — optimized density field
- `outputs/final_density.csv` — numerical density matrix
- `outputs/history.csv` — convergence history
- `outputs/summary.json` — configuration and final metrics
- `tests/` — numerical and artifact checks
- `notebooks/topology_simp_colab.ipynb` — privacy-clean Colab launcher
- `docs/index.html` — minimal GitHub Pages project page
- `REPORT.md` and `PROVENANCE.md` — method, limits, and evidence boundary

## 実行 / Run

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
./run_example.sh
```

任意の解像度でも実行できます。

```bash
python3 topology_simp.py --nelx 24 --nely 8 --volfrac 0.4 --output outputs
```

## テスト / Test

```bash
python3 -m unittest discover -s tests -v
```

The tests verify matrix symmetry and semidefiniteness, volume conservation, objective improvement, density bounds, and artifact generation.

## ライセンス / License

このリポジトリに含まれるクリーンルーム実装と成果物は [MIT License](./LICENSE) で公開します。監査時に参照した第三者資料や PDF は含まれず、このライセンスの対象外です。詳しい範囲は [ライセンス範囲 / License Scope](./LICENSE_STATUS.md) を参照してください。

The clean-room implementation and artifacts included in this repository are released under the [MIT License](./LICENSE). Third-party audit materials and PDFs are not included and are outside this license. See [License Scope](./LICENSE_STATUS.md) for details.

## 監査上の注意 / Audit note

Drive の形状最適化関連フォルダには参考 PDF と研究計画メモがありましたが、実行可能な本人作成トポロジー最適化モデルはありませんでした。数理最適化アーカイブは講義資料のみです。別の最適化課題アーカイブは第三者に明示的に帰属するため、本実装には一切流用していません。詳細は `PROVENANCE.md` を参照してください。

The inspected Drive sources contained reference material but no verifiable user-authored topology-optimization implementation. Third-party coursework was quarantined and not reused. See `PROVENANCE.md`.
