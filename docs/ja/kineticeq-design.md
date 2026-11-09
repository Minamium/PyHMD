[English](../kineticeq-design.md) | [日本語](kineticeq-design.md) | [ドキュメント一覧](index.md)

# kineticEQから採用した設計

参照したのは[kineticEQ](https://github.com/Minamium/kineticEQ)の
commit `cb61018319d4543ca44bfb63004577d82523965a`です。
BGKの方程式・collision・implicit solverはMHDに流用していません。

| kineticEQの実装・考え方 | PyHMDでの具体化 |
| --- | --- |
| `api/config.py`: frozen Config、明示的なmodel/scheme/backend | frozen Configと各設定群。未対応名を拒否。実行環境をRuntimeConfigへ分離 |
| `api/engine.py`: orchestrationをEngineへ集約 | 設定/結果管理はPython、stepper全体はnativeが所有 |
| `api/run_spec.py`: Python/JSONを同じ型に正規化 | RunSpecがJSONと既存INIを同じ経路へ変換 |
| `core/schemes/registry.py`: 実装のfactory表 | 実在する4native backendだけをregistryに登録 |
| `distributed/config.py`: 分割・rank・backendの契約 | MPI有無/rank/launcherを厳密に検査、Pythonの多重MPI起動を拒否 |
| `distributed/context.py`: local rankとGPU割当の区別 | nativeのlocal rank対応、割当mask保持、UUID重複検査を維持 |
| `results/writer.py`, `results/io.py`: manifestとartifact参照 | 原子的result.json、相対artifact参照、失敗/途中終了の明示 |
| 保存済み結果をofflineで描く | Result.load、render、validateが再積分せず動く |

kineticEQのTorch依存やJITのすべてをコピーすることは目的にしていません。
今回のC++/CUDAエンジンは単体でMPI・時間積分を完結できるため、subprocess境界を選びました。
その結果、CPUの設定・起動管理は標準ライブラリだけ、可視化はoptional dependencyにできます。
同じ手法の別実装を増やさず、既存の検証をnativeの責務のまま残しています。

[2つのagent skill](skills.md)もこの分離に合わせています。開発用は変更の責務と検証へ、
利用用は既存の設定・実行・保存結果の解析へ案内します。
