[English](../validation.md) | [日本語](validation.md) | [ドキュメント一覧](index.md)

# 検証

この構造整理では数値式を維持し、数値の移行検証とPython/native境界の検証を分けています。
既存コードとの一致は、元の物理モデルが論文と一致する証明にはなりません。

## CPUで確認した範囲

2026-11-09、macOS arm64、Apple Clang 21.0.0、Python 3.13.2で確認しました。

| 検査 | 確認結果 |
| --- | --- |
| `make test` | 既存のCPU回帰テスト一式に合格 |
| `make test-python` | 設定拒否、plan、結果所有、終了状態、正確な再開の7テストに合格 |
| 旧版との移行比較 | 全6ケースでCSV、保存変数U、磁束Ψ、checkpoint clockがbit一致 |
| パッケージ | wheel/sdist生成、checkout外のwheelからnativeビルド・計算・数値検査に合格 |
| 可視化 | 保存結果からoverview/history PNGとVTK生成、PNGを目視確認 |

既存のnative suiteは一様状態、Alfvén波、抵抗拡散、Brio–Wu、CT不変量、熱回復、拒否ステップ、
metric再構成、球面抵抗演算子、境界、checkpoint/再開、輸送台帳を含みます。
極端な低β衝撃波などに残る精度上の制限を解消したという主張ではありません。

移行比較には一様状態、Alfvén波、抵抗拡散、Brio–Wu、64×32円盤、600×128円盤を使いました。
600×128は100受理ステップ、到達時刻0.1053304632541587までです。
このケースは意図したmax_steps打切りで、長時間積分の検証ではありません。
詳細は[migration_cpu.json](../validation/migration_cpu.json)を参照してください。

参照する凍結版のsource指紋は[reference_import.json](../reference_import.json)に記録しています。
bit一致の適用範囲は同じcompiler/platformでの比較です。異なるGPU、compiler、MPI分割に拡張しません。

```bash
make test
make test-python
python3 scripts/verify_migration.py \
  --reference /path/to/frozen/pymhd_cpu \
  --candidate build/pyhmd_cpu --target-grid \
  --report build/migration_validation.json
```

移行比較には別途保存した旧版binaryが必要です。通常のCIはリポジトリ内で完結するCPUテストと
installed-wheelの契約を実行し、旧版との比較を実施したようには報告しません。

## 今回未確認の範囲

この新しい配置でのCUDA/MPIビルド、GPU性能、複数GPU結果はローカル環境で再検証していません。
既存nativeの対応コードとビルドtargetを維持していますが、CPU検証の結果で代用しません。
SF2014の定量再現も未達で、[科学的な検証範囲](scientific-status.md)に既知の差を記録しています。

GitHub ActionsではCPU suiteとPython 3.10のinstalled-wheel検査を設定しています。
各commitの実行結果は[Actions](https://github.com/Minamium/PyHMD/actions)で確認できます。

## Linux CIと文書・skillの確認

commit `0141fdd`では、LinuxのCPU/APIテスト、円盤の実行・検証・描画、Python 3.10でのwheel検証も合格しました。
[記録されたActions実行](https://github.com/Minamium/PyHMD/actions/runs/42583012368)と
[各commitのActions](https://github.com/Minamium/PyHMD/actions)を参照してください。
scalar/batched PPMの比較は入力ステンシルの大きさに対するFP64丸め誤差を許容します。
cache/retry/fluxの厳密一致が指定された検査は維持しています。

```bash
python3 scripts/check_docs.py
```

この確認はローカルMarkdownリンク、英日ページの存在、言語切替リンクを検査します。
外部サイト、翻訳の正しさ、表示レイアウトまで保証するものではありません。
skill作成時は作成環境のskill-creator validatorでfrontmatterも検証します。
