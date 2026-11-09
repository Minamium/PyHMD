[English](../results.md) | [日本語](results.md) | [ドキュメント一覧](index.md)

# 結果と実行状態

| 保存物 | 内容 |
| --- | --- |
| `result.json` | schema、設定、runtime、native能力、binary SHA、状態、artifact参照 |
| `inputs/run.json` | 正規化された再実行用spec |
| `inputs/native.ini` | nativeへ実際に渡した全設定 |
| `inputs/pyhmd_*` | 実行に使ったbinaryの固定コピー |
| `run.log` | native/MPIのstdoutとstderr |
| `metadata.json`, `diagnostics.csv`, `status.json` | nativeの設定・物理量・到達/停止状態 |
| `snapshot_*.csv` | セルデータ、時刻・step |
| `checkpoint_step*/` | U/Ψ、clock、補正・面流束台帳。輸送診断時はsidecar |

Engineは実行前に出力先を確認し、`result.json`を排他的に作成して所有を確定します。
manifestは同じディレクトリ内のtemporary fileからrenameします。
新しい計算を既存の非空結果へ混ぜません。
キーボード中断/例外時は、そのEngineが起動したprocess群だけを終了処理してmanifestへ記録します。

| Native終了 | Result状態 | 意味 |
| --- | --- | --- |
| 0 + 正しいt_end到達 | completed | 積分完了。科学的成功とは独立 |
| 3 + max_steps | incomplete | step上限で途中終了 |
| 4 + persistent_retries | stopped | 深いretryの連続で意図的停止 |
| 他、またはstatus契約の不一致 | failed | SimulationError.resultから結果を読める |

`Result.load`は旧metadata/status形式も読みます。plottingや数値配列の読込は要求時だけimportします。
`Result.validate`は既存のstrict numerical integrity validatorで、floorがあれば不合格となります。
`scientific_validation`は自動で成功に変更しません。

## 再開

格子、rank分割、出力間隔、物理/数値設定、algorithm/source署名が一致する必要があります。
新しいout_dirを指定し、RuntimeConfig.restartへcheckpointディレクトリを渡します。
signatureに含まれないt_end/max_stepsなどは、保存時刻/stepより先に進める値へ変更できます。
main binary形式version1、旧magic、metadata version2を維持していますが、今回のsource構成変更でfingerprintは変わっています。
旧source世代のcheckpointを新版で黙って再開することはありません。
将来の互換変換は通常restartとは分けて設計・検証する必要があります。

## 実行中の状態とCLI

準備・実行・中断では`created`、`running`、`interrupted`になる場合もあります。
Engineの`run`は同期実行で、nativeの標準出力・標準エラーは`run.log`へ保存します。
legacy結果の検証では期待時刻を明示します。例：`pyhmd validate runs/legacy --expect-time 1`。
Python CLIには`--restart`オプションがないため、Pythonの`RuntimeConfig.restart`またはINIの`restart`を使います。
[READMEの再開例](../../README.ja.md)も参照してください。
