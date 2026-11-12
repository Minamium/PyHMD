# Phase 0：計測・検証基盤

[English](../phase0.md) | [日本語](phase0.md) | [一覧](index.md)

[元の計画書](../plans/2026-11-09-high-order-grmhd-roadmap.md)は日本語のまま保存しています。
**G0は未達です。** 単一ノードV100×4の初回ジョブは対照・候補をビルドし、
最初のGPU等価性試験で停止しました。 [実験結果](experiment-results.md)

## 検証で修正した前提

- 現存CPUバイナリと設定を先に凍結しました。strict FPと未保有のプラットフォームの比較用バイナリは、
  checksum確認済みの`671fdf6`アーカイブから新しく構築する対照群です。過去の実在バイナリとは区別します。
- 非零の平均Byを持つ周期MPI境界では、Ψの準周期ジャンプ補正が必要です。
  独立した数値修正commitでwrap ghostを補正し、周期端の物理頂点を保持します。
  交換だけを行う平均磁場probeで旧2/4 rankの失敗と候補の成功を検証します。
- 既存GPUテストの大域nxによる局所配列走査を修正し、失敗をrank間で集約しました。
  このテスト自体はcache有無の比較であり、CPU/CUDA比較は別probeで行います。
- 点値・体積平均・面平均を混同せず、t=0の表現誤差、最終の体積平均誤差、CT表現誤差、
  sincで求めた面平均誤差を分けます。1/4周期も確認します。
- 線形波は基本量に固有ベクトル摂動を加えたため、保存量変換後にO(A²)の項を含みます。
  A≤1e−3に制限し、非線形誤差の頭打ちを記録します。
- ENTRは有限性検査と流束でも使います。未使用成分の省略や冗長通信削除は測定後まで保留します。
- GRの計量・逆計量・行列式・接続を19箇所に密保存すると、幾何だけで14,744 B/cell、192³で約104 GBです。
  1.5 KB/cellの見積もりとは両立せず、共有・対称性・再計算方針を定めて再見積もりが必要です。

AoSや2.5Dは設計上の制約であり、GRMHDそのものが不可能という意味ではありません。
軸対称のGR試験は有用ですが、継続的乱流/EHT比較という目標には後段の3D化が必要です。
後段の各手法は実装前に原典確認を続けます。

## 追加した機能と設定

`oblique_alfven`、4モードの`linear_wave`、`field_loop`、`rj2a`、`orszag_tang`、`blast_low_beta`を追加しました。
既存の問題・本番スキームは維持しています。`physics.wave_mode`は線形波のモード、
`noise_amplitude`/`noise_seed`は固定圧力での円盤初期密度摂動です。振幅0は旧経路を保ちます。

`result.midplane_diagnostics=true`は既存の受理出力・ログ時刻で`midplane_metrics.csv`を保存します。
追加のdt切詰めはせず、両CT方式で使えます。赤道faceのBθと**全エネルギー由来のraw pressure**からμを求め、
r=[1.1,1.5]のΨ差は端点を線形補間します。

`pymhd build --fp-mode strict`または`make cpu FP_MODE=strict`は別名の`_strict`実行ファイルを作ります。
CPUは`-ffp-contract=off`、CUDA device側は`-fmad=false`です。既定は変えず、
全てのcompiler/platform差を消すものでもありません。`--binary`で明示選択してください。
新版のstrict/default checkpoint署名は区別されます。

checkpoint読込検査はschema別に処理し、全台帳とchecksumを確認します。
未知のversionは拒否し、ファイルを書き換えたり、未定義のv2を推測で読んだりしません。

## 検証と実行

```bash
make test TEST_OUTPUT=runs/my_fresh_cpu_verification
make test-python
make test-convergence CONVERGENCE_OUTPUT=runs/my_fresh_orders
make benchmark
python3 scripts/verify_migration.py --reference "$FROZEN_ROOT/pyhmd_cpu" \
  --reference-configs "$FROZEN_ROOT/configs" --candidate build/pymhd_cpu \
  --report build/migration.json
# 許可されたGPU割当内のみ。freeze.json/source.tarも検証します。
make test-gpu FROZEN_SOURCE="$FROZEN_ROOT/source" GPU_OUTPUT=runs/my_fresh_gpu_audit
```

移行比較はghostを含むU/Ψ、補正・流束台帳、輸送sidecar、試行回数、時計、snapshot、
wall time以外の診断を対象にします。retryは実際に棄却を起こし、各バイナリが自分の世代のcheckpointを再開します。
U・熱台帳・流束台帳・counterだけを変えた故意変更の検出も試験します。

CPU/CUDAや分割間の共通入力には、ghostも含む**初期状態probe専用形式**を使います。
INI bytes、FP方針、格子、checksumを確認し、生成元も記録します。
受理済み時計や履歴を持たないため、通常の再開には使いません。
RHS差は微小な残差ではなく、入力規模×初期特性rateで規格化します。100step差と時間は別に報告します。

## V100ジョブと未完了事項

ローカルCPUの証拠は、[全状態の移行比較](../validation/phase0_migration_cpu.json)、
[斜め伝播波・線形波の実測次数](../validation/phase0_orders_cpu.json)、
[短時間ベンチマークとFP署名検査](../validation/phase0_smokes_cpu.json)、
[小格子のアンサンブル実行確認](../validation/phase0_noise_harness_cpu.json) に分けて記録しました。
最後の結果は600×128・t=200で予定しているノイズ幅の測定ではありません。
CPUの全native試験とPython契約試験9件も通過しました。

[バッチ雛形](../../scripts/phase0_gpu4.pbs.example)にサイト環境と凍結ソースの場所を外から渡します。
単一物理ノードとscheduler割当のV100×4を検証し、空いているGPUを勝手に選びません。
比較群は別ディレクトリでビルドし、積分前に実行ファイルのSHA256を固定します。

ジョブ245473は1200×256、4 rankで最初のキャッシュ有無のRHS比較に失敗しました。
失敗判定はログにあり、差の実測値と位置は未記録です。 [ジョブ記録](../validation/rokko_experiments_20261112.json)

後続の予定は周期平均磁場、CPU移行1/2/4 rank、共通入力CPU/CUDA・分割・世代比較、
GPU反復・再開・retry、斜め/線形波の収束、通常/strict × face_hll/corner_hllの4条件の測定です。
各条件は600×128・t=200、密度摂動1e−13を16本と無摂動対照1本です。
不正なbaseline/memberはノイズ幅を認定せず、タイムアウト時は次のGPUケースへ進みません。

最終時刻の分位は同じ物理時刻の測定値です。途中の共通時刻の分位は受理stepサンプルからの**補間**と明記します。
格子・スキーム・FP方針・摂動・積分区間に条件付けた感度であり、普遍的な誤差上限ではありません。

Brio–Wuの解の枝を含む独立した衝撃波参照、GPU結果、ノイズ幅の解釈、保留したWP0.5の最適化が残ります。
RJ2a/Brio–Wuの実行確認だけで波形精度を保証しません。GitHub CIはCPU・API・収束・文書・wheelで、CUDA実機は対象外です。
波の定義は[Felker & Stone 2018 §5.1–5.7](https://arxiv.org/html/1711.07439v2)で確認しました。
円偏光波はn=(1,2)/√5、e=(−2,1)/√5、速度+1で、速度摂動と磁場摂動は逆符号です。
将来のGR設計は現行ニュートンMHDとは別で、[White, Stone & Gammie 2016](https://arxiv.org/html/1511.00943v2)も参照します。
