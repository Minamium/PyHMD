[English](../configuration.md) | [日本語](configuration.md) | [ドキュメント一覧](index.md)

# 設定

`Config`はproblemと数値モデル、`RuntimeConfig`は実行環境、`ResultOptions`は記録方針です。
すべてfrozen dataclassで、JSONの未知キー・重複キー・非有限値・boolを整数として渡す指定を拒否します。

```json
{
  "schema_version": 1,
  "config": {
    "problem": "thin_disk",
    "grid": {"nx1": 64, "nx2": 32},
    "time": {"t_end": 0.1}
  },
  "runtime": {"backend": "cpu", "mpi": false, "ranks": 1},
  "result": {
    "out_dir": "runs/disk",
    "output_interval": 0.05,
    "checkpoint_interval": 0.05
  }
}
```

pathは実行時のcwdを基準に解決します。specファイルの所在に暗黙にchdirしません。
`--out`は出力先を上書き指定します。launcherはJSON配列であり、shell文字列として評価しません。
INIの`#`/`;`/改行はnative parserの制限があるため、nativeへ渡すpath内では使用できません。

## 設定グループ

- `grid`: geometry、nx1/nx2、radial_grid、x1/x2のmin/max。
- `physics`: gamma、gm、eta_model/eta、h/corona、mu0/mu_ref/b0、bend、alpha_m/chi、cooling、initial_equilibrium、wave_amplitude。
- `numerics`: rho/p_floor、thermal_recovery、internal_energy_sync、dual_energy_eta、entropy閾値、CT補正、retry、cache、metric_ppm、ct_emf_scheme、cfl_staggered。
- `boundary`: lower1/upper1/lower2/upper2、inner_normal_field。
- `time`: cfl、t_end、max_steps。

選択肢は`pymhd.config`のdataclassが公開仕様です。
model=`resistive_mhd`、scheme=`ppm_hll_ssprk3`、dtype=`float64`のみを登録しています。
`thermal_recovery`など既存の比較optionを維持しています。HLLD、3D、FP32などの名前だけの追加はしていません。

辞書の部分指定はproblem presetへmergeします。dataclassを直接渡す場合は、そのdataclass全体を明示指定したものとして扱います。
設定の型と基本範囲はPython、組合せ制約は`Engine.plan()`のnative検証が担当します。
例えばCartesianのaxis境界、片側だけのperiodic、対応外の輸送診断は実行前に拒否します。

## Nativeとの互換性

既存INIはsectionごとの同名キーを区別しない、nativeと同じflat-key契約です。
未知キーや重複設定を拒否し、既存のentropy_switch/log_r等はcanonical fieldへ変換します。
`RunSpec.load`は出力間隔・checkpoint・restartも取り込みます。
`result.output_interval`と`log_interval`を省略するとproblem presetの値を使用します。

単純な数値integrity検査、輸送診断、科学的比較の設定を混同しません。
特にtransport_diagnosticsはsphere/log disk・face_hll・evolveなどのnative制約を持ちます。

設定の定義と選択肢は[config.py](../../src/pymhd/config.py)を参照してください。
輸送台帳などの保存は実行前に選びます。保存していない面流束を、後から実測の輸送診断として復元することはできません。
出力間隔はdtを切り詰めて数値軌道を変えるため、計算結果や性能の比較では固定します。

## Phase 0の計測設定

追加試験問題、`wave_mode`、円盤限定の`noise_amplitude`/`noise_seed`、`ResultOptions.midplane_diagnostics`は[Phase 0](phase0.md)を参照してください。赤道診断は既存の出力/ログ時刻にraw pressureの統計を保存し、追加のdt切詰めを行いません。
