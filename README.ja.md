# PyMHD

[English](README.md) | [日本語](README.ja.md)

Pythonから設定・実行・解析できる、C++/CUDAによる有限体積法の抵抗性MHDソフトウェアです。
現在の対象は二次元Cartesianと球面軸対称2.5D、FP64、PPM/HLL/SSP-RK3/CTです。

kineticEQの **Config / Engine / Result** と、モデル・実行環境・保存結果を分ける設計を参考にしています。
時間積分はネイティブ実行ファイル内で行い、PythonやPyTorchを各RK段に介在させません。
インポートや`plan`はビルド・ジョブ投入・GPU選択を行いません。

## インストールと最初の計算

CPU版にはPython 3.10以降、C++17コンパイラ、GNU Makeが必要です。

```bash
git clone https://github.com/Minamium/PyMHD.git
cd PyMHD
python3 -m pip install -e '.[viz]'
pymhd build --backend cpu
pymhd init --problem alfven --output run.json
pymhd plan run.json --out runs/alfven_demo
pymhd run run.json --out runs/alfven_demo
pymhd inspect runs/alfven_demo
pymhd validate runs/alfven_demo
pymhd render runs/alfven_demo
```

設定・実行管理と基本的な検証はPython標準ライブラリだけで動きます。
`[analysis]`はNumPy、`[viz]`はNumPy・Matplotlib・Pillowを追加します。
wheelにもネイティブソースとpresetを含むため、インストール後も明示的にビルドできます。

## Python API

```python
from pymhd import Config, Engine, ResultOptions

cfg = Config(
    problem="thin_disk",
    grid={"nx1": 64, "nx2": 32},
    time={"t_end": 0.1},
)
result = Engine(
    cfg,
    result_out_dir="runs/disk_demo",
    result_options=ResultOptions(output_interval=0.05, checkpoint_interval=0.05),
).run()

print(result.status, result.complete)
result.validate()  # 数値的整合性の検査。論文再現の判定ではありません。
result.render(xlim=(0, 10), ylim=(0, 10))
```

辞書による部分指定は問題presetへ適用され、凍結dataclassに正規化されます。
未知のキー・未対応のmodel/scheme/dtypeを拒否します。
設定の組合せは実行前にネイティブ側でも検証します。
既存INIも`RunSpec.load("configs/thin_disk.ini")`、または`pymhd run configs/thin_disk.ini --out ...`で使用できます。

## 現在の対応範囲

| 項目 | 対応 |
| --- | --- |
| 物理 | 古典圧縮性MHD、理想気体、点質量重力、定数/標準SF2014/strong抵抗、Ohmic熱除去 |
| 座標 | 二次元Cartesian、球面軸対称2.5D。軸と赤道の対応する向きには制約あり |
| 数値法 | metric PPM、HLL、可変dtの陽的SSP-RK3、二次元磁束関数CT |
| 安定化 | CTエネルギー補正、選択式の熱回復、拒否ステップのrollback/retry、補正台帳 |
| CPU | 単一プロセス、MPI半径分割 |
| GPU | CUDA単GPU、MPI+CUDAで1 rank/GPU、FP64 |
| 保存 | CSV、署名・checksum付きbinary checkpoint、設定・実行ファイルの記録 |
| 解析 | Result再読込、Matplotlib、GIF/MP4、VTK、実離散輸送の任意診断 |

三次元・相対論・AMR・HLLD・陰解法・GPU-aware MPIは未実装です。
対応していない機能の名前だけを登録してフォールバックすることはありません。
今回の構造整理は物理モデルや境界の修正ではなく、SF2014の定量再現も未達です。
[科学的な検証範囲と既知の差](docs/ja/scientific-status.md)を参照してください。

## CUDA・MPI

CUDAは`nvcc`、MPIは`mpicxx`/`mpirun`を使用します。HPCではサイトの割当・module規則に従ってください。

```bash
# V100=70、A30=80、GH200=90
pymhd build --backend cuda --cuda-sm 70
pymhd run examples/thin_disk_cuda.json --out runs/gpu_demo

# MPI版。Pythonの外側でmpirunを重ねず、RuntimeConfigがnative ranksを起動します。
pymhd build --backend cuda --mpi --cuda-sm 70
pymhd run examples/thin_disk_cuda4.json --out runs/gpu4_demo
```

`CUDA_VISIBLE_DEVICES`を保持し、各native rankが割当内のlocal rankに対応するGPUを使います。
MPIの分割は半径方向だけ、haloと中央面表はhost経由です。
MPI実行時のエラーはnative側で協調処理します。
[HPCの実行方法](docs/ja/hpc.md)を参照してください。

`make cpu`、`make cuda-single CUDA_SM=70`、`make cuda CUDA_SM=70`も利用できます。
実行ファイルは`build/pymhd_*`です。rokkoの環境とPBSジョブは[HPCガイド](docs/ja/hpc.md)にまとめています。

旧名`pyhmd`でeditable installしていた場合は、このcheckoutを再インストールし、
importとCLIに`pymhd`を使ってください。既存結果は引き続き読み込めます。
以前のcheckpointから再開するときは、ソースと再開署名が一致する凍結済み実行ファイルを選びます。

## 結果と再開

`result.json`は設定、source fingerprint、binary SHA256、実行コマンド、状態、artifactを記録します。
nativeの`metadata.json`、`diagnostics.csv`、`status.json`と輸送台帳も保持します。
実行したbinaryを結果の`inputs/`へ保存し、後の再ビルドと区別します。
既存の非空出力先を拒否し、manifestは原子的に更新します。

```python
from dataclasses import replace
from pymhd import Result, RunSpec

result = Result.load("runs/disk_demo")
checkpoint = result.out_dir / result.checkpoints[-1]["path"]
spec = RunSpec.load("examples/thin_disk_cpu.json")
spec = replace(spec, config=replace(spec.config, time=replace(spec.config.time, t_end=0.2)))
continued = spec.with_runtime(restart=str(checkpoint)).engine(out_dir="runs/continued").run()
```

max_steps打切りは`incomplete`/終了コード3、再試行停止は`stopped`/4、異常は`failed`/2です。
`complete`やvalidator合格は科学的再現の成功を意味しません。
旧PyMHDのCSV結果は読み込めますが、旧ビルドcheckpointの通常再開はsource署名の違いで拒否します。

## 検証・設計

```bash
make test
make test-python
python3 scripts/verify_migration.py --reference /path/to/baseline/pymhd_cpu \
  --reference-configs /path/to/frozen/configs \
  --candidate build/pymhd_cpu --target-grid --report build/migration.json
```

- [設計とモジュール境界](docs/ja/architecture.md)
- [kineticEQから採用した設計](docs/ja/kineticeq-design.md)
- [設定とJSON仕様](docs/ja/configuration.md)
- [結果・状態・再開の契約](docs/ja/results.md)
- [検証方法と実測範囲](docs/ja/validation.md)
- [GPU実験・t=500の実測・Phase 0の結果](docs/ja/experiment-results.md)


## 言語とagent skills

英語を主版とし、`README.md`と`docs/*.md`に配置します。日本語は`README.ja.md`と`docs/ja/*.md`です。
API名・設定キー・コマンド・機械可読ファイルは英語のまま、対応機能と制限事項を両版で揃えます。
[日本語ドキュメント一覧](docs/ja/index.md)と[英語ドキュメント一覧](docs/index.md)を入口にしてください。

- [`$pymhd-develop`](.agents/skills/pymhd-develop/SKILL.md)：実装変更、数値/APIの契約、検証、配布、文書。
- [`$pymhd-use`](.agents/skills/pymhd-use/SKILL.md)：設定、実行、再開、ジョブ状態、目標時刻までの所要時間、結果解析・可視化。

skillは`.agents/skills/`に同梱しています。[利用例](docs/ja/skills.md)を参照してください。
グローバルなskill配置や個人設定は変更しません。

`plan`は検証とコマンド表示だけで、出力先を作りません。`run`は同期実行し、native出力を`run.log`へ保存します。
CUDA+MPIのMake targetはOpen MPI互換の`mpicxx --showme`を使います。
Engineのlauncherには`-np N`が付くため、Slurmの`srun`は[直接native実行](docs/ja/hpc.md)で使います。

文書・skillだけを変更した場合の確認は`python3 scripts/check_docs.py`です。
数値本体/APIを変更した場合は、それぞれに対応するsolver/APIテストも実行します。

License: MIT. 元のPyMHDの式・対応範囲を維持した構造整理として公開しています。

Phase 0の試験問題、strict FP、受動的な赤道診断とGPU検証の未完了事項は[Phase 0](docs/ja/phase0.md)にまとめています。
