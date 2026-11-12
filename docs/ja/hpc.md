[English](../hpc.md) | [日本語](hpc.md) | [ドキュメント一覧](index.md)

# HPCとnative実行

CPU/CUDA、MPIの有無はビルド時に明示します。Make方式を維持し、新しいCMakeを要求しません。

```bash
make cpu
make cpu-mpi MPICXX=mpicxx
make cuda-single CUDA_SM=70 NVCC=nvcc
make cuda CUDA_SM=70 NVCC=nvcc MPICXX=mpicxx
```

Python APIの`build_native`/`pymhd build`はbundled sourceから同じtargetを作ります。
`CXX`、`NVCC`、`MPICXX`、`CXXFLAGS`、`CUDAFLAGS`でtoolchainを明示できます。
build先はcheckoutでは`build/`、wheelではXDG cache下、または`--build-dir`/`PYMHD_BUILD_DIR`です。
自動コンパイル・SSH・qsub/sbatchは行いません。HPCでのビルド/計算はサイトが許可した割当内で行ってください。

## EngineからのMPI実行

```python
import os
from pymhd import RunSpec

spec = RunSpec.load("examples/thin_disk_cuda4.json")
spec = spec.with_runtime(launcher=("mpirun", "--hostfile", os.environ["PBS_NODEFILE"]))
result = spec.engine(out_dir="runs/my_allocated_job").run()
```

Engine自身を`mpirun -np 4 python ...`で多重起動しません。
native MPI executableが各rankのCUDA deviceとMPI lifetimeを所有します。
RuntimeConfigはGPU番号を上書きしません。schedulerが指定したCUDA_VISIBLE_DEVICESを保持します。
径方向の各rankに最低3セル必要です。

## Nativeの直接実行

Python orchestrationを使わず、既存手順と同じように実行できます。
この場合は`result.json`やbinaryの固定コピーは作られません。

```bash
mpirun --hostfile "$PBS_NODEFILE" -np 4 build/pymhd_cuda \
  --config configs/thin_disk.ini --t-end 1 \
  --checkpoint-interval 0.1 --output runs/native_example
```

`examples/pbs_gpu4.pbs.example`は単一ノードの雛形です。queue/group/moduleはサイトごとに設定してください。
MPI haloはhost stagingで、GPU-aware MPIや通信計算overlapを実装したとは扱いません。
今回の新構成はCPUで検証しており、CUDA/MPIの再ビルドと実測は別途必要です。

## MPI wrapperとSlurm

CUDA+MPI targetのflag取得はOpen MPI互換の`mpicxx --showme:compile`と`--showme:link`を前提とします。
別のwrapperを使う場合はMakefileの適合が必要です。SM指定例はV100=70、A30=80、GH200=90です。
Engineはlauncher配列に`-np N`を付けます。対応するMPI launcherを使い、`srun`は直接nativeを起動してください。

```bash
# サイトのMPI連携が有効な、確保済みSlurm割当内。
srun --ntasks=4 build/pymhd_cuda \
  --config configs/thin_disk.ini --t-end 1 --output runs/slurm_native_example
```

[単一ノードPBS雛形](../../examples/pbs_gpu4.pbs.example)はqueue/group/moduleをサイトに合わせて設定します。
性能比較ではGPU型番・rank数・格子・受理step数・進んだ物理時刻・wall timeに加え、retry/補正量と出力設定も揃えて記録します。

## Strict FPの対照実験

`pymhd build --backend cuda --mpi --cuda-sm 70 --fp-mode strict`で`pymhd_cuda_strict`を作り、`--binary`で明示選択します。既定のcompiler方針は変えません。凍結ソースの扱いとV100×4試験は[Phase 0](phase0.md)を参照してください。

## rokkoのPBSジョブ

rokko上のPyMHDチェックアウトから投入します。環境スクリプトは`python/3.12.1`と
`nvhpc/23.1`を読み込み、GCC 12.3、Open MPI、V100向けSM 70、ホスト側の単スレッド設定を
選びます。PBSが割り当てたGPUの可視性を保持します。 [環境スクリプト](../../scripts/rokko/env.sh)

割当内で手動ビルドするときはBashで`source scripts/rokko/env.sh`を実行してから、
`PYTHONPATH=src python3 -m pymhd build --backend cuda --mpi --cuda-sm 70`を使えます。

```bash
# V100×4の短時間円盤試走
qsub scripts/rokko/disk.pbs

# V100×1で自分のJSON/INI設定を使う
qsub -l select=1:ngpus=1:ncpus=8:mpiprocs=1:mem=64gb \
  -v PYMHD_SPEC=/path/to/run.json scripts/rokko/disk.pbs

# Phase 0。freeze.jsonとsource.tarに隣接する展開済み凍結ソースを渡す
qsub -v PYMHD_FROZEN_SOURCE=/path/to/frozen/source scripts/rokko/phase0_gpu4.pbs
```

[円盤ジョブ](../../scripts/rokko/disk.pbs)はPBSのノードファイルからMPI rank数を取得します。
GPUを2、4、8台にするときは単一ノードで`ngpus`と`mpiprocs`を揃えてください。
格子、終了時刻、checkpoint間隔は`PYMHD_SPEC`で指定します。既定は専用ディレクトリ内でビルドします。
`PYMHD_BINARY`でビルド済みCUDA+MPIバイナリ、`PYMHD_FP_MODE=strict`でstrictビルドを選べます。
確保済みの割当内ではチェックアウトから`bash scripts/rokko/run_disk.sh`でも実行できます。

PBSログには`runs/`または`PYMHD_RUN_BASE`配下のジョブディレクトリが表示されます。
円盤の結果はその中の`result/`、[Phase 0ジョブ](../../scripts/rokko/phase0_gpu4.pbs)の
進捗は`suite/progress.json`です。Phase 0の最初のGPU試験には未解決の差が残っています。
[実験結果](experiment-results.md)
