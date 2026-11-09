[English](../hpc.md) | [日本語](hpc.md) | [ドキュメント一覧](index.md)

# HPCとnative実行

CPU/CUDA、MPIの有無はビルド時に明示します。Make方式を維持し、新しいCMakeを要求しません。

```bash
make cpu
make cpu-mpi MPICXX=mpicxx
make cuda-single CUDA_SM=70 NVCC=nvcc
make cuda CUDA_SM=70 NVCC=nvcc MPICXX=mpicxx
```

Python APIの`build_native`/`pyhmd build`はbundled sourceから同じtargetを作ります。
`CXX`、`NVCC`、`MPICXX`、`CXXFLAGS`、`CUDAFLAGS`でtoolchainを明示できます。
build先はcheckoutでは`build/`、wheelではXDG cache下、または`--build-dir`/`PYHMD_BUILD_DIR`です。
自動コンパイル・SSH・qsub/sbatchは行いません。HPCでのビルド/計算はサイトが許可した割当内で行ってください。

## EngineからのMPI実行

```python
import os
from pyhmd import RunSpec

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
mpirun --hostfile "$PBS_NODEFILE" -np 4 build/pyhmd_cuda \
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
srun --ntasks=4 build/pyhmd_cuda \
  --config configs/thin_disk.ini --t-end 1 --output runs/slurm_native_example
```

[単一ノードPBS雛形](../../examples/pbs_gpu4.pbs.example)はqueue/group/moduleをサイトに合わせて設定します。
性能比較ではGPU型番・rank数・格子・受理step数・進んだ物理時刻・wall timeに加え、retry/補正量と出力設定も揃えて記録します。
