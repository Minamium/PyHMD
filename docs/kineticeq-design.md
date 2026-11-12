# Design adopted from kineticEQ

[English](kineticeq-design.md) | [日本語](ja/kineticeq-design.md) | [Index](index.md)

The reference is [kineticEQ](https://github.com/Minamium/kineticEQ), commit
`cb61018319d4543ca44bfb63004577d82523965a`. Its BGK equations, collision operator,
and implicit solvers are not reused as MHD algorithms.

| kineticEQ implementation / idea | PyMHD adaptation |
| --- | --- |
| `api/config.py`: frozen Config, explicit model/scheme/backend | Frozen groups, unsupported-name rejection, separate RuntimeConfig |
| `api/engine.py`: centralized orchestration | Python manages setup/results; the whole stepper remains native |
| `api/run_spec.py`: Python/JSON normalization | RunSpec normalizes JSON and existing INI through one path |
| `core/schemes/registry.py`: implemented factories | Only four real native targets in the registry |
| `distributed/config.py`: decomposition/rank/backend contracts | MPI/rank/launcher checks; multiple Python Engine ranks rejected |
| `distributed/context.py`: local rank versus GPU assignment | Local-rank mapping, allocation masks, UUID-duplication checks retained |
| `results/writer.py`, `results/io.py`: manifests and artifact references | Atomic result.json, relative artifacts, explicit failed/partial states |
| Offline rendering | Result.load, render, validate without reintegration |

Torch and all JIT machinery are not copied. C++/CUDA already owns MPI and time
integration, so a subprocess boundary fits the engine. Configuration/launch use
the standard library and visualization is optional. Numerical implementations
and tests keep one native owner instead of acquiring parallel Python copies.

The [two agent skills](skills.md) preserve that separation: one guides changes
and their validation, the other operates existing interfaces and saved results.
