---
name: pyhmd-develop
description: Implement or review PyHMD solver, CUDA/MPI, Python API, tests, packaging, or documentation. Use for code changes and regression diagnosis. Use pyhmd-use for operating existing simulations.
---

# PyHMD development

Read the checkout's [AGENTS.md](../../../AGENTS.md) and the contract affected by
the change. Python owns typed configuration, launch, provenance and results.
Native C++/CUDA owns state, MPI, RK integration, retries and accepted-step
accounting. Keep one canonical numerical implementation.

| Change | Entry point and verification |
| --- | --- |
| Configuration or API | Change `config.py`, INI/preset conversion and native validation together. Run `make test-python`. [Configuration](../../../docs/configuration.md) |
| Solver or physics | Inspect native `physics/`, `numerics/`, `geometry/`, `boundaries/` and `runtime/` as applicable. Run `make test` and the relevant analytic, conservation or refinement checks. [Validation](../../../docs/validation.md) |
| CUDA/MPI | Check `backend/`, `distributed/` and the failed stage log. Compile and test the affected backend, grid and rank count on allocated hardware. [HPC](../../../docs/hpc.md) |
| Outputs or restart | Preserve output ownership, signatures, schema compatibility, and exclusion of rejected trials from ledgers. [Results](../../../docs/results.md) |
| Bundled native sources or build | Build wheel/sdist and check an installed wheel's native build and run outside the checkout. [Architecture](../../../docs/architecture.md) |
| Docs or skills | Update paired English/Japanese docs and run `python3 scripts/check_docs.py`. For skill edits, validate frontmatter and review UI metadata. [Skills](../../../docs/skills.md) |

For refactors, compare full state, clocks and ledgers against a separately frozen
executable and settings. Match the output schedule because it can clip dt.
For deliberate numerical changes, use an independent physical or analytic
criterion and record expected differences. [Migration gate](../../../scripts/verify_migration.py)

When an operator comparison fails, isolate it from the rest of the suite.
Start both paths from identical input, compare local owned cells and shared CT
vertices with the correct global coordinates, and report the absolute difference,
component, location, and scale. Use an input/stencil scale when the reference
residual is near zero. Preserve the original failing result and validate any
revised criterion against a known discrepancy. [Phase 0 evidence](../../../docs/phase0.md)

For validation reports, identify the source generation, platform, grid, settings
and measured stage. Keep numerical changes, CPU checks, GPU build/run checks and
scientific conclusions identifiable. Commit compact evidence and exclude raw
runs, binaries and personal cluster details. [Experiment record](../../../docs/experiment-results.md)
