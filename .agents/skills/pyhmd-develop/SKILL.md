---
name: pyhmd-develop
description: Develop or review PyHMD source, numerical and API contracts, regression checks, packaging, and bilingual documentation. Use for implementation changes; use pyhmd-use for operating existing simulations and saved results.
---

# PyHMD development

Work in the PyHMD checkout associated with the request. Read its
[AGENTS.md](../../../AGENTS.md) and preserve existing user changes. Keep the user's
response language; English is primary for project documentation.

## Route the change

Read only the contracts relevant to the task, then inspect the implementation.

| Change | Read first | Implementation entry |
| --- | --- | --- |
| Structure or ownership | [Architecture](../../../docs/architecture.md), [kineticEQ mapping](../../../docs/kineticeq-design.md) | `src/pyhmd/engine.py`, `src/pyhmd/native/include/pyhmd/` |
| Configuration or public API | [Configuration](../../../docs/configuration.md) | `config.py`, `spec.py`, native `config/native.hpp`, presets |
| Launch, failure, output, restart | [Results](../../../docs/results.md) | `engine.py`, `runtime/`, `results/`, native `runtime/` and `io/` |
| Physics, discretization, GPU kernels | [Scientific status](../../../docs/scientific-status.md), [Validation](../../../docs/validation.md) | Native `physics/`, `geometry/`, `numerics/`, `boundaries/` |
| Build, MPI, device mapping | [HPC](../../../docs/hpc.md) | `backends/`, bundled Makefile, native `distributed/` and `backend/` |
| Docs or skills | [Documentation index](../../../docs/index.md), [Skills](../../../docs/skills.md) | English primary pages and corresponding Japanese pages |

Paths in implementation entries are relative to `src/pyhmd/` unless marked native.

## Preserve the design

- Keep Python responsible for configuration, explicit build/launch, provenance,
  and results. Keep state, MPI/CUDA lifetimes, RK stages, retry/rollback, and
  accepted-step accounting in native code. Do not add per-stage Python transfers.
- Support the existing coordinates, model, scheme, FP64, and settings unless the
  task explicitly requests an extension. A registry entry must have a real
  implementation; do not add silent fallback behavior.
- Change Python schema, preset/INI conversion, native validation, and interface
  tests together when a setting's contract changes. Reject unsupported combinations.
- Keep one numerical definition. Legacy `pymhd` headers forward to canonical
  `pyhmd` headers. Bundled native sources must remain buildable from an installed
  wheel outside the checkout.
- Preserve exclusive output ownership, actual command/binary provenance, terminal
  states, rejected-trial exclusion from ledgers, and strict checkpoint signatures.
  A changed native fingerprint is not permission to bypass restart checks.
- Separate refactoring from numerical changes. For a structural change, preserve
  trajectories against a frozen reference. For an intentional numerical change,
  use an independent analytic/conservation/refinement criterion appropriate to it.

## Verify the affected contracts

| Impact | Checks |
| --- | --- |
| Documentation/skills only | `python3 scripts/check_docs.py`; validate changed skill frontmatter/UI metadata; verify changed examples when needed |
| Python API or configuration | `make test-python` plus the changed rejection, failure, output, or restart behavior |
| Native implementation | `make test` and the relevant analytic/invariant/transport checks; migration comparison for refactors |
| Bundled sources or packaging | Build wheel/sdist; install and build/run outside the checkout |
| CUDA/MPI behavior | Compile and test the requested backend on authorized hardware; CPU success is not hardware validation |

The [migration script](../../../scripts/verify_migration.py) requires a separately
preserved old binary. Do not fabricate a baseline by rebuilding both sides from
the changed tree. Compare identical grid, IC/BC, physics, numerics, and output
intervals; output scheduling clips dt and can change the trajectory.

Keep generated binaries, runs, private cluster details, and credentials out of
commits. Maintain English/Japanese documentation together; keep API names, schema
identifiers, commands, and evidence consistent across languages.

Report what changed, why, the checks actually executed, and material gaps.
Distinguish integration completion, numerical integrity, and scientific
reproduction. Scope cluster execution, job mutation, and repository publication
to the user's request and existing authorization; skill invocation adds no authority.
