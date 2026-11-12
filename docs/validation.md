# Validation

[English](validation.md) | [日本語](ja/validation.md) | [Index](index.md)

The structural migration preserves numerical expressions. Migration comparison
and Python/native boundary checks are separate from scientific validation:
agreement with the previous code does not establish agreement with the paper.

## Observed CPU checks

Local checks on 2026-11-09 used macOS arm64, Apple Clang 21.0.0, Python 3.13.2.

| Check | Observed result |
| --- | --- |
| `make test` | Existing CPU regression suite passed |
| `make test-python` | Seven tests passed: configuration rejection, plans, output ownership, terminal states, exact restart |
| Migration | Six cases matched CSV, conservative U, flux Ψ, and checkpoint clock bits |
| Packaging | Wheel/sdist built; installed wheel built native code, ran and validated outside checkout |
| Visualization | Overview/history PNG and VTK generated; PNG inspected visually |
| Linux CI | CPU/API suite, disk run/validation/render, Python 3.10 installed-wheel checks passed at `0141fdd` |

See [baseline Actions run 42583012368](https://github.com/Minamium/PyMHD/actions/runs/42583012368)
and [current Actions](https://github.com/Minamium/PyMHD/actions).

Native tests cover uniform states, Alfvén waves, diffusion, Brio–Wu, CT invariants,
thermal recovery, rejected steps, metric reconstruction, spherical resistive
operators, boundaries, checkpoints/restart, and transport ledgers. They do not
resolve known extreme-low-β shock accuracy limits.

Migration cases are uniform, Alfvén, diffusion, Brio–Wu, a 64×32 disk, and a
600×128 disk. The latter runs 100 accepted steps to `t=0.1053304632541587`,
intentionally stopping at `max_steps`; it is not a long-time check. Evidence:
[migration_cpu.json](validation/migration_cpu.json),
[reference_import.json](reference_import.json). Bit equality applies to the same
compiler/platform comparison, not different GPUs, compilers, or decompositions.

```bash
make test
make test-python
python3 scripts/verify_migration.py \
  --reference /path/to/frozen/pymhd_cpu \
  --reference-configs /path/to/frozen/configs \
  --candidate build/pymhd_cpu --target-grid \
  --report build/migration_validation.json
```

Migration needs a separately preserved old binary. Normal CI uses self-contained
CPU and installed-wheel checks and does not claim to rerun that external gate.
Scalar/batched PPM comparison allows FP64 roundoff relative to input stencil
scale; cache/retry/flux checks retain exact equality where specified in the tests.

## Documentation checks

```bash
python3 scripts/check_docs.py
```

This checks local Markdown targets, English/Japanese page presence, and language
navigation. It does not verify websites, translation quality, or rendered layout.
Skill frontmatter is also checked with the skill-creator validator when authored.

## Hardware and scientific limits

CUDA/MPI builds, GPU performance, and multi-GPU results have not been revalidated
in this layout. CPU checks do not substitute for hardware tests. Quantitative
SF2014 reproduction remains unverified; see [scientific status](scientific-status.md).

[Phase 0](phase0.md) records the expanded migration/observer/probe/convergence checks. CUDA remains outside hosted CI; real-device acceptance is pending the authorized GPU job and independent shock-profile references.
