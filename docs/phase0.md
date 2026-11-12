# Phase 0: measurement and verification

[English](phase0.md) | [日本語](ja/phase0.md) | [Index](index.md)

This work implements the measurement foundation of the
[original planning input](plans/2026-11-09-high-order-grmhd-roadmap.md), preserved
in Japanese. **G0 remains open.** The first single-node V100×4 job built the
controls and candidate, then stopped at its first GPU equivalence check.
[Experiment results](experiment-results.md)

## Review findings and adjustments

| Finding | Resolution |
| --- | --- |
| Existing binaries versus new strict-FP controls | Preserve the actual pre-change CPU executable and configs first. Construct new strict/platform controls only from the checksum-verified `671fdf6` source archive; label them as constructed controls, not historical binaries. |
| MPI periodic Ψ lacked its net-transverse-flux jump | A separate numerical-fix commit adjusts wrapped ghosts and preserves the physical periodic endpoints. An exchange-only mean-field probe must fail on the old 2/4-rank path and pass on the candidate. |
| GPU cache test indexed global nx into local buffers | Use local mesh sizes and a collective error check. That test is cached/uncached equivalence, not CPU/CUDA validation. |
| Point values, volume means, and face means differ | Report t=0 representation error, final analytic volume error, CT-representation error, and analytic sinc face error separately; include a quarter-period check. |
| Linear waves are not finite-amplitude exact solutions | Primitive eigenvector perturbations contain nonlinear O(A²) terms after conversion. Use A≤1e−3, sweep modes, and record the nonlinear floor scale. |
| Removing ENTR reconstruction is not automatically structural | Finite-state checks and entropy fluxes still use it. Defer removal and redundant-communication deletion until dedicated measured equivalence is available. |
| Future GR geometry storage contradicts 1.5 KB/cell | Literal dense storage of metric, inverse, determinant, and connections at 19 locations is 14,744 B/cell, about 104 GB at 192³ before state/workspace. Deduplicated/symmetric or recomputed layouts need a fresh estimate. This is not implemented here. |

AoS and 2.5D are engineering constraints, not proofs that GRMHD is impossible.
Axisymmetric GR tests remain useful; the proposed sustained-turbulence/EHT
comparison requires the later 3D work. Future phases still require original-paper
verification before implementation.

## Implemented interfaces

- Benchmarks: `oblique_alfven`, `linear_wave` (fast/Alfvén/slow/entropy),
  `field_loop`, `rj2a`, `orszag_tang`, `blast_low_beta`. Existing initial problems
  and default production schemes remain available.
- `physics.wave_mode` selects a linear mode. `noise_amplitude` and `noise_seed`
  perturb initial disk density at fixed pressure using a global-index hash;
  zero amplitude preserves the original initialization path.
- `result.midplane_diagnostics=true` writes `midplane_metrics.csv` at existing
  accepted-output/log samples, without additional dt clipping. It supports both
  CT variants and records authoritative equatorial Bθ, **raw total-energy
  pressure**, μ, and a linearly interpolated Ψ difference across r=[1.1,1.5].
- `pymhd build --fp-mode strict` / `make cpu FP_MODE=strict` build separate
  `_strict` executables. CPU uses `-ffp-contract=off`, CUDA device code uses
  `-fmad=false`. Defaults are unchanged; this does not remove all compiler or
  platform variability. Select strict binaries explicitly with `--binary`.
  New strict/default native checkpoint signatures are distinct.
- The checkpoint inspector dispatches by schema, verifies checksums, reads every
  state/ledger/sidecar and rejects unknown versions. It does not rewrite files or
  implement a hypothetical v2 layout.

## Test commands

```bash
make test TEST_OUTPUT=runs/my_fresh_cpu_verification
make test-python
make test-convergence CONVERGENCE_OUTPUT=runs/my_fresh_orders
make benchmark
python3 scripts/verify_migration.py --reference "$FROZEN_ROOT/pyhmd_cpu" \
  --reference-configs "$FROZEN_ROOT/configs" --candidate build/pymhd_cpu \
  --report build/migration.json
# Only in an authorized GPU allocation, with a verified freeze.json/source.tar:
make test-gpu FROZEN_SOURCE="$FROZEN_ROOT/source" GPU_OUTPUT=runs/my_fresh_gpu_audit
```

The migration gate compares full U/Ψ including ghosts, all correction/flux ledgers,
transport sidecars, counters, clocks, snapshots, and diagnostics except wall time.
Retry drivers must actually reject a trial. Each executable restarts its own
checkpoint generation; cross-generation signature checks are never bypassed.
Fault injection tests change U, a thermal ledger, a flux ledger, or counters and
confirm detection, including valid checksums on deliberately changed test files.

A separate **initial-state probe format** shares complete frozen U/Ψ including
ghosts between CPU/CUDA and rank decompositions. It checks identical INI bytes,
FP policy, dimensions, and checksum and records producer source identity. It has
no accepted clock or history and is never a normal restart path. RHS differences
are scaled by input magnitude times the initial characteristic rate, not a nearly
zero residual; 100-step differences and timings are reported separately.

## V100 job and remaining gates

Local CPU evidence is recorded separately: [full-state migration](validation/phase0_migration_cpu.json),
[measured oblique/linear orders](validation/phase0_orders_cpu.json),
[short benchmark and FP-signature checks](validation/phase0_smokes_cpu.json), and
[small-grid ensemble-harness smoke](validation/phase0_noise_harness_cpu.json).
The latter is not the requested 600×128, t=200 noise-floor measurement.
The full CPU native suite and nine Python contract tests also passed.

The [batch template](../scripts/phase0_gpu4.pbs.example) needs a site environment
and frozen-source directory supplied externally. It verifies one physical node
and exactly four scheduler-assigned V100 devices; it does not choose free GPUs.
The suite builds controls/candidates in separate directories and freezes SHA256s
before integration.

Job 245473 stopped at the first cached/uncached RHS check with the 1200×256
grid and four ranks. The test logged the failed assertion and omitted the
measured difference and location. [Job record](validation/rokko_experiments_20261112.json)

After that check, the planned sequence covers periodic mean fields, CPU
migration at 1/2/4 ranks, common-input CPU/CUDA and decomposition/generation
comparisons, GPU repeat/restart/forced retry, oblique/linear convergence,
and conditional noise ensembles. Each noise group uses 16 density perturbations of amplitude 1e−13 plus
an unperturbed control, grid 600×128, t=200, and fixed output scheduling. The four
groups are default/strict × face_hll/corner_hll. Invalid baselines/members cannot
qualify a noise floor. A timed-out case stops the suite before another GPU case.

Final-time quantiles are measured at the same physical time. Interior time-grid
quantiles are explicitly **interpolated** from native accepted-step samples.
This is empirical sensitivity conditional on the specified grid, scheme, profile,
perturbation, and integration interval, not an unconditional error bound.

Independent shock-profile references (including the Brio–Wu branch), complete
GPU results, quantitative noise interpretation, and deferred WP0.5 optimizations
remain required before declaring G0 complete. RJ2a/Brio–Wu execution checks do
not certify exact-profile accuracy. CUDA remains outside GitHub-hosted CI;
CPU shared code, API, convergence, documentation and wheel checks run there.

Wave and benchmark conventions were checked against
[Felker & Stone 2018, sections 5.1–5.7](https://arxiv.org/html/1711.07439v2).
The circular wave uses n=(1,2)/√5, transverse e=(−2,1)/√5 and propagation speed +1;
its velocity perturbation is opposite to its magnetic perturbation. The future
GR design is distinct from the current Newtonian implementation; compare
[White, Stone & Gammie 2016](https://arxiv.org/html/1511.00943v2).
