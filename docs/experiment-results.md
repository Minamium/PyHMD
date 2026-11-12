# Experiment results through 2026-11-12 JST

[English](experiment-results.md) | [日本語](ja/experiment-results.md) | [Index](index.md)

The SF2014 reference reproduction is unverified. The earlier solver reached
t=10,000 on V100×4 and showed persistent inner-disk magnetic flux loss and
net mass injection at the inner disk boundary. Phase 0 passed the local CPU
checks described below. Its first V100×4 verification stopped at the
cached/uncached RHS comparison. [Scientific scope](scientific-status.md)

## Measured time to t=500

Rokko job 245189 ran the spherical 2.5D standard SF2014 disk with CFL 0.4,
face-HLL CT, transport diagnostics, and four MPI ranks on one node with four
V100s. These are the exact accepted rows at t=500 and t=1000 from the completed
historical run. Wall time is the cumulative solver diagnostic and excludes
queue wait and batch setup. The source fingerprint starts with `a9eabd4f`.
[Measured rows, settings and checksums](validation/rokko_experiments_20261112.json)

| Grid | Steps at t=500 | Wall at t=500 | Steps at t=1000 | Wall at t=1000 |
| --- | ---: | ---: | ---: | ---: |
| 600×128 | 216,723 | 460.83 s (7 min 41 s) | 417,041 | 887.66 s |
| 1200×256 | 473,567 | 2,227.30 s (37 min 7 s) | 878,164 | 4,145.73 s |

The Phase 0 suite specifies t=200 per member. At t=500, multiplying the
historical 600×128 sample by 68 members gives about 8 h 42 min of solver time.
This projection uses one historical configuration for all members. CT choice,
FP mode, output cadence and adaptive dt can change the runtime.

## Historical GPU runs and performance

- On Miyabi GH200, the early 600×128 run intended for t=100 stopped after the
  last recorded sample near t=53.14. A saved-state investigation linked
  repeated pressure corrections to low-β thermal losses from the mismatch
  between the CT magnetic update and the cell energy/momentum update. The
  precise first failing RK stage remains unknown. A later GH200 1200×256 case
  failed near t=17.013.
- The revised solver added pressure recovery and retry, CT/energy accounting,
  face-aware CFL, and CUDA caching and reductions. On pomme A30, final
  generation `6d4132b4…` completed 1200×256 to t=30 and 600×128 to t=100,
  both with zero floors and rejected trials. A separate 600×128 continuation
  was stopped by request and retained a t=400 checkpoint, with later accepted
  progress logged through t=440. A matched 1200×256 A30 benchmark at block
  size 128 measured 14.8403 → 12.2547 ms per step, a 1.211× speedup.
- The table is the older, single-node Rokko strong-scaling test with a common
  binary and inputs. It reports median ms per accepted step in the timed
  window, including MPI, with initialization and endpoint output excluded.

| Grid | 1 V100 | 2 V100 | 4 V100 | 8 V100 |
| --- | ---: | ---: | ---: | ---: |
| 600×128 | 3.9855 | 2.8936 | 2.1776 | 1.8500 |
| 1200×256 | 13.3686 | 8.1468 | 5.0061 | 3.4915 |

The 4-to-8 GPU speedup was 1.18× and 1.43× respectively. The measurements cover
one node. Phase 0 solver performance timings remain pending.

## SF2014 disk and transport results

| Experiment | Observed result |
| --- | --- |
| 600×128 to t=10,000, V100×4 | 4,010,933 accepted steps, zero rejected steps, 2 h 25 min 44 s batch wall. The integration reached the target and the post-run zero-floor check failed. Inner equatorial magnetization μ fell from 0.01002 to 1.18×10⁻⁶. Accepted HLL flux at the inner disk boundary gave cumulative mass outflow of −1.1498 per radian, indicating net inflow into the domain. |
| 600×128 and 1200×256 to t=1000, V100×4 | Both runs completed with zero rejected steps. Mean equatorial μ over r≈1.1–1.5 fell to 7.98×10⁻⁶ and 9.26×10⁻⁶. In the annulus containing the inner boundary, the signed cumulative UCT jump contributed about −95% and −87% of the initial flux difference. The contribution was below 0.04% in magnitude when the inner boundary was excluded at r≈1.1. |
| Repetition, restart and transport diagnostics at t=100/200 | Same-binary repeats and physical U/Ψ/clock with diagnostics on/off matched at the checked outputs. The historical continuous/restart snapshots also matched. The historical old/new binary pair differed by about 4.29×10⁻¹³ in maximum absolute U at 1,000 steps and by about 2.06% in inner μ at t=200. The cause of the generation difference is open. |
| CFL 0.4 versus 0.2 at t=200, 600×128 | Accepted steps rose from 96,551 to 193,357. Inner μ was 2.1046×10⁻⁴ and 2.0872×10⁻⁴. Both retained an inner-boundary UCT contribution of about −82% of the initial annular flux difference. |

The UCT ledger records a signed update to the flux function. The individual
contributions and the net HLL mass flux motivate tests of the inner ghost field,
electric field, and sink prescription. The initial resistivity audit was
consistent with the standard model's selected coefficients and anisotropy.
Paper-specific PLUTO inputs and a like-for-like boundary comparison remain open.
[Model and fidelity gaps](scientific-status.md)

## Phase 0 checks and remaining work

- The local CPU suite and nine Python contract tests passed. The frozen CPU
  comparison matched full state, clocks and ledgers across eight cases,
  including 600×128 for 100 steps, restarts and a forced three-rejection
  case. The common-input CPU probe also matched RHS, EMF and the 100-step
  trajectory. [Migration](validation/phase0_migration_cpu.json) ·
  [Common-input probe](validation/phase0_common_input_cpu.json)
- All 42 CPU wave cases completed, including the quarter-period propagation
  check. For the oblique Alfvén wave at amplitude 0.1 with face-HLL CT,
  refinement from 32×16 to 64×32 gave observed L1/L∞ orders of 2.04/1.54.
  Short non-wave benchmarks, strict/default checkpoint isolation and a
  small-grid CPU ensemble smoke test also completed. The wheel's bundled native
  source built and ran outside the checkout. [Wave orders](validation/phase0_orders_cpu.json) ·
  [Smoke checks](validation/phase0_smokes_cpu.json) ·
  [CPU ensemble smoke](validation/phase0_noise_harness_cpu.json)
- Rokko job 245473 verified four V100s and passed four main build commands and
  25 probe compile commands for the frozen-source controls and candidate. It
  exited with status 1 after 16 min 17 s of batch wall. The first GPU check,
  `cached_uncached_cuda_rhs` at 1200×256 with four ranks, tripped its 2×10⁻¹²
  component-normalized threshold. The log retains the assertion. Error
  magnitude, component and location remain unknown. Subsequent GPU checks and
  the production-grid ensemble were unreached. [Job evidence](validation/rokko_experiments_20261112.json)

G0 remains open. The next hardware check needs the cache comparison isolated
with absolute and scaled RHS/EMF differences and their locations recorded.
Independent shock profiles, the production-grid noise measurement, and the
deferred component/communication changes also remain open.
[Phase 0 methods and scope](phase0.md)

The earlier detailed reports remain in the legacy PyMHD workspace. The records
used here are `gh200_failure_investigation_20261105.md`,
`thermal_repair_20261106.md`, `rokko_v100_scaling_1to8_20261106.md`,
`rokko_reference_analysis_20261108.md`, `transport_resolution_analysis_20261108.md`,
`parallel_reproducibility_cfl_analysis_20261109.md`, and
`model_consistency_followup_20261109.md`.
