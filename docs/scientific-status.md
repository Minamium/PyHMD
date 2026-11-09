# Numerical and scientific scope

[English](scientific-status.md) | [日本語](ja/scientific-status.md) | [Index](index.md)

PyHMD is a research reimplementation targeting SF2014. Structural migration does
not expand or replace existing physical settings, coordinates, or operators.

The implemented system is classical compressible resistive MHD, an ideal gas,
and spherical point-mass gravity. It uses PPM/HLL/SSP-RK3, flux-function Ψ CT in
axisymmetry, and explicit anisotropic Ohmic terms. Overall convergence order
needs independent evaluation including boundaries, resistivity, and thermal
recovery; PPM/RK3 names alone do not establish third-order convergence.

## Known fidelity gaps

- Original paper-specific PLUTO 4.0 `init.c`, `definitions.h`, `pluto.ini`, and
  `res_eta.c` have not been obtained.
- The inner cell-centered magnetic angle and current Ψ boundary prescription
  are different discrete conditions.
- Ghost-velocity signs do not determine actual HLL mass-flux direction; compare
  the intended sink with its discrete behavior.
- Hydrostatic/magnetohydrostatic initial-equilibrium interpretations differ.
- Paper equivalence of compatible Ohmic cooling, CT correction, and internal-energy
  recovery is unconfirmed.
- Standard and strong diffusivity models are different experiments. Long runs
  do not trigger automatic model switching.
- CPU checks do not establish CUDA correctness; single-node MPI does not establish
  multi-node behavior.

Completion, small divB, and budgets closed with corrections do not alone establish
the physical solution. `completed` and validator success do not certify SF2014 reproduction.

References: [Stepanovs & Fendt 2014](https://arxiv.org/html/1407.3136),
[metric PPM geometry](https://arxiv.org/abs/1404.0537),
[UCT construction](https://doi.org/10.1016/j.jcp.2020.109748).

The [Phase 0 update](phase0.md) adds measurement problems and corrects wrapped periodic MPI flux handling. It does not replace the SF2014 physical boundary prescription or resolve its scientific reproduction gaps.
