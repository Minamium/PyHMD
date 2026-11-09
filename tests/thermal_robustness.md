# Low-beta thermal regression

`thermal_robustness.cpp` exercises pressure and entropy in a genuinely
two-dimensional magnetic field. It complements the existing CT divergence,
weak-field diffusion and short disk smoke checks. None of these local tests
certifies a long-time disk calculation or measured CUDA performance.

Build and run from the repository root:

```sh
c++ -O3 -std=c++17 -Iinclude tests/thermal_robustness.cpp \
  -o build/thermal_robustness
build/thermal_robustness
```

The physical test has periodic x,y in [0,1], rho=1, p=1e-8,
vx=0.2, vy=0.13, and k=2pi. Its flux function and toroidal field are

```
Psi = A (cos(kx) + cos(ky)), A=1/k
Bx = dPsi/dy, By = -dPsi/dx, Bz=k Psi
```

Thus curl B=k B and J cross B=0. In ideal MHD this pattern translates
at the imposed constant velocity with constant gas pressure and density.
With constant eta it also decays by exp(-eta k²t). When all Ohmic heat
is immediately radiated, gas entropy is still constant. Initial magnetic
fields are averaged consistently with the flux-function representation;
Bz uses the analytic cell average. Finite-grid force imbalance is expected
to converge away, and is monitored through density and magnetic-mode errors.

The 32×32 and 64×64 cases advance to t=0.02 with ideal MHD and with eta=0.01 plus
Ohmic cooling. The production acceptance cases explicitly select auxiliary
internal-energy dual mode (`entropy_switch=3`) and CT energy correction. Acceptance requires raw positive pressure, no pressure or
density floor, specific entropy within 1% of its initial value, density
within 2% of rho=1, magnetic-mode translation/decay error below 0.3%, and
CT flux divergence at roundoff. These are stated accuracy tolerances, not
a claim that a dual energy or entropy switch conserves total energy exactly.
The periodic total energy change must equal the recorded entropy, auxiliary internal-energy, CT,
boundary and floor corrections minus the signed compatible cooling integral.
The absolute sum of the thermal corrections must decrease under refinement.
This tests an explicitly accounted modification of conservative energy;
it does not label the corrections physically exact.

An additional instantaneous test isolates eta-on minus eta-off thermal
RHS for Bz=sin(2pi x), v=0, rho=p=1. With full Ohmic cooling the continuum
resistive thermal derivative is zero. The finite-grid residual is reported
at 32 and 64 radial cells and must decrease under refinement or already be
below a roundoff-compatible absolute threshold. This test
does **not** accept a nonzero residual as exact discrete heat cancellation;
it distinguishes truncation residual from the pressure recovery applied
after a Runge–Kutta stage.

## Baseline that exposes the original failure

Before thermal recovery was added, the original 8-variable implementation
produced the following for the ideal 32×32 force-free case:

| Diagnostic | Value |
|---|---:|
| Steps to t=0.02 | 7 |
| Pressure/density floor events | 2960 |
| Minimum specific entropy / initial entropy | 9.97858e-5 |
| Maximum specific entropy / initial entropy | 24502.5 |
| Injected floor energy | 1.55874e-5 |
| Maximum density error | 0.00254021 |
| Relative magnetic-mode error | 3.82167e-5 |
| Relative CT flux divergence | 7.33791e-17 |

The original source SHA256 values were
`model.hpp=6e94ee4738640b50f6c338d6f6baad68f36d70d2f20d2564ac2945274172d6a1`
and `solver.hpp=d15d76f6737b7220d4c266e3a4a95181e693f5b74540c255aa93020d8f18a13c`.
The test therefore detects severe thermal damage even when the magnetic
mode and divergence look accurate.

The isolated Ohmic thermal residual had L1=0.00439273 at 32 cells and
0.00102831 at 64 cells, ratio 0.23409. It converged but was not identically
zero. A pressure floor or entropy recovery cannot by itself establish an
exact conservative exchange between magnetic and gas energies.

## Checks after recovery and compatible cooling

The updated implementation passed these local CPU tests:

| Case | Grid | Floors | Max specific-entropy error / initial | Total signed correction | Energy-budget defect |
|---|---|---:|---:|---:|---:|
| Ideal force-free translation | 32² | 0 | 2.94e-3 | −7.36e-5 | 3.12e-15 |
| Ideal force-free translation | 64² | 0 | 1.02e-3 | −3.85e-6 | 5.73e-16 |
| Resistive force-free translation, full cooling | 32² | 0 | 3.03e-3 | −7.60e-5 | 2.77e-15 |
| Resistive force-free translation, full cooling | 64² | 0 | 1.05e-3 | −4.21e-6 | 3.61e-15 |

The correction magnitudes fell by factors 0.0524 and 0.0554 under refinement.
The magnetic translation/decay mode error was at most 3.82e-5. The isolated
Ohmic thermal residual fell to L1≈3–4e-17 at both resolutions. These results
cover the stated smooth low-beta benchmarks, rather than a disk launch,
strong shock, or long-time accretion solution.

Additional new C++ tests use the same build pattern:

- `pressure_recovery.cpp`: 10,004 finite pressure reconstruction cases spanning
  disparate energy scales, plus invalid input rejection. This includes the
  cancellation counterexample that defeated the previous floor helper.
- `thermal_recovery.cpp`: smooth entropy selection, shock thermal-energy
  preservation, low-beta contact energy selection, nonoverlapping entropy/floor energy budgets, and rejection of
  invalid density, entropy, pressure and nonfinite candidates.
- `step_retry.cpp`: deliberately oversized double-rarefaction step, accepted
  reduced timestep, exact agreement with a fresh flat/HLL step, rejection
  accounting rollback, and exhausted-step state restoration plus failure data.
  Its exhausted-step subtest intentionally prints a failed-stage diagnostic.
- `contact_advection.cpp`: p=0.01, rho=1 versus 0.001, uniform velocity
  (0.2,0.13), weak constant Bphi=1e-4, periodic 128×4 domain, t=0.1.
  This is a pressure-balanced contact with unequal initial specific entropy,
  rather than the constant-entropy force-free wave above. Dual mode must retain
  pressure and velocity, conserve periodic mass/energy, use no entropy recovery,
  and require no floors or retries. A one-step full smooth-cell entropy
  diagnostic is reported separately and does not determine acceptance.
- `strong_contact.cpp`: Bphi=10 and 100 low-beta contacts, with dual and
  total-energy reference solutions through t=0.001, pressure/velocity accuracy,
  absence of floors/retries, and correction budgets normalized to gas energy.
- `adiabatic_flow.cpp`: low-beta expansion along a uniform magnetic field,
  testing rho=rho0/(1+a0t), p=p0/(1+a0t)^gamma and vx=a0*x/(1+a0t)
  in a central region remote from the boundaries. Expansion and compression
  instantaneous RHS both test uAux_dot=-gamma*uAux*a0. A compression trajectory
  is deliberately excluded because OUTFLOW ghosts prohibit its required inflow.

## Contact mixing check for the dual-mode choice

In the pressure-balanced contact test, dual mode completed 171 steps with
maximum relative pressure error 1.30e-14, velocity error 1.53e-15, no entropy
recovery events, no floors and no retries. Its relative energy change was
3.68e-15. Before the contact safeguard, the same initial state with full smooth-cell entropy selection gave
max(p)/p0=3.71394 after just one step. The conservative entropy tracer is
nonlinear in density at fixed pressure; numerically mixing its contact values
and then reconstructing pressure can add spurious heat. Restricting recovery
to small thermal-energy fractions avoids that error in this high-thermal-
fraction contact. Magnetically dominated contacts require the additional
safeguard tested below.

## Strong-field contact counterexample without the safeguard

`contact_mixing_diagnostic.cpp` is deliberately a diagnostic executable,
**not a passing accuracy test**. It repeats the pressure-balanced rho jump
with strong uniform Bphi. With no magnetic gradients, the CT energy correction
is zero and the exact pressure/velocity should remain constant. Total-energy
mode retains that state; the low-thermal-energy dual mode without the contact safeguard can still mix entropy
and overheat it:

| Bphi | beta | Mode | max(p)/p0 after 1 step | max(p)/p0 after 20 steps |
|---|---:|---|---:|---:|
| 10 | 2e-4 | Total energy | 1 (error <1e-12) | 1 (error <7e-12) |
| 10 | 2e-4 | Dual entropy | 4.15246 | 3.21143 |
| 100 | 2e-6 | Total energy | 1 (error <9e-11) | 1 (error <4e-10) |
| 100 | 2e-6 | Dual entropy | 4.15599 | 353.540 |

The Bphi=100 dual run reaches 20 steps at t=1.58104e-5, with velocity error
0.429 and recorded entropy-energy injection 0.01621. Initial total gas internal
energy is only 0.001875 in this domain, so this is about 8.65 times that value.
It triggers no pressure floors or rejected steps. A small correction relative
to the much larger magnetic energy would therefore conceal a serious thermal
error. Restricting entropy recovery by a thermal-energy ratio alone does not
resolve magnetically dominated contact mixing.

The unprotected counterexample remains reproducible with
`build/contact_mixing_diagnostic 0 --unguarded`; zero selects the 20-step
diagnostic. Supplying a positive first argument instead selects an end time.
The second option directly disables the contact sensor in this controlled
diagnostic, without modifying production source.

## Contact-safeguard acceptance

The added contact sensor protects pressure-balanced density jumps from entropy
recovery, including a buffer around detected cells. It uses total energy and
synchronizes the advected entropy there. This is a project-specific safeguard,
not a claim that these thresholds are PLUTO defaults. A negative/nonfinite
total-energy candidate in protected cells must still reject its stage.

`strong_contact.cpp` advances both total-energy reference and protected dual
cases to t=0.001. All cases retain positive raw pressure with zero floors and
zero rejected steps:

| Bphi | beta | Dual steps | Dual max pressure relative error | Total-E reference error | Dual velocity error | Entropy-energy correction |
|---|---:|---:|---:|---:|---:|---:|
| 10 | 2e-4 | 127 | 4.03e-11 | 4.12e-11 | 5.32e-13 | −8.12e-14 |
| 100 | 2e-6 | 1265 | 1.31e-8 | 6.15e-9 | 2.09e-12 | +1.24e-12 |

The Bphi=100 energy-budget defect is 8.20e-12, consistent with accumulated
double-precision total-energy rounding relative to magnetic energy. The test
separately limits pressure error to 5e-8, velocity error to 1e-9, and entropy
correction to 1e-8 times initial **gas** internal energy, so a large magnetic
energy cannot conceal the previous contact-heating failure. CT correction is
zero for the uniform toroidal field. The low-beta force-free tests and the
smooth/shock/contact recovery tests also pass with this safeguard enabled.

These cases establish behavior for this pressure-balanced contact family;
they do not certify all contacts, turbulence, shocks or the long-time disk.

## Auxiliary internal-energy mode acceptance

The current production acceptance uses mode 3: an auxiliary gas internal-energy
**density**, conservatively advected with flux uAux*v, with source -p div(v)
and the compatible Ohmic heating when it is not radiated. This avoids the
nonlinear entropy/density conversion responsible for the contact counterexample.
It is synchronized from total energy in trusted or shocked cells and supplies
pressure in the low-thermal-energy smooth cells. Entropy remains a separately
synchronized diagnostic in this mode. Both channel selection and all energy
corrections are explicitly tested.

Mode 3 passed the same strong contacts without relying on the entropy contact
safeguard:

| Bphi | Steps to t=0.001 | Max pressure relative error | Max velocity error | Internal-energy correction | Energy-budget defect |
|---|---:|---:|---:|---:|---:|
| 10 | 127 | 3.80e-14 | 5.59e-14 | −9.98e-14 | 9.71e-14 |
| 100 | 1265 | 2.18e-11 | 4.40e-12 | −4.14e-12 | 9.83e-12 |

The 2D force-free table above now reports mode 3. Its entropy error is below
0.31% at 32² and 0.11% at 64²; it is a convergent pressure-work truncation error,
rather than the nearly roundoff entropy preservation enforced by the entropy
equation. This distinction is retained in the report instead of reusing the
older mode-2 numbers.

For homologous expansion with a0=0.4 and beta=2e-8, the mode-3 trajectory through
t=0.05 gives density relative error 2.81e-11, pressure relative error 3.94e-9,
entropy relative error 3.99e-9 and velocity error 8.87e-12 in the assessed
central region, with zero floors or retries. Initial RHS errors for expansion
and compression are 4.00e-15 and 2.66e-15. These tests check the -p div(v)
coupling and adiabatic exponent independently of constant-pressure contact
advection. They do not claim a complete compression boundary test.
