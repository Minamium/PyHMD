#!/usr/bin/env python3
"""Snapshot estimates of SF2014 Sections 3.3/4.1 and Appendix B fluxes.

Not the solver's discrete HLL flux ledger. See numerical_model.md for geometry,
normal conventions, moving-surface limitations and normalization.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from pyhmd.plotting._viz import load_snapshot, snapshot_files


def interp(x, xp, fp):
    """Linear interpolation with linear endpoint extrapolation."""
    x = np.asarray(x)
    i = np.clip(np.searchsorted(xp, x) - 1, 0, len(xp) - 2)
    return fp[i] + (x - xp[i]) / (xp[i + 1] - xp[i]) * (fp[i + 1] - fp[i])


def sample(a, r, theta, rc, tc):
    i = int(np.clip(np.searchsorted(rc, r) - 1, 0, len(rc) - 2))
    j = int(np.clip(np.searchsorted(tc, theta) - 1, 0, len(tc) - 2))
    x, y = (r-rc[i])/(rc[i+1]-rc[i]), (theta-tc[j])/(tc[j+1]-tc[j])
    return (1-y)*((1-x)*a[j, i]+x*a[j, i+1])+y*((1-x)*a[j+1, i]+x*a[j+1, i+1])


def fields(s):
    d, m = s.data, s.metadata
    if np.any(d['rho'] <= 0) or np.any(d['pressure'] <= 0):
        raise ValueError('Flux diagnostics require positive density and pressure')
    r, t = d['x1'], d['x2']
    rc, tc = r[0], t[:, 0]
    # Centered finite differences of cell fields, not the staggered CT curl.
    dr = lambda a: np.gradient(a, rc, axis=1, edge_order=2)
    dt = lambda a: np.gradient(a, tc, axis=0, edge_order=2)
    jr = dt(np.sin(t)*d['Bphi'])/(r*np.sin(t))
    jt = -dr(r*d['Bphi'])/r
    jp = (dr(r*d['Btheta'])-dt(d['Br']))/r
    # Midplane reflection implies Br=0. Face Btheta is absent from snapshots.
    rho, pressure = d['rho'][-1], d['pressure'][-1]
    bp = abs(d['Btheta'][-1])
    va, cs, mu = bp/np.sqrt(rho), np.sqrt(pressure/rho), bp**2/(2*pressure)
    h = 2*m['epsilon']*r*np.sin(t)
    z = r*np.cos(t)
    profile = np.exp(-2*np.maximum(z/h-1, 0)**2)
    kind = m['eta_model']
    if kind == 0:
        ep = np.full_like(r, m['eta'])
    elif kind == 1:
        ep = m['alpha_m']*np.maximum(interp(r*np.sin(t), rc, va), 0)*h*profile
    else:
        if 'mu_ref' not in m:
            raise ValueError('Strong diffusivity requires mu_ref in snapshot metadata')
        ref = m['mu_ref']
        ep = (m['alpha_m']*np.sqrt(2*ref)*(np.maximum(interp(r*np.sin(t), rc, mu), 0)/ref)**2
              *np.maximum(interp(r*np.sin(t), rc, cs), 0)*h*profile)
    et = ep if kind == 0 else m['chi']*ep
    v = np.stack([d[k] for k in ('vr', 'vtheta', 'vphi')], axis=-1)
    b = np.stack([d[k] for k in ('Br', 'Btheta', 'Bphi')], axis=-1)
    electric = -np.cross(v, b)+np.stack([et*jr, et*jt, ep*jp], axis=-1)
    poynting = np.cross(electric, b)
    speed2 = np.sum(v*v, axis=-1)
    cyl = r*np.sin(t)
    flux = {}
    for k, name in enumerate(('r', 'theta')):
        flux['mass_'+name] = d['rho']*v[..., k]
        flux['j_adv_'+name] = cyl*d['rho']*d['vphi']*v[..., k]
        flux['j_mag_'+name] = -cyl*d['Bphi']*b[..., k]
        flux['e_kin_'+name] = .5*speed2*d['rho']*v[..., k]
        flux['e_grav_'+name] = -m['gm']/r*d['rho']*v[..., k]
        flux['e_enthalpy_'+name] = m['gamma']/(m['gamma']-1)*d['pressure']*v[..., k]
        flux['e_poynting_'+name] = poynting[..., k]
    flux['density'] = d['rho']
    flux['angular_density'] = cyl*d['rho']*d['vphi']
    flux['energy_density'] = (d['pressure']/(m['gamma']-1)+.5*d['rho']*speed2
                              +.5*np.sum(b*b, axis=-1)-m['gm']/r*d['rho'])
    flux['cooling'] = (et*(jr*jr+jt*jt)+ep*jp*jp)*m['cool_ohmic']
    if kind == 0:
        eta_mid = np.full_like(rc, m['eta'])
    elif kind == 1:
        eta_mid = m['alpha_m']*va*(2*m['epsilon']*rc)
    else:
        eta_mid = m['alpha_m']*np.sqrt(2*ref)*(mu/ref)**2*cs*(2*m['epsilon']*rc)
    return flux, mu, eta_mid, rc, tc


TERMS = ('mass', 'j_adv', 'j_mag', 'e_kin', 'e_grav', 'e_enthalpy', 'e_poynting')


def surface(s, mode, angle):
    rc, tc = s.data['x1'][0], s.data['x2'][:, 0]
    if mode == 'fixed':
        return np.full_like(rc, angle)
    result = []
    for v in s.data['vr'].T:
        crossings = np.flatnonzero((v[:-1] >= 0) & (v[1:] < 0))
        if not len(crossings):
            result.append(np.nan)
        else:
            j = crossings[-1]  # first outward crossing when walking up from disk
            result.append(tc[j] - v[j]*(tc[j+1]-tc[j])/(v[j+1]-v[j]))
    return np.asarray(result)


def evaluate(s, args):
    m = s.metadata
    if m['geometry'] != 'spherical' or abs(m['x2_faces'][-1]-np.pi/2) > 1e-12:
        raise ValueError('Requires a spherical upper hemisphere ending at the equator')
    if min(s.shape) < 3:
        raise ValueError('At least three cells per direction are needed for current estimates')
    rf, tf = np.asarray(m['x1_faces']), np.asarray(m['x2_faces'])
    rlo, rhi = args.inner, args.outer
    if not rf[0] <= rlo < rhi <= rf[-1]:
        raise ValueError('Requested control radii must lie inside the mesh')
    angle = np.pi/2-2*m['epsilon'] if args.theta is None else args.theta
    if not tf[0] < angle < tf[-1]:
        raise ValueError('Disk surface angle must be inside the angular domain')
    f, mu, eta, rc, tc = fields(s)
    surf = surface(s, args.surface, angle)
    knots = np.unique(np.r_[rlo, rf[(rf > rlo) & (rf < rhi)], rc[(rc > rlo) & (rc < rhi)], rhi])
    theta = interp(knots, rc, surf)
    if not np.isfinite(theta).all():
        raise ValueError(f'{s.path.name}: vr-zero surface missing in control volume; use --surface fixed for early times')
    if np.any((theta <= tf[0]) | (theta >= tf[-1])):
        raise ValueError('Extracted surface leaves angular domain')
    factor = 4*np.pi if args.bipolar else 2*np.pi
    row = {'time': s.time, 'step': s.step, 'inner_radius': rlo, 'outer_radius': rhi,
           'azimuth_hemisphere_factor': factor, 'surface_theta_min': theta.min(), 'surface_theta_max': theta.max()}
    # Radial caps: exact solid-angle weights on cut cells; center evaluation.
    for radius, th0, label in ((rlo, theta[0], 'inner'), (rhi, theta[-1], 'outer')):
        for term in TERMS:
            total = 0.
            for a, b in zip(tf[:-1], tf[1:]):
                low = max(a, th0)
                if b > low:
                    total += sample(f[term+'_r'], radius, (low+b)/2, rc, tc)*radius**2*(np.cos(low)-np.cos(b))
            row[f'{term}_acc_{label}'] = -factor*total
    # Upper surface theta(r): outward vector area = (r² sinθ θ', -r sinθ) dr dphi.
    # Gauss integration includes the sloping surface's radial flux component.
    for term in TERMS:
        total = 0.
        for a, b, ta, tb in zip(knots[:-1], knots[1:], theta[:-1], theta[1:]):
            slope = (tb-ta)/(b-a)
            for q in (-1/np.sqrt(3), 1/np.sqrt(3)):
                r = (a+b)/2+q*(b-a)/2
                t = ta+slope*(r-a)
                total += (b-a)/2*np.sin(t)*(r*r*slope*sample(f[term+'_r'], r, t, rc, tc)
                                                        -r*sample(f[term+'_theta'], r, t, rc, tc))
        row[term+'_jet'] = factor*total
        row[term+'_acc_net'] = row[term+'_acc_outer']-row[term+'_acc_inner']
    # Cut-cell volumes use theta(r_center); exact for fixed-angle control volumes.
    rv = np.maximum(np.minimum(rf[1:], rhi)**3-np.maximum(rf[:-1], rlo)**3, 0)/3
    limit = np.maximum(tf[:-1, None], surf[None, :])
    angular = np.where(tf[1:, None] > limit, np.cos(limit)-np.cos(tf[1:, None]), 0)
    volume = factor*angular*rv[None, :]
    for name, key in (('disk_mass', 'density'), ('disk_angular_momentum', 'angular_density'),
                      ('disk_energy_with_gravity', 'energy_density'), ('ohmic_cooling_power', 'cooling')):
        row[name] = float(np.sum(np.where(volume > 0, volume*f[key], 0)))
    # Sec.4 inner disk interval 1.1..1.5; radial-length average explicitly declared.
    ml, mh = args.mu_inner, args.mu_outer
    if not rf[0] <= ml < mh <= rf[-1]:
        raise ValueError('Magnetization averaging interval must lie inside the mesh')
    weight = np.maximum(np.minimum(rf[1:], mh)-np.maximum(rf[:-1], ml), 0)
    row['mu_mid_radial_mean'] = float(np.sum(mu*weight)/sum(weight))
    row['eta_p_mid_radial_mean'] = float(np.sum(eta*weight)/sum(weight))
    row['mass_ejection_fraction'] = row['mass_jet']/row['mass_acc_outer'] if row['mass_acc_outer'] else np.nan
    return row


def budgets(rows, moving):
    time = np.array([r['time'] for r in rows])
    for name, stock, terms in (('mass', 'disk_mass', ('mass',)),
                               ('angular', 'disk_angular_momentum', ('j_adv', 'j_mag')),
                               ('energy', 'disk_energy_with_gravity', ('e_kin', 'e_grav', 'e_enthalpy', 'e_poynting'))):
        rate = np.array([sum(r[t+'_acc_net']-r[t+'_jet'] for t in terms)
                         -(r['ohmic_cooling_power'] if name == 'energy' else 0) for r in rows])
        integral = np.r_[0., np.cumsum(.5*(rate[1:]+rate[:-1])*np.diff(time))]
        initial = rows[0][stock]
        for k, row in enumerate(rows):
            row[name+'_snapshot_integrated_supply'] = integral[k]
            row[name+'_snapshot_budget_residual'] = np.nan if moving else row[stock]-initial-integral[k]


def plot(rows, path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    t = [r['time'] for r in rows]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for key in ('mass_acc_outer', 'mass_acc_inner', 'mass_jet'):
        axes[0, 0].plot(t, [r[key] for r in rows], label=key)
    for key in ('j_adv_jet', 'j_mag_jet'):
        axes[0, 1].plot(t, [r[key] for r in rows], label=key)
    for key in ('e_kin_jet', 'e_grav_jet', 'e_enthalpy_jet', 'e_poynting_jet'):
        axes[1, 0].plot(t, [r[key] for r in rows], label=key)
    for key in ('mu_mid_radial_mean', 'eta_p_mid_radial_mean'):
        axes[1, 1].plot(t, [r[key] for r in rows], label=key)
    for ax in axes.flat:
        ax.set_xlabel('time [r0/vK0]'); ax.grid(alpha=.25); ax.legend(fontsize=8)
    fig.suptitle('SF2014 control-volume snapshot estimates (not discrete solver fluxes)')
    fig.savefig(path, dpi=160); plt.close(fig)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run_dir', type=Path)
    p.add_argument('--output', type=Path, help='Output directory (default RUN/paper_diagnostics)')
    p.add_argument('--inner', type=float, default=1.)
    p.add_argument('--outer', type=float, default=10.)
    p.add_argument('--start-time', type=float, default=0., help='Ignore earlier snapshots (useful for vr-zero surfaces)')
    p.add_argument('--surface', choices=('fixed', 'vr-zero'), default='fixed')
    p.add_argument('--theta', type=float, help='Fixed angle in radians (default pi/2-2 epsilon)')
    p.add_argument('--mu-inner', type=float, default=1.1)
    p.add_argument('--mu-outer', type=float, default=1.5)
    p.add_argument('--bipolar', action='store_true', help='Multiply upper-hemisphere 2pi integrals by two, as in paper')
    args = p.parse_args()
    if args.surface != 'fixed' and args.theta is not None:
        p.error('--theta applies only to --surface fixed')
    snapshots = (load_snapshot(args.run_dir, str(path.resolve())) for path in snapshot_files(args.run_dir))
    rows = [evaluate(s, args) for s in snapshots if s.time >= args.start_time]
    if not rows:
        raise ValueError('No snapshots remain after --start-time selection')
    if any(b['time'] <= a['time'] for a, b in zip(rows, rows[1:])):
        raise ValueError('Snapshot times must be strictly increasing')
    budgets(rows, args.surface != 'fixed')
    out = args.output or args.run_dir/'paper_diagnostics'
    out.mkdir(parents=True, exist_ok=True)
    with (out/'fluxes.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    metadata = {'source': 'https://arxiv.org/html/1407.3136', 'sections': '3.1, 3.3, 4.1, Appendix B',
                'surface': args.surface, 'fixed_theta': rows[0]['surface_theta_min'] if args.surface == 'fixed' else None, 'normalization': '4pi bipolar' if args.bipolar else '2pi upper hemisphere',
                'mu_average': 'radial length weighted, last-row pressure and Btheta, Br=0 by symmetry',
                'mu_interval': [args.mu_inner, args.mu_outer], 'control_radii': [args.inner, args.outer],
                'limitations': ['Cell snapshot interpolation/extrapolation; no authoritative face fields or HLL fluxes',
                                'Resistive current from centered cell-field finite differences',
                                'Fixed angle approximates paper vr-zero surface',
                                'Time integrals use snapshot trapezoids; floors and CT boundary energy adjustments absent',
                                'Moving-surface budget residual disabled: no Reynolds transport correction',
                                'Equatorial flux set to zero by symmetry; no zero-area pole contribution']}
    (out/'method.json').write_text(json.dumps(metadata, indent=2)+'\n')
    plotted = plot(rows, out/'fluxes.png')
    print(f'Wrote {out / "fluxes.csv"} and method.json'+(' and fluxes.png' if plotted else ' (matplotlib unavailable)'))


if __name__ == '__main__':
    main()
