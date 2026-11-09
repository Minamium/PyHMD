"""Immutable configuration for the native capabilities actually implemented."""
from __future__ import annotations
from dataclasses import asdict, dataclass, fields, replace
from importlib.resources import files
import math
from pathlib import Path
from typing import Mapping

PROBLEMS = ('thin_disk', 'uniform', 'alfven', 'diffusion', 'brio_wu',
            'oblique_alfven', 'linear_wave', 'field_loop', 'rj2a', 'orszag_tang', 'blast_low_beta')

def _problem(value):
    if not isinstance(value,str): raise ValueError('problem must be a string')
    value=value.lower()
    return {'disk':'thin_disk','brio-wu':'brio_wu'}.get(value,value)

def _choice(value, choices, name):
    if not isinstance(value, str) or value not in choices:
        raise ValueError(f'{name} must be one of {choices}; got {value!r}')

def _int(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')

def _number(value, name, minimum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'{name} must be a finite number')
    if minimum is not None and value < minimum:
        raise ValueError(f'{name} must be >= {minimum}')

def _typed(cls, value, defaults=None):
    if isinstance(value, cls):
        return value
    if value is None:
        value = {}
    if not isinstance(value, Mapping):
        raise TypeError(f'{cls.__name__} must be an object')
    allowed = {f.name for f in fields(cls)}
    unknown = set(value) - allowed
    if unknown:
        raise ValueError(f'unknown {cls.__name__} fields: {sorted(unknown)}')
    return cls(**{**(defaults or {}), **dict(value)})

@dataclass(frozen=True)
class GridConfig:
    geometry: str = 'spherical'
    nx1: int = 600
    nx2: int = 128
    radial_grid: str = 'log'
    x1_min: float = 1.
    x1_max: float = 1500.
    x2_min: float = 0.
    x2_max: float = math.pi/2

    def __post_init__(self):
        _choice(self.geometry, ('cartesian', 'spherical'), 'geometry')
        _choice(self.radial_grid, ('uniform', 'log'), 'radial_grid')
        for name in ('nx1', 'nx2'): _int(getattr(self, name), name, 3)
        for name in ('x1_min', 'x1_max', 'x2_min', 'x2_max'): _number(getattr(self, name), name)
        if self.x1_max <= self.x1_min or self.x2_max <= self.x2_min:
            raise ValueError('grid upper bounds must exceed lower bounds')
        if (self.geometry == 'spherical' or self.radial_grid == 'log') and self.x1_min <= 0:
            raise ValueError('spherical/log grid requires x1_min > 0')

@dataclass(frozen=True)
class PhysicsConfig:
    gamma: float = 5/3
    gm: float = 1.
    eta_model: str = 'sf2014'
    eta: float = .001
    h: float = .1
    corona: float = .001
    mu0: float = .01
    mu_ref: float = .01
    b0: float | None = None
    bend: float = .5
    alpha_m: float = 1.65
    chi: float = .5
    cooling: bool = True
    initial_equilibrium: str = 'magnetohydrostatic'
    wave_amplitude: float = .1
    wave_mode: str = 'fast'
    noise_seed: int = 0
    noise_amplitude: float = 0.

    def __post_init__(self):
        for f in fields(self):
            if f.name not in ('eta_model', 'cooling', 'initial_equilibrium', 'b0', 'wave_mode'):
                _number(getattr(self, f.name), f.name, 0)
        if self.b0 is not None: _number(self.b0, 'b0')
        if self.gamma <= 1 or not 0 < self.h < 1 or min(self.corona, self.mu_ref, self.bend) <= 0:
            raise ValueError('invalid gamma, h, corona, mu_ref or bend')
        _choice(self.eta_model, ('constant', 'sf2014', 'strong'), 'eta_model')
        _choice(self.initial_equilibrium, ('hydrostatic', 'magnetohydrostatic'), 'initial_equilibrium')
        if not isinstance(self.cooling, bool): raise ValueError('cooling must be bool')
        _choice(self.wave_mode, ('fast','alfven','slow','entropy'), 'wave_mode')
        _int(self.noise_seed, 'noise_seed')
        if self.noise_amplitude>1e-6: raise ValueError('noise_amplitude must be <= 1e-6')

@dataclass(frozen=True)
class NumericsConfig:
    rho_floor: float = 1e-10
    p_floor: float = 1e-12
    thermal_recovery: str = 'internal_energy'
    internal_energy_sync: str = 'modern'
    dual_energy_eta: float = 1e-3
    entropy_ratio: float = 1e-3
    entropy_shock_threshold: float = .05
    entropy_contact_threshold: float = .1
    ct_energy_correction: bool = True
    max_step_retries: int = 8
    retry_flat_same_dt: bool = True
    negative_pressure_policy: str = 'retry'
    reconstruction_cache: bool = True
    metric_ppm: bool = True
    ct_emf_scheme: str = 'face_hll'
    cfl_staggered: bool = True
    retry_streak_limit: int = 64
    retry_streak_threshold: int = 4

    def __post_init__(self):
        for name in ('rho_floor', 'p_floor', 'dual_energy_eta', 'entropy_ratio', 'entropy_shock_threshold', 'entropy_contact_threshold'):
            _number(getattr(self, name), name)
            if getattr(self, name) <= 0: raise ValueError(f'{name} must be positive')
        for name in ('ct_energy_correction', 'retry_flat_same_dt', 'reconstruction_cache', 'metric_ppm', 'cfl_staggered'):
            if not isinstance(getattr(self, name), bool): raise ValueError(f'{name} must be bool')
        for name in ('max_step_retries', 'retry_streak_limit', 'retry_streak_threshold'): _int(getattr(self, name), name)
        _choice(self.thermal_recovery, ('off', 'selective', 'dual_energy', 'internal_energy'), 'thermal_recovery')
        _choice(self.internal_energy_sync, ('legacy', 'modern'), 'internal_energy_sync')
        _choice(self.negative_pressure_policy, ('retry', 'floor'), 'negative_pressure_policy')
        _choice(self.ct_emf_scheme, ('face_hll', 'corner_hll'), 'ct_emf_scheme')

@dataclass(frozen=True)
class BoundaryConfig:
    lower1: str = 'paper_inner'
    upper1: str = 'paper_outer'
    lower2: str = 'axis'
    upper2: str = 'equator'
    inner_normal_field: str = 'evolve'

    def __post_init__(self):
        choices = ('outflow', 'periodic', 'reflect', 'axis', 'equator', 'paper_inner', 'paper_outer')
        for name in ('lower1', 'upper1', 'lower2', 'upper2'): _choice(getattr(self, name), choices, name)
        _choice(self.inner_normal_field, ('evolve', 'project'), 'inner_normal_field')

@dataclass(frozen=True)
class TimeConfig:
    cfl: float = .4
    t_end: float = 100.
    max_steps: int = 100000000

    def __post_init__(self):
        _number(self.cfl, 'cfl'); _number(self.t_end, 't_end'); _int(self.max_steps, 'max_steps', 1)
        if not 0 < self.cfl <= .5 or self.t_end <= 0: raise ValueError('require 0 < cfl <= .5 and t_end > 0')

GROUPS = {'grid': GridConfig, 'physics': PhysicsConfig, 'numerics': NumericsConfig, 'boundary': BoundaryConfig, 'time': TimeConfig}
FIELD_GROUP = {f.name: name for name, cls in GROUPS.items() for f in fields(cls)}

def read_ini(path):
    """Match the native flat-key INI contract, including duplicate rejection."""
    values = {}
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        line = line.split('#', 1)[0].split(';', 1)[0].strip()
        if not line: continue
        if line.startswith('['):
            if not line.endswith(']'): raise ValueError(f'{path}:{number}: malformed section')
            continue
        if '=' not in line: raise ValueError(f'{path}:{number}: expected key=value')
        key, value = (q.strip() for q in line.split('=', 1)); key = key.lower()
        if not key or not value or key in values: raise ValueError(f'{path}:{number}: empty/duplicate key {key!r}')
        values[key] = value[1:-1] if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'" else value
    return values

def _base(problem):
    data = {name: asdict(cls()) for name, cls in GROUPS.items()}
    if problem != 'thin_disk':
        data['grid'].update(geometry='cartesian', nx1=128, nx2=4, radial_grid='uniform', x1_min=0., x1_max=1., x2_max=.125)
        data['physics'].update(gm=0., eta=0., eta_model='constant', cooling=False)
        data['numerics'].update(thermal_recovery='off', internal_energy_sync='legacy', ct_energy_correction=False, ct_emf_scheme='corner_hll')
        data['boundary'].update(lower1='outflow', upper1='outflow', lower2='periodic', upper2='periodic')
    return data

def _convert(key, value, default):
    if isinstance(default, bool):
        choices = {'yes': True, 'true': True, '1': True, 'no': False, 'false': False, '0': False, 'none': False, 'on': True, 'off': False}
        if str(value).lower() not in choices: raise ValueError(f'{key} must be a boolean')
        return choices[str(value).lower()]
    if isinstance(default, int): return int(value)
    if isinstance(default, float) or key == 'b0': return float(value)
    return str(value).lower()

def _canonical(flat):
    flat = dict(flat)
    for alias, key in [('entropy_switch', 'thermal_recovery')]:
        if alias in flat:
            if key in flat: raise ValueError(f'{alias} and {key} duplicate the same setting')
            flat[key] = flat.pop(alias)
    if 'log_r' in flat:
        if 'radial_grid' in flat: raise ValueError('log_r and radial_grid duplicate the same setting')
        flat['radial_grid'] = 'log' if _convert('log_r', flat.pop('log_r'), True) else 'uniform'
    aliases = {'thermal_recovery': {'0':'off','no':'off','1':'selective','yes':'selective','2':'dual_energy','3':'internal_energy'},
               'eta_model': {'0':'constant','1':'sf2014','paper':'sf2014','2':'strong'},
               'initial_equilibrium': {'0':'hydrostatic','1':'magnetohydrostatic','lorentz_balanced':'magnetohydrostatic'},
               'negative_pressure_policy': {'0':'retry','strict':'retry','1':'floor'},
               'radial_grid': {'logarithmic':'log'}}
    for key, mapping in aliases.items():
        if key in flat: flat[key] = mapping.get(str(flat[key]).lower(), str(flat[key]).lower())
    return flat

def _apply_flat(data, flat, problem):
    for key, value in _canonical(flat).items():
        if key == 'case' or key in RESULT_KEYS or key == 'restart': continue
        if key not in FIELD_GROUP: raise ValueError(f'unknown native configuration key: {key}')
        group = FIELD_GROUP[key]; data[group][key] = _convert(key, value, data[group][key])
    if problem=='thin_disk' and 'b0' in flat and 'mu0' not in flat:
        p = data['physics']; p['mu0'] = p['b0']**2/(2*p['h']**2*p['gm'])

def preset(problem):
    problem=_problem(problem)
    _choice(problem, PROBLEMS, 'problem')
    flat = read_ini(files('pyhmd.presets').joinpath(problem+'.ini'))
    data = _base(problem); _apply_flat(data, flat, problem)
    return data, flat

@dataclass(frozen=True)
class Config:
    problem: str = 'thin_disk'
    model: str = 'resistive_mhd'
    scheme: str = 'ppm_hll_ssprk3'
    dtype: str = 'float64'
    grid: GridConfig | Mapping | None = None
    physics: PhysicsConfig | Mapping | None = None
    numerics: NumericsConfig | Mapping | None = None
    boundary: BoundaryConfig | Mapping | None = None
    time: TimeConfig | Mapping | None = None

    def __post_init__(self):
        object.__setattr__(self,'problem',_problem(self.problem))
        _choice(self.problem, PROBLEMS, 'problem')
        _choice(self.model, ('resistive_mhd',), 'model')
        _choice(self.scheme, ('ppm_hll_ssprk3',), 'scheme')
        _choice(self.dtype, ('float64',), 'dtype')
        defaults, _ = preset(self.problem)
        for name, cls in GROUPS.items(): object.__setattr__(self, name, _typed(cls, getattr(self, name), defaults[name]))
        if self.problem == 'thin_disk' and self.grid.geometry != 'spherical': raise ValueError('thin_disk requires spherical geometry')
        if self.physics.noise_amplitude and self.problem != 'thin_disk': raise ValueError('noise is supported only for thin_disk')
        if self.physics.wave_mode != 'fast' and self.problem != 'linear_wave': raise ValueError('wave_mode applies only to linear_wave')
        if self.problem=='linear_wave' and not 0 < self.physics.wave_amplitude <= 1e-3: raise ValueError('linear_wave requires 0 < wave_amplitude <= 1e-3')

    @classmethod
    def from_dict(cls, value): return _typed(cls, value)
    def as_dict(self): return asdict(self)
    def native_settings(self):
        settings = {'case': self.problem}
        for name in GROUPS: settings.update(asdict(getattr(self, name)))
        return {k:v for k,v in settings.items() if v is not None}

@dataclass(frozen=True)
class ResultOptions:
    out_dir: str | Path | None = None
    output_interval: float | None = None
    log_interval: int | None = None
    checkpoint_interval: float = 0.
    checkpoint_on_retry: bool = True
    transport_diagnostics: bool = False
    midplane_diagnostics: bool = False

    def __post_init__(self):
        if self.out_dir is not None and (not isinstance(self.out_dir, (str, Path)) or not str(self.out_dir)):
            raise ValueError('out_dir must be a nonempty path')
        if self.output_interval is not None:
            _number(self.output_interval, 'output_interval')
            if self.output_interval <= 0: raise ValueError('output_interval must be positive')
        if self.log_interval is not None: _int(self.log_interval, 'log_interval', 1)
        _number(self.checkpoint_interval, 'checkpoint_interval', 0)
        for key in ('checkpoint_on_retry', 'transport_diagnostics', 'midplane_diagnostics'):
            if not isinstance(getattr(self, key), bool): raise ValueError(f'{key} must be bool')

    def resolved(self, problem):
        _, flat = preset(problem)
        return replace(self, output_interval=self.output_interval if self.output_interval is not None else float(flat.get('output_interval',10)),
                       log_interval=self.log_interval if self.log_interval is not None else int(flat.get('log_interval',100)))

RESULT_KEYS = {'output', 'output_interval', 'log_interval', 'checkpoint_interval', 'checkpoint_on_retry', 'transport_diagnostics', 'midplane_diagnostics'}

def config_from_ini(path):
    flat = read_ini(path); problem = _problem(flat.get('case', 'thin_disk'))
    data, _ = preset(problem); _apply_flat(data, flat, problem)
    cfg = Config(problem=problem, **data)
    options = {}
    for key in RESULT_KEYS & flat.keys():
        target = 'out_dir' if key == 'output' else key
        defaults = {'output': '', 'output_interval': 10., 'log_interval': 100, 'checkpoint_interval': 0., 'checkpoint_on_retry': True, 'transport_diagnostics': False, 'midplane_diagnostics': False}
        options[target] = flat[key] if key == 'output' else _convert(key, flat[key], defaults[key])
    return cfg, ResultOptions(**options), flat.get('restart')

def native_ini(settings):
    def value(x):
        if isinstance(x, bool): return 'yes' if x else 'no'
        return '"'+x+'"' if isinstance(x,str) else str(x)
    for key, item in settings.items():
        if any(c in str(item) for c in '\n\r#;'):
            raise ValueError(f'{key} contains characters unsupported by the native INI format')
    return '# Generated by PyHMD; fully resolved native settings.\n[run]\n'+''.join(f'{k} = {value(v)}\n' for k,v in settings.items())
