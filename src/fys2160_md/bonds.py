"""Explicit intramolecular bonds; lengths and energies are in reduced units."""
from dataclasses import dataclass, asdict
import math


def _positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive.")
    return value


@dataclass(frozen=True)
class HarmonicBond:
    """U = stiffness / 2 * (r - length)**2."""
    length: float = 0.7
    stiffness: float = 100.0

    def __post_init__(self):
        object.__setattr__(self, 'length', _positive(self.length, 'length'))
        object.__setattr__(self, 'stiffness', _positive(self.stiffness, 'stiffness'))


@dataclass(frozen=True)
class RigidBond:
    """Fixed separation, enforced with position and velocity constraints."""
    length: float = 0.7

    def __post_init__(self):
        object.__setattr__(self, 'length', _positive(self.length, 'length'))


@dataclass(frozen=True)
class Class2Bond:
    """Legacy anharmonic bond: U = k2*q**2 + k3*q**3 + k4*q**4."""
    length: float = 0.7
    k2: float = 35.0
    k3: float = -40.0
    k4: float = 59.0

    def __post_init__(self):
        for name in ('length', 'k2', 'k4'):
            object.__setattr__(self, name, _positive(getattr(self, name), name))
        if not math.isfinite(float(self.k3)):
            raise ValueError('k3 must be finite.')
        object.__setattr__(self, 'k3', float(self.k3))


def resolve_model(model=None, molecule=None, bond=None):
    if model is not None and molecule is not None:
        raise ValueError('Specify molecule and bond, or the legacy model, not both.')
    if bond is not None and not isinstance(bond, (HarmonicBond, RigidBond, Class2Bond)):
        raise TypeError('bond must be HarmonicBond, RigidBond or Class2Bond.')
    if model is not None:
        if model not in ('atomic', 'diatomic-flexible', 'diatomic-rigid'):
            raise ValueError('Unknown model.')
        if bond is None:
            bond = {'atomic': None, 'diatomic-flexible': Class2Bond(), 'diatomic-rigid': RigidBond()}[model]
        expected = 'atomic' if bond is None else ('diatomic-rigid' if isinstance(bond, RigidBond) else 'diatomic-flexible')
        if model != expected:
            raise ValueError('model and bond disagree.')
        return model, bond
    if molecule not in (None, 'diatomic'):
        raise ValueError('molecule must be None (atoms) or "diatomic".')
    if molecule is None:
        if bond is not None:
            raise ValueError('Set molecule="diatomic" when specifying a bond.')
        return 'atomic', None
    bond = HarmonicBond() if bond is None else bond
    return ('diatomic-rigid' if isinstance(bond, RigidBond) else 'diatomic-flexible'), bond


def bond_config(bond):
    return None if bond is None else {'kind': type(bond).__name__, **asdict(bond)}


def read_bond(config):
    if config is None:
        return None
    args = config.copy()
    kind = args.pop('kind')
    return {'HarmonicBond': HarmonicBond, 'RigidBond': RigidBond, 'Class2Bond': Class2Bond}[kind](**args)
