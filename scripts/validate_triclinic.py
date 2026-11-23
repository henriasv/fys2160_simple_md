"""Compare a tilted Lennard-Jones box with LAMMPS: energy, forces and pressure tensor, at rest.

Usage: python scripts/validate_triclinic.py [path to lmp_serial]

LAMMPS is a reference for behaviour only; nothing is copied from it. The atoms are
at rest, so the pressure tensor is the pair virial alone, which both codes report
for the same shifted 12-6 potential with cutoff 2.5.
"""
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np
from fys2160_md import FCC, MDSimulation

lammps = sys.argv[1] if len(sys.argv) > 1 else 'lmp_serial'
crystal = FCC(N=500, rho=0.95)
crystal.strain(xy=0.21, xz=-0.13, yz=0.17)
rng = np.random.default_rng(7)
crystal.positions += rng.normal(scale=0.05, size=crystal.positions.shape)  # leave the lattice a little
crystal.velocities = np.zeros_like(crystal.positions)

engine = MDSimulation(pair_potential='lj', ensemble='nve', temperature=1.0, timestep=0.005, cutoff=2.5)._engine(crystal)
ours = engine.observables
engine._advance(0)
forces = engine._f.copy()

(lx, ly, lz), (xy, xz, yz) = crystal.box.lengths, crystal.tilt
x = engine._x
with tempfile.TemporaryDirectory() as folder:
    folder = Path(folder)
    lines = ['tilted LJ', '', f'{len(x)} atoms', '1 atom types', '',
             f'0 {lx:.17g} xlo xhi', f'0 {ly:.17g} ylo yhi', f'0 {lz:.17g} zlo zhi', f'{xy:.17g} {xz:.17g} {yz:.17g} xy xz yz', '',
             'Masses', '', '1 1.0', '', 'Atoms # atomic', '']
    lines += [f'{i+1} 1 {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}' for i, p in enumerate(x)]
    (folder/'data').write_text('\n'.join(lines)+'\n')
    (folder/'in').write_text('\n'.join([
        'units lj', 'atom_style atomic', 'boundary p p p', 'read_data data',
        'pair_style lj/cut 2.5', 'pair_coeff 1 1 1.0 1.0', 'pair_modify shift yes',
        'neighbor 0.4 bin', 'velocity all set 0 0 0',
        'thermo_style custom pe pxx pyy pzz pxy pxz pyz', 'thermo_modify format float %.15g norm no',
        'dump forces all custom 1 forces id fx fy fz', 'dump_modify forces sort id format float %.15g',
        'run 0']))
    output = subprocess.run([lammps, '-in', 'in', '-log', 'log', '-screen', 'none'], cwd=folder,
                            capture_output=True, text=True, check=True)
    log = (folder/'log').read_text().splitlines()
    header = next(i for i, line in enumerate(log) if line.split()[:1] == ['PotEng'])
    pe, pxx, pyy, pzz, pxy, pxz, pyz = (float(v) for v in log[header+1].split())
    dump = (folder/'forces').read_text().splitlines()
    theirs = np.array([[float(v) for v in line.split()[1:]] for line in dump[9:]])

print(f'energy: ours {engine._potential:.12f}, LAMMPS {pe:.12f}, relative difference '
      f'{abs(engine._potential-pe)/abs(pe):.1e}')
print(f'forces: largest difference {np.abs(forces-theirs).max():.1e} (largest force {np.abs(theirs).max():.1f})')
for name, value in zip(('pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz'), (pxx, pyy, pzz, pxy, pxz, pyz)):
    print(f'{name}: ours {ours[name]: .12f}, LAMMPS {value: .12f}, difference {abs(ours[name]-value):.1e}')
