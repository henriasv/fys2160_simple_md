"""Compare BKS alpha-quartz in a tilted box with LAMMPS: energy, forces and pressure tensor.

Usage: python scripts/validate_quartz.py [path to lmp_serial with KSPACE]

The crystal is sheared (all three tilts), its atoms are moved off their sites and
left at rest, so the pressure tensor is the virial alone. Both codes use the same
Buckingham terms with the dispersion beyond the cutoff added as a tail correction
(no shift), real-space Ewald pairs without tables, and PPPM of order 5 with our mesh
and splitting parameter.
LAMMPS is a reference for behaviour only; nothing is copied from it.

Expect the energy to agree to about 1e-14, forces to about 1e-5 eV/Å, the off-diagonal
pressures to well under a bar and the diagonal ones to differ by about 3 bar (of tens
of kbar). That difference is in the real-space Coulomb force: both codes use the same
approximation of erfc (Abramowitz and Stegun 7.1.26, LAMMPS without tables), but
LAMMPS pairs it with the derivative of the exact erfc, while our force is the exact
derivative of the energy we report. Summing that difference over the pairs predicts
the 3 bar. The PPPM virial itself differs from the strain derivative of the energy
by a few bar in both codes, as the influence function is held fixed.
"""
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np
from fys2160_md import AlphaQuartz, Buckingham, MDSimulation

lammps = sys.argv[1] if len(sys.argv) > 1 else 'lmp_serial'
bks = {('O', 'O'): Buckingham(A=1388.7730, rho=0.362319, C=175.0),
       ('Si', 'O'): Buckingham(A=18003.7572, rho=0.205205, C=133.5381)}
charge = {'Si': 2.4, 'O': -1.2}
crystal = AlphaQuartz(cells=(5, 3, 4))
crystal.strain(xy=0.07, xz=-0.05, yz=0.06, xx=0.01)
rng = np.random.default_rng(11)
crystal.positions += rng.normal(scale=0.03, size=crystal.positions.shape)
crystal.velocities = np.zeros_like(crystal.positions)

sim = MDSimulation(units='metal', pair_potential=bks, charge=charge, coulomb='pppm', ensemble='nve',
                   temperature=300, timestep=0.001, cutoff=8.0, kspace_accuracy=1e-6, lj_shift=False, tail=True)
engine = sim._engine(crystal)
ours = engine.observables
forces = engine._f.copy()
mesh, g = engine._pppm_state[:3].astype(int), engine._pppm_state[3]

(lx, ly, lz), (xy, xz, yz) = crystal.box.lengths, crystal.tilt
x, species = engine._x, crystal.atoms.species
types = np.where(species == 'Si', 1, 2)
with tempfile.TemporaryDirectory() as folder:
    folder = Path(folder)
    lines = ['tilted BKS quartz', '', f'{len(x)} atoms', '2 atom types', '',
             f'0 {lx:.17g} xlo xhi', f'0 {ly:.17g} ylo yhi', f'0 {lz:.17g} zlo zhi', f'{xy:.17g} {xz:.17g} {yz:.17g} xy xz yz', '',
             'Masses', '', '1 28.0855', '2 15.9994', '', 'Atoms # charge', '']
    lines += [f'{i+1} {t} {charge["Si" if t == 1 else "O"]} {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}' for i, (t, p) in enumerate(zip(types, x))]
    (folder/'data').write_text('\n'.join(lines)+'\n')
    (folder/'in').write_text('\n'.join([
        'units metal', 'atom_style charge', 'boundary p p p', 'read_data data',
        'pair_style buck/coul/long 8.0', 'pair_coeff 1 1 0.0 1.0 0.0',
        'pair_coeff 1 2 18003.7572 0.205205 133.5381', 'pair_coeff 2 2 1388.7730 0.362319 175.0',
        'pair_modify shift no tail yes table 0',
        f'kspace_style pppm 1e-6', f'kspace_modify mesh {mesh[0]} {mesh[1]} {mesh[2]} order 5 gewald {g:.17g}',
        'neighbor 2.0 bin', 'velocity all set 0 0 0',
        'thermo_style custom pe pxx pyy pzz pxy pxz pyz', 'thermo_modify format float %.15g norm no',
        'dump forces all custom 1 forces id fx fy fz', 'dump_modify forces sort id format float %.15g',
        'run 0']))
    subprocess.run([lammps, '-in', 'in', '-log', 'log', '-screen', 'none'], cwd=folder,
                   capture_output=True, text=True, check=True)
    log = (folder/'log').read_text().splitlines()
    header = next(i for i, line in enumerate(log) if line.split()[:1] == ['PotEng'])
    pe, pxx, pyy, pzz, pxy, pxz, pyz = (float(v) for v in log[header+1].split())
    dump = (folder/'forces').read_text().splitlines()
    theirs = np.array([[float(v) for v in line.split()[1:]] for line in dump[9:]])

print(f'{len(x)} atoms, tilt {np.round(crystal.tilt, 3)}, PPPM mesh {mesh.tolist()}, g {g:.6f}')
print(f'energy: ours {engine._potential:.10f} eV, LAMMPS {pe:.10f} eV, relative difference '
      f'{abs(engine._potential-pe)/abs(pe):.1e}')
print(f'forces: largest difference {np.abs(forces-theirs).max():.1e} eV/Å (largest force {np.abs(theirs).max():.1f})')
for name, value in zip(('pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz'), (pxx, pyy, pzz, pxy, pxz, pyz)):
    print(f'{name}: ours {ours[name]: .6f} bar, LAMMPS {value: .6f} bar, difference {abs(ours[name]-value):.1e}')
