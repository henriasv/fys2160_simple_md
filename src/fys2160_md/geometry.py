"""Periodic geometry constructors, independent of the simulation runner."""
import numpy as np
from .bonds import resolve_model, bonds_for_species, RigidWater
from .species import species_counts, parameter, per_atom
from .system import System, MODELS, _integer, _positive


def _box(N, rho, box, model):
    N=_integer(N,'N',2)
    if model not in MODELS:
        raise ValueError(f'model must be one of {tuple(MODELS)}')
    if N*(1 if model=='atomic' else 2)>20000:
        raise ValueError('This course solver supports at most 20000 atoms.')
    if rho is not None and box is not None:
        raise ValueError('Specify rho or box, not both; N fixes the remaining quantity.')
    if box is None:
        rho=_positive(.001 if rho is None else rho,'rho')
        lengths=np.full(3,(N/rho)**(1/3))
    else:
        lengths=np.broadcast_to(np.asarray(box,dtype=float),(3,)).copy()
        if not np.isfinite(lengths).all() or np.any(lengths<=0):
            raise ValueError('box must be a positive side length or three positive lengths.')
    return N,lengths


def _rng(seed):
    seed=_integer(seed,'seed')
    if seed>=2**64:raise ValueError('seed must be below 2**64.')
    return np.random.default_rng(seed)


def _atoms(centers, lengths, model, rng, bond_lengths):
    if model=='atomic':return centers
    directions=rng.normal(size=(len(centers),3))
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    half_length=np.asarray(bond_lengths).reshape(-1,1)/2
    return np.stack((centers-half_length*directions,centers+half_length*directions),axis=1).reshape(-1,3)%lengths


def FCC(*, N=None, rho=None, box=None, molecule=None, bond=None, model=None, species=None, masses=1., seed=87287):
    """Perfect periodic FCC sites at rho+N or box+N; the default box is cubic.

    A scalar box is a cube side length; a triple gives rectangular box lengths.
    N counts atoms, or molecular centres for a diatomic model. Cubes require
    N=4*n**3. Rectangular boxes must fit integer conventional cells of equal size.
    """
    N, labels = species_counts(N, species)
    model, bond = resolve_model(model, molecule, bond)
    N,lengths=_box(N,rho,box,model)
    specs = () if model=='atomic' else bonds_for_species(bond, labels)
    if any(b.length >= min(lengths)/2 for b in specs):
        raise ValueError('Bond length must be less than half the shortest box side.')
    lattice_constant=(4*np.prod(lengths)/N)**(1/3)
    ratios=lengths/lattice_constant
    cells=np.rint(ratios).astype(int)
    if np.any(cells<1) or 4*np.prod(cells)!=N or not np.allclose(ratios,cells,rtol=1e-10,atol=1e-10):
        if np.allclose(lengths,lengths[0],rtol=1e-10,atol=1e-10):
            lower=max(1,int(np.floor(np.cbrt(N/4))))
            raise ValueError(f'A cubic FCC crystal needs N=4*n^3; N={N} is incompatible. Nearby counts: {4*lower**3} or {4*(lower+1)**3}.')
        raise ValueError('A rectangular FCC box must contain integer numbers of equal-sized conventional cells, with N=4*nx*ny*nz. Adjust N or the box aspect ratio.')
    rng=_rng(seed)
    if len(set(labels)) > 1: rng.shuffle(labels)
    bond_lengths = [] if model=='atomic' else [b.length for b in bonds_for_species(bond, labels)]
    labels = np.repeat(labels, 1 if model == 'atomic' else 2)
    mass_values = per_atom(parameter(masses, tuple(dict.fromkeys(labels.tolist())), 'masses'), labels)
    bases=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
    sites=(np.indices(tuple(cells)).reshape(3,-1).T[:,None,:]+bases).reshape(-1,3)
    centers=(sites+.25)*lattice_constant
    return System(_atoms(centers,lengths,model,rng,bond_lengths),lengths,model=model,bond=bond,species=labels,masses=mass_values)


def _water_molecules(centers, water, rng, masses):
    """Randomly oriented rigid waters with centres of mass at the given points."""
    half = np.radians(water.hoh_angle/2)
    template = np.array([[0, 0, 0], [np.sin(half), np.cos(half), 0], [-np.sin(half), np.cos(half), 0]])*water.oh_length
    weights = np.array([masses['O'], masses['H'], masses['H']])
    template -= weights @ template/weights.sum()
    q = rng.normal(size=(len(centers), 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    w, x, y, z = q.T
    rotation = np.stack([np.stack([1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)], -1),
                         np.stack([2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)], -1),
                         np.stack([2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)], -1)], 1)
    return (centers[:, None, :] + np.einsum('nij,aj->nai', rotation, template)).reshape(-1, 3)


def Random(*, N=None, rho=None, box=None, min_distance=None, max_attempts=1000,
           molecule=None, bond=None, model=None, species=None, masses=None, seed=87287):
    """Self-avoiding random placement with periodic minimum-image distances.

    min_distance is the minimum separation between atoms of different molecules,
    in the system's length unit (default 0.9 reduced units, or 1.8 Å for water).
    Failed packing raises an error rather than reducing this distance.
    max_attempts bounds trials for each particle.

    molecule="water" places randomly oriented rigid waters (species "O" and "H")
    in LAMMPS real units: N molecules, box in Å or rho in molecules per Å^3
    (0.0334 is about 1 g/cm^3), masses in g/mol.
    """
    if molecule == 'water' or model == 'water-rigid':
        return _random_water(N=N, rho=rho, box=box, min_distance=1.8 if min_distance is None else min_distance,
                             max_attempts=max_attempts, model=model, molecule=molecule, bond=bond, species=species,
                             masses=masses, seed=seed)
    min_distance = .9 if min_distance is None else min_distance
    masses = 1. if masses is None else masses
    N, labels = species_counts(N, species)
    model, bond = resolve_model(model, molecule, bond)
    N,lengths=_box(N,rho,box,model)
    specs = () if model=='atomic' else bonds_for_species(bond, labels)
    if any(b.length >= min(lengths)/2 for b in specs):
        raise ValueError('Bond length must be less than half the shortest box side.')
    distance=_positive(min_distance,'min_distance')
    attempts=_integer(max_attempts,'max_attempts')
    rng=_rng(seed)
    if len(set(labels)) > 1: rng.shuffle(labels)
    bond_lengths = [] if model=='atomic' else [b.length for b in bonds_for_species(bond, labels)]
    labels = np.repeat(labels, 1 if model == 'atomic' else 2)
    mass_values = per_atom(parameter(masses, tuple(dict.fromkeys(labels.tolist())), 'masses'), labels)
    per_particle=1 if model=='atomic' else 2
    positions=np.empty((N*per_particle,3),dtype=float)
    for i in range(N):
        for _ in range(attempts):
            center=rng.uniform(0,lengths,size=(1,3))
            candidate=_atoms(center,lengths,model,rng,bond_lengths[i:i+1])
            dr=candidate[:,None,:]-positions[None,:i*per_particle,:]
            dr-=lengths*np.rint(dr/lengths)
            if np.all(np.sum(dr*dr,axis=2)>=distance**2):
                positions[i*per_particle:(i+1)*per_particle]=candidate
                break
        else:
            raise ValueError(f'Could only place {i} of {N} particles with min_distance={distance:g}. Lower rho, enlarge the box, reduce min_distance explicitly, or use FCC for a dense crystal.')
    return System(positions,lengths,model=model,bond=bond,species=labels,masses=mass_values)


def Plummer(*, N=256, scale_radius=1., mass=1., G=1., softening=.05, seed=87287):
    """An isolated, approximately stationary spherical cluster, including velocities.

    Sample rho(r) proportional to (1+r²/scale_radius²)^(-5/2) and the
    isotropic Plummer speed distribution. Remove COM motion and rescale speeds
    so 2K=-W for this finite sample and the chosen softened force. This is an
    initial condition, not a guarantee of equilibrium; relax before measuring.
    G and softening must match the simulation. All particles have the given mass.
    The box is a viewing frame of side 12*scale_radius, never a physical boundary.
    The distribution has an untruncated tail; particles outside the view remain
    in the simulation. Direct gravity supports at most 4096 particles.
    """
    from . import _core
    N=_integer(N, 'N', 2)
    if N>4096:raise ValueError('Plummer supports at most 4096 particles.')
    a=_positive(scale_radius, 'scale_radius')
    mass=_positive(mass, 'mass');G=_positive(G, 'G')
    softening=_positive(softening, 'softening', zero=True)
    rng=_rng(seed)
    r=a/np.sqrt(rng.uniform(np.finfo(float).eps, 1., N)**(-2/3)-1)
    direction=rng.normal(size=(N,3));direction/=np.linalg.norm(direction,axis=1)[:,None]
    x=r[:,None]*direction
    # q=v/v_escape has density proportional to q²(1-q²)^(7/2).
    q=np.empty(N);pending=np.arange(N)
    while len(pending):
        trial=rng.uniform(size=len(pending))
        accept=rng.uniform(0,.1,len(pending))<trial**2*(1-trial**2)**3.5
        q[pending[accept]]=trial[accept];pending=pending[~accept]
    speed=q*np.sqrt(2*G*N*mass/np.sqrt(r*r+a*a))
    direction=rng.normal(size=(N,3));direction/=np.linalg.norm(direction,axis=1)[:,None]
    v=speed[:,None]*direction;v-=v.mean(axis=0)
    x-=x.mean(axis=0);x+=6*a
    masses=np.full(N,mass);f=np.zeros_like(x)
    _,W,_,_=_core.advance_gravity(x,v,f,masses,0,.001,G,softening,0.)
    v*=np.sqrt(-W/np.sum(masses[:,None]*v*v))
    return System(x,12*a,velocities=v,masses=masses,boundary='open')


_ICE_BASIS = np.array([[0, 0, 0], [2/3, 1/3, .5], [0, 0, 3/8], [2/3, 1/3, 7/8]])
WATER_MASSES = {'O': 15.9994, 'H': 1.00794}


def _ice_cells(N):
    """Most nearly cubic nx*ny*nz = N/8 with at least two cells in every direction."""
    shape = np.array([1, np.sqrt(3), np.sqrt(8/3)])
    best = None
    if N % 8 == 0:
        cells = N//8
        for nx in range(2, cells+1):
            for ny in range(2, cells//nx+1):
                if cells % (nx*ny):
                    continue
                nz = cells//(nx*ny)
                if nz < 2:
                    continue
                lengths = shape*[nx, ny, nz]
                ratio = lengths.max()/lengths.min()
                if best is None or ratio < best[0]:
                    best = (ratio, (nx, ny, nz))
    return best


def IceIh(*, N, density=0.92, molecule='water', bond=None, masses=None, seed=87287):
    """Proton-disordered hexagonal ice in LAMMPS real units (Å, g/mol).

    N water molecules fill nx*ny*nz orthorhombic cells of 8 molecules, chosen as
    nearly cubic as possible; N=768 gives 6x4x4 cells and N=2880 gives 10x6x6.
    density is in g/cm^3 and fixes the lattice constant (ideal c/a = sqrt(8/3)).
    Hydrogens obey the ice rules (two donated and two accepted per oxygen), from
    a seeded Euler circuit on the oxygen network, and loop flips reduce the net
    dipole. Each molecule has the exact RigidWater geometry, oriented along its
    two donor bonds. Species are "O" and "H".
    """
    if molecule != 'water':
        raise ValueError('IceIh builds water; use molecule="water".')
    water = RigidWater() if bond is None else bond
    if not isinstance(water, RigidWater):
        raise TypeError('IceIh takes bond=RigidWater(...).')
    N = _integer(N, 'N', 16)
    if 3*N > 20000:
        raise ValueError('This course solver supports at most 20000 atoms.')
    density = _positive(density, 'density')
    masses = parameter(WATER_MASSES if masses is None else masses, ('O', 'H'), 'masses')
    best = _ice_cells(N)
    if best is None or best[0] > 1.3:
        nearby = sorted({8*nx*ny*nz for nx in range(2, 16) for ny in range(2, 12) for nz in range(2, 12)
                         if max(nx, ny*np.sqrt(3), nz*np.sqrt(8/3))/min(nx, ny*np.sqrt(3), nz*np.sqrt(8/3)) <= 1.3},
                        key=lambda count: abs(count-N))[:4]
        raise ValueError(f'Ice Ih needs N=8*nx*ny*nz in a nearly cubic box; N={N} is incompatible. '
                         f'Nearby counts: {", ".join(map(str, sorted(nearby)))}.')
    cells = np.array(best[1])
    molar_mass = masses['O'] + 2*masses['H']
    a = (8*molar_mass/(0.602214076*density*np.sqrt(8)))**(1/3)
    c = a*np.sqrt(8/3)
    unit = np.array([a, np.sqrt(3)*a, c])
    hexagonal = np.array([[a, 0, 0], [-a/2, np.sqrt(3)*a/2, 0], [0, 0, c]])
    base = np.array([((f + [0, j, 0]) @ hexagonal) % unit for j in range(2) for f in _ICE_BASIS])
    offsets = np.indices(cells).reshape(3, -1).T*unit
    oxygen = (offsets[:, None, :] + base[None, :, :]).reshape(-1, 3)
    box = cells*unit
    # Oxygen network: four tetrahedral neighbours at a*sqrt(3/8).
    dr = oxygen[:, None, :] - oxygen[None, :, :]
    dr -= box*np.rint(dr/box)
    distance = np.linalg.norm(dr, axis=2)
    i, j = np.nonzero(np.triu((distance > 0) & (distance < 1.1*a*np.sqrt(3/8))))
    edges = np.stack((i, j), axis=1)
    if np.any(np.bincount(edges.ravel(), minlength=len(oxygen)) != 4):
        raise RuntimeError('Ice lattice construction failed.')
    rng = _rng(seed)
    adjacency = [[] for _ in oxygen]
    for e, (u, v) in enumerate(edges):
        adjacency[u].append((v, e)); adjacency[v].append((u, e))
    for neighbours in adjacency:
        rng.shuffle(neighbours)
    # Euler circuit (Hierholzer): every vertex gets in-degree = out-degree = 2.
    donor = np.full(len(edges), -1)
    used = np.zeros(len(edges), dtype=bool)
    pointer = [0]*len(oxygen)
    stack = [0]
    while stack:
        vertex = stack[-1]
        while pointer[vertex] < 4 and used[adjacency[vertex][pointer[vertex]][1]]:
            pointer[vertex] += 1
        if pointer[vertex] == 4:
            stack.pop()
            continue
        neighbour, e = adjacency[vertex][pointer[vertex]]
        used[e] = True
        donor[e] = vertex
        stack.append(neighbour)
    if np.any(donor < 0):
        raise RuntimeError('Ice proton assignment failed.')
    # Unit bond vectors from each edge's first oxygen to its second.
    vectors = dr[edges[:, 1], edges[:, 0]]
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]
    def dipole_change(e):
        return -2*vectors[e] if donor[e] == edges[e, 0] else 2*vectors[e]
    def outgoing(vertex):
        return [e for _, e in adjacency[vertex] if donor[e] == vertex]
    dipole = sum(vectors[e] if donor[e] == edges[e, 0] else -vectors[e] for e in range(len(edges)))
    # Flipping a directed loop keeps the ice rules; accept flips that reduce |dipole|.
    for _ in range(40*len(oxygen)):
        start = int(rng.integers(len(oxygen)))
        path, seen, vertex = [], {start: 0}, start
        while True:
            e = outgoing(vertex)[int(rng.integers(2))]
            path.append(e)
            vertex = int(edges[e, 1] if edges[e, 0] == vertex else edges[e, 0])
            if vertex in seen:
                loop = path[seen[vertex]:]
                break
            seen[vertex] = len(path)
        change = sum(dipole_change(e) for e in loop)
        if np.linalg.norm(dipole + change) < np.linalg.norm(dipole):
            dipole = dipole + change
            for e in loop:
                donor[e] = edges[e, 1] if donor[e] == edges[e, 0] else edges[e, 0]
    half = np.radians(water.hoh_angle/2)
    positions = np.empty((3*N, 3))
    for vertex in range(N):
        bonds = [vectors[e] if edges[e, 0] == vertex else -vectors[e] for e in outgoing(vertex)]
        bisector = bonds[0] + bonds[1]; bisector /= np.linalg.norm(bisector)
        across = bonds[0] - bonds[1]; across /= np.linalg.norm(across)
        positions[3*vertex] = oxygen[vertex]
        positions[3*vertex+1] = oxygen[vertex] + water.oh_length*(np.cos(half)*bisector + np.sin(half)*across)
        positions[3*vertex+2] = oxygen[vertex] + water.oh_length*(np.cos(half)*bisector - np.sin(half)*across)
    labels = np.tile(['O', 'H', 'H'], N)
    return System(positions, box, molecule='water', bond=water, species=labels,
                  masses=per_atom(masses, labels), units='real')


def _random_water(*, N, rho, box, min_distance, max_attempts, model, molecule, bond, species, masses, seed):
    if species is not None:
        raise ValueError('Water species are always "O" and "H"; omit species.')
    model, water = resolve_model(model, molecule, bond)
    N = _integer(N, 'N', 2)
    if 3*N > 20000:
        raise ValueError('This course solver supports at most 20000 atoms.')
    if rho is not None and box is not None:
        raise ValueError('Specify rho or box, not both; N fixes the remaining quantity.')
    if box is None:
        lengths = np.full(3, (N/_positive(.0334 if rho is None else rho, 'rho'))**(1/3))
    else:
        lengths = np.broadcast_to(np.asarray(box, dtype=float), (3,)).copy()
        if not np.isfinite(lengths).all() or np.any(lengths <= 0):
            raise ValueError('box must be a positive side length or three positive lengths.')
    masses = parameter(WATER_MASSES if masses is None else masses, ('O', 'H'), 'masses')
    distance = _positive(min_distance, 'min_distance')
    attempts = _integer(max_attempts, 'max_attempts')
    rng = _rng(seed)
    positions = np.empty((3*N, 3))
    for i in range(N):
        for _ in range(attempts):
            candidate = _water_molecules(rng.uniform(0, lengths, size=(1, 3)), water, rng, masses)
            dr = candidate[:, None, :] - positions[None, :3*i, :]
            dr -= lengths*np.rint(dr/lengths)
            if np.all(np.sum(dr*dr, axis=2) >= distance**2):
                positions[3*i:3*i+3] = candidate
                break
        else:
            raise ValueError(f'Could only place {i} of {N} water molecules with min_distance={distance:g} Å. '
                             'Enlarge the box, reduce min_distance explicitly, or start from IceIh.')
    labels = np.tile(['O', 'H', 'H'], N)
    return System(positions % lengths, lengths, molecule='water', bond=water, species=labels,
                  masses=per_atom(masses, labels), units='real')
