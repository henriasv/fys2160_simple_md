"""Periodic geometry constructors, independent of the simulation runner."""
import numpy as np
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


def _atoms(centers, lengths, model, rng):
    if model=='atomic':return centers
    directions=rng.normal(size=(len(centers),3))
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    return np.stack((centers-.35*directions,centers+.35*directions),axis=1).reshape(-1,3)%lengths


def FCC(*, N=500, rho=None, box=None, model='atomic', seed=87287, **settings):
    """Perfect periodic FCC sites at rho+N or box+N; the default box is cubic.

    A scalar box is a cube side length; a triple gives rectangular box lengths.
    N counts atoms, or molecular centres for a diatomic model. Cubes require
    N=4*n**3. Rectangular boxes must fit integer conventional cells of equal size.
    """
    N,lengths=_box(N,rho,box,model)
    lattice_constant=(4*np.prod(lengths)/N)**(1/3)
    ratios=lengths/lattice_constant
    cells=np.rint(ratios).astype(int)
    if np.any(cells<1) or 4*np.prod(cells)!=N or not np.allclose(ratios,cells,rtol=1e-10,atol=1e-10):
        if np.allclose(lengths,lengths[0],rtol=1e-10,atol=1e-10):
            lower=max(1,int(np.floor(np.cbrt(N/4))))
            raise ValueError(f'A cubic FCC crystal needs N=4*n^3; N={N} is incompatible. Nearby counts: {4*lower**3} or {4*(lower+1)**3}.')
        raise ValueError('A rectangular FCC box must contain integer numbers of equal-sized conventional cells, with N=4*nx*ny*nz. Adjust N or the box aspect ratio.')
    rng=_rng(seed)
    bases=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
    sites=(np.indices(tuple(cells)).reshape(3,-1).T[:,None,:]+bases).reshape(-1,3)
    centers=(sites+.25)*lattice_constant
    return System(_atoms(centers,lengths,model,rng),lengths,model=model,seed=seed,_initial_rng=rng,**settings)


def Random(*, N=500, rho=None, box=None, min_distance=.9, max_attempts=1000,
           model='atomic', seed=87287, **settings):
    """Self-avoiding random placement with periodic minimum-image distances.

    min_distance is the minimum atom separation in sigma units (excluding the
    two partners within a molecule). Failed packing raises an error rather
    than reducing this distance. max_attempts bounds trials for each particle.
    """
    N,lengths=_box(N,rho,box,model)
    distance=_positive(min_distance,'min_distance')
    attempts=_integer(max_attempts,'max_attempts')
    rng=_rng(seed)
    per_particle=1 if model=='atomic' else 2
    positions=np.empty((N*per_particle,3),dtype=float)
    for i in range(N):
        for _ in range(attempts):
            center=rng.uniform(0,lengths,size=(1,3))
            candidate=_atoms(center,lengths,model,rng)
            dr=candidate[:,None,:]-positions[None,:i*per_particle,:]
            dr-=lengths*np.rint(dr/lengths)
            if np.all(np.sum(dr*dr,axis=2)>=distance**2):
                positions[i*per_particle:(i+1)*per_particle]=candidate
                break
        else:
            raise ValueError(f'Could only place {i} of {N} particles with min_distance={distance:g}. Lower rho, enlarge the box, reduce min_distance explicitly, or use FCC for a dense crystal.')
    return System(positions,lengths,model=model,seed=seed,_initial_rng=rng,**settings)
