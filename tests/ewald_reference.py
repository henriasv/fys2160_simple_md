"""Direct Ewald sum for point charges in a periodic box, rectangular or tilted (tests only).

E = C [1/2 sum_{i,j,n}' q_i q_j erfc(a r)/r + 2 pi/V sum_k exp(-k^2/4a^2)/k^2 |S(k)|^2
       - a/sqrt(pi) sum q^2], with S(k) = sum_j q_j exp(i k.r_j).
"""
from math import erfc, pi, sqrt
import numpy as np


def ewald(positions, charges, box, *, alpha=0.35, real_cutoff=None, kmax=None, coulomb=1.0, exclusions=(), tilt=None,
          tensor=False):
    """Energy, forces and virial trace, and with tensor=True the virial tensor (xx, yy,
    zz, xy, xz, yz) as well. exclusions lists (i, j) pairs whose bare Coulomb
    interaction (nearest image) is removed, as for rigid molecules. tilt = (xy, xz,
    yz) makes the box tilted, with edges (lx, 0, 0), (xy, ly, 0) and (xz, yz, lz)."""
    x = np.asarray(positions, float)
    q = np.asarray(charges, float)
    box = np.asarray(box, float)
    xy, xz, yz = (0., 0., 0.) if tilt is None else tilt
    edges = np.array([[box[0], 0, 0], [xy, box[1], 0], [xz, yz, box[2]]])  # rows a, b, c
    reciprocal = 2*pi*np.linalg.inv(edges).T  # rows: k for unit wave numbers along a, b, c
    width = 1/np.linalg.norm(np.linalg.inv(edges), axis=0)  # distances between opposite faces
    n = len(x)
    real_cutoff = real_cutoff or 6.0/alpha
    reach = [int(np.ceil(real_cutoff/w))+1 for w in width]
    images = [np.array([i, j, k])@edges for i in range(-reach[0], reach[0]+1)
              for j in range(-reach[1], reach[1]+1) for k in range(-reach[2], reach[2]+1)]
    energy = 0.0
    forces = np.zeros_like(x)
    virial = 0.0
    W = np.zeros((3, 3))
    for i in range(n):
        for j in range(n):
            for image in images:
                if i == j and not image.any():
                    continue
                dr = x[i] - x[j] + image
                r = np.linalg.norm(dr)
                if r > real_cutoff:
                    continue
                energy += 0.5*q[i]*q[j]*erfc(alpha*r)/r
                magnitude = q[i]*q[j]*(erfc(alpha*r)/r + 2*alpha/sqrt(pi)*np.exp(-alpha*alpha*r*r))/(r*r)
                forces[i] += magnitude*dr
                virial += 0.5*magnitude*r*r
                W += 0.5*magnitude*np.outer(dr, dr)
    volume = np.prod(box)
    # |k| below the Gaussian's reach has wave numbers up to |k| |edge|/(2 pi) along each edge.
    kmax = kmax or int(np.ceil(3.2*alpha*np.linalg.norm(edges, axis=1).max()/pi))+2
    for a in range(-kmax, kmax+1):
        for b in range(-kmax, kmax+1):
            for c in range(-kmax, kmax+1):
                if a == b == c == 0:
                    continue
                k = np.array([a, b, c])@reciprocal
                k2 = k @ k
                weight = np.exp(-k2/(4*alpha*alpha))/k2
                if weight < 1e-18:
                    continue
                phase = x @ k
                S = np.sum(q*np.exp(1j*phase))
                term = 2*pi/volume*weight*abs(S)**2
                energy += term
                virial += term*(1 - k2/(2*alpha*alpha))
                W += term*(np.eye(3) - 2*np.outer(k, k)*(1/k2 + 1/(4*alpha*alpha)))
                forces += (4*pi/volume*weight*q*np.imag(np.conj(S)*np.exp(1j*phase)))[:, None]*k
    energy -= alpha/sqrt(pi)*np.sum(q*q)
    for i, j in exclusions:
        dr = x[i] - x[j]
        dr -= np.rint(dr@np.linalg.inv(edges))@edges
        r = np.linalg.norm(dr)
        energy -= q[i]*q[j]/r
        magnitude = q[i]*q[j]/r**3
        forces[i] -= magnitude*dr
        forces[j] += magnitude*dr
        virial -= magnitude*r*r
    if tensor:
        return coulomb*energy, coulomb*forces, coulomb*virial, coulomb*W[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]
    return coulomb*energy, coulomb*forces, coulomb*virial
