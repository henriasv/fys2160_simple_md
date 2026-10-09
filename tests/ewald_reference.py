"""Direct Ewald sum for point charges in an orthorhombic periodic box (tests only).

E = C [1/2 sum_{i,j,n}' q_i q_j erfc(a r)/r + 2 pi/V sum_k exp(-k^2/4a^2)/k^2 |S(k)|^2
       - a/sqrt(pi) sum q^2], with S(k) = sum_j q_j exp(i k.r_j).
"""
from math import erfc, pi, sqrt
import numpy as np


def ewald(positions, charges, box, *, alpha=0.35, real_cutoff=None, kmax=None, coulomb=1.0, exclusions=()):
    """Energy, forces and virial trace. exclusions lists (i, j) pairs whose bare
    Coulomb interaction (nearest image) is removed, as for rigid molecules."""
    x = np.asarray(positions, float)
    q = np.asarray(charges, float)
    box = np.asarray(box, float)
    n = len(x)
    real_cutoff = real_cutoff or 6.0/alpha
    images = [np.array([i, j, k]) for i in range(-int(np.ceil(real_cutoff/box[0])), int(np.ceil(real_cutoff/box[0]))+1)
              for j in range(-int(np.ceil(real_cutoff/box[1])), int(np.ceil(real_cutoff/box[1]))+1)
              for k in range(-int(np.ceil(real_cutoff/box[2])), int(np.ceil(real_cutoff/box[2]))+1)]
    energy = 0.0
    forces = np.zeros_like(x)
    virial = 0.0
    for i in range(n):
        for j in range(n):
            for image in images:
                if i == j and not image.any():
                    continue
                dr = x[i] - x[j] + image*box
                r = np.linalg.norm(dr)
                if r > real_cutoff:
                    continue
                energy += 0.5*q[i]*q[j]*erfc(alpha*r)/r
                magnitude = q[i]*q[j]*(erfc(alpha*r)/r + 2*alpha/sqrt(pi)*np.exp(-alpha*alpha*r*r))/(r*r)
                forces[i] += magnitude*dr
                virial += 0.5*magnitude*r*r
    volume = np.prod(box)
    kmax = kmax or int(np.ceil(3.2*alpha*box.max()/pi))+2
    for a in range(-kmax, kmax+1):
        for b in range(-kmax, kmax+1):
            for c in range(-kmax, kmax+1):
                if a == b == c == 0:
                    continue
                k = 2*pi*np.array([a, b, c])/box
                k2 = k @ k
                weight = np.exp(-k2/(4*alpha*alpha))/k2
                if weight < 1e-18:
                    continue
                phase = x @ k
                S = np.sum(q*np.exp(1j*phase))
                energy += 2*pi/volume*weight*abs(S)**2
                virial += 2*pi/volume*weight*abs(S)**2*(1 - k2/(2*alpha*alpha))
                forces += (4*pi/volume*weight*q*np.imag(np.conj(S)*np.exp(1j*phase)))[:, None]*k
    energy -= alpha/sqrt(pi)*np.sum(q*q)
    for i, j in exclusions:
        dr = x[i] - x[j]
        dr -= box*np.rint(dr/box)
        r = np.linalg.norm(dr)
        energy -= q[i]*q[j]/r
        magnitude = q[i]*q[j]/r**3
        forces[i] -= magnitude*dr
        forces[j] += magnitude*dr
        virial -= magnitude*r*r
    return coulomb*energy, coulomb*forces, coulomb*virial
