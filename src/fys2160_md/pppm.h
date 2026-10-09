/* Particle-particle particle-mesh (PPPM) electrostatics with ik differentiation
   and the Hockney-Eastwood optimal influence function, as in LAMMPS kspace pppm. */
#ifndef FYS2160_PPPM_H
#define FYS2160_PPPM_H
#include "fft.h"

typedef struct {
    int order, n[3];
    double g;              /* Ewald splitting parameter */
    double box[3];         /* box for which greens and k vectors are valid */
    double *greens;        /* optimal influence function, one value per mesh point */
    double *k[3];          /* wave numbers along each axis */
    double *rho, *work;    /* complex meshes (interleaved re, im) */
    double *weights;       /* order weights per axis per site */
    int *start;            /* first stencil index per axis per site */
    int sites;
    FFT3D fft;
} PPPM;

/* LAMMPS-style choice of g and mesh from relative accuracy. qqr2e*accuracy is the
   absolute force tolerance. Returns 0 on success, writes estimated relative accuracy. */
int pppm_select(const double box[3], int atoms, double qsqsum, double qqr2e, double cutoff,
                double accuracy, int order, int mesh[3], double *g, double *estimate);
int pppm_init(PPPM *p, const int mesh[3], double g, int order, int sites);
void pppm_free(PPPM *p);
/* Adds k-space forces to f (sites x 3); returns energy (including self and
   neutralizing background) and the k-space virial trace. */
int pppm_compute(PPPM *p, const double box[3], int sites, const double *x, const double *q,
                 double qqr2e, double *f, double *energy, double *virial);

#endif
