/* Charged atoms and rigid three-site water with optional TIP4P-type M site.
   Units are supplied by the caller: mvv2e converts mass*velocity^2 to energy,
   qqr2e is the Coulomb constant, and temperature is passed as kB*T. */
#include "molecular.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "pppm.h"

static double uniform01(uint64_t *state) {
    uint64_t x = *state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    *state = x;
    return ((double)((x * UINT64_C(2685821657736338717)) >> 11) + 0.5)
           * (1.0 / 9007199254740992.0);
}
static double gaussian(uint64_t *state) {
    double u = uniform01(state), v = uniform01(state);
    return sqrt(-2.0 * log(u)) * cos(6.2831853071795864769 * v);
}
static double minimum_image(double x, double length) {
    return x - length * nearbyint(x / length);
}

/* Constraint c joins atoms PAIR[c][0] and PAIR[c][1] of a water molecule. */
static const int PAIR[3][2] = {{0, 1}, {0, 2}, {1, 2}};

typedef struct {
    int n, size, molecules, types, coulomb;
    const double *mass, *charge, *epsilon, *sigma;
    const int *type;
    double *shift;
    double cutoff, coulomb_cutoff, extent, alpha, qqr2e, mvv2e, g;
    double constraint[3];
} Model;

typedef struct {
    int *pairs;
    size_t size, capacity;
    double *reference, box[3], radius;
} MoleculeList;

typedef struct { double lj, coulomb, kspace, virial, molecular_virial; } Terms;
/* CPU seconds per part, reported for benchmarks. */
typedef struct { double list, pair, kspace, total; } Timings;
static double seconds(void) { return (double)clock()/CLOCKS_PER_SEC; }

/* Abramowitz-Stegun 7.1.26 erfc (|error| < 1.5e-7, as LAMMPS without tables) and
   its exact derivative, so real-space forces are gradients of the energies used. */
static double erfc_approx(double x, double *derivative) {
    const double p = 0.3275911, a1 = 0.254829592, a2 = -0.284496736,
                 a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
    double t = 1/(1 + p*x), e = exp(-x*x);
    double poly = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))));
    double slope = a1 + t*(2*a2 + t*(3*a3 + t*(4*a4 + t*5*a5)));
    *derivative = e*(-p*t*t*slope - 2*x*poly);
    return poly*e;
}

/* Real-space Ewald pair, A*(erfc(g r) - excluded)/r. coefficient = -(dE/dr)/r. */
static double coulomb_term(double r2, double A, double g, double excluded, double *coefficient) {
    double r = sqrt(r2), slope, complement = erfc_approx(g*r, &slope);
    double dEdr = A*(g*slope/r - (complement - excluded)/r2);
    *coefficient = -dEdr/r;
    return A*(complement - excluded)/r;
}

static void sites_and_offsets(const Model *m, const double *x, const double *box, double *site, double *offset) {
    if (m->size == 1) {
        memcpy(site, x, (size_t)m->n*3*sizeof(double));
        memset(offset, 0, (size_t)m->n*3*sizeof(double));
        return;
    }
    for (int I = 0; I < m->molecules; ++I) {
        int o = 3*I;
        for (int d = 0; d < 3; ++d) {
            double origin = x[3*o+d];
            double h1 = minimum_image(x[3*(o+1)+d] - origin, box[d]);
            double h2 = minimum_image(x[3*(o+2)+d] - origin, box[d]);
            offset[3*o+d] = 0; offset[3*(o+1)+d] = h1; offset[3*(o+2)+d] = h2;
            site[3*o+d] = origin + m->alpha*(h1 + h2);
            site[3*(o+1)+d] = origin + h1;
            site[3*(o+2)+d] = origin + h2;
        }
    }
}

static int build_list(MoleculeList *nb, const Model *m, const double *x, const double *u, const double *box) {
    nb->size = 0;
    double r2max = nb->radius*nb->radius;
    for (int I = 0; I < m->molecules; ++I) for (int J = I+1; J < m->molecules; ++J) {
        int a = m->size*I, b = m->size*J;
        double r2 = 0;
        for (int d = 0; d < 3; ++d) {
            double dx = minimum_image(x[3*a+d] - x[3*b+d], box[d]);
            r2 += dx*dx;
        }
        if (r2 >= r2max) continue;
        if (nb->size == nb->capacity) {
            size_t capacity = nb->capacity ? 2*nb->capacity : 4096;
            int *pairs = realloc(nb->pairs, 2*capacity*sizeof(int));
            if (!pairs) { PyErr_NoMemory(); return -1; }
            nb->pairs = pairs; nb->capacity = capacity;
        }
        nb->pairs[2*nb->size] = I; nb->pairs[2*nb->size+1] = J; ++nb->size;
    }
    memcpy(nb->reference, u, (size_t)m->n*3*sizeof(double));
    memcpy(nb->box, box, 3*sizeof(double));
    return 0;
}

static int needs_rebuild(const MoleculeList *nb, int n, const double *u, const double *box, double skin) {
    double scale = box[0]/nb->box[0];
    for (int d = 1; d < 3; ++d) scale = fmin(scale, box[d]/nb->box[d]);
    double allowance = 0.5*(scale*nb->radius - (nb->radius - skin));
    if (allowance <= 0) return 1;
    for (int i = 0; i < n; ++i) {
        double displacement2 = 0;
        for (int d = 0; d < 3; ++d) {
            double dx = u[3*i+d] - nb->reference[3*i+d]*box[d]/nb->box[d];
            displacement2 += dx*dx;
        }
        if (displacement2 > allowance*allowance) return 1;
    }
    return 0;
}

static int forces(const Model *m, const MoleculeList *nb, PPPM *pppm, const double *x, const double *box,
                  double *site, double *offset, double *site_force, double *f, Terms *t, Timings *clock_) {
    double start = seconds();
    const int n = m->n, size = m->size;
    memset(f, 0, (size_t)n*3*sizeof(double));
    memset(site_force, 0, (size_t)n*3*sizeof(double));
    memset(t, 0, sizeof *t);
    sites_and_offsets(m, x, box, site, offset);
    const double rc2 = m->cutoff*m->cutoff, rcc2 = m->coulomb_cutoff*m->coulomb_cutoff;
    for (size_t p = 0; p < nb->size; ++p) {
        int a0 = size*nb->pairs[2*p], b0 = size*nb->pairs[2*p+1];
        double shift[3];
        for (int d = 0; d < 3; ++d) shift[d] = box[d]*nearbyint((x[3*a0+d] - x[3*b0+d])/box[d]);
        for (int a = a0; a < a0+size; ++a) for (int b = b0; b < b0+size; ++b) {
            int pair = m->type[a]*m->types + m->type[b];
            double epsilon = m->epsilon[pair];
            if (epsilon == 0) continue;
            double dr[3], r2 = 0;
            for (int d = 0; d < 3; ++d) {
                dr[d] = x[3*a0+d] + offset[3*a+d] - x[3*b0+d] - offset[3*b+d] - shift[d];
                r2 += dr[d]*dr[d];
            }
            if (r2 >= rc2) continue;
            if (r2 < 1e-16 || !isfinite(r2)) {
                PyErr_SetString(PyExc_ValueError, "Overlapping atoms or unstable trajectory; check positions and timestep.");
                return -1;
            }
            double s2 = m->sigma[pair]*m->sigma[pair]/r2, inv6 = s2*s2*s2;
            double coefficient = epsilon*24*inv6*(2*inv6 - 1)/r2;
            t->lj += epsilon*4*inv6*(inv6 - 1) - m->shift[pair];
            t->virial += coefficient*r2;
            for (int d = 0; d < 3; ++d) { f[3*a+d] += coefficient*dr[d]; f[3*b+d] -= coefficient*dr[d]; }
        }
        if (!m->coulomb) continue;
        for (int s = a0; s < a0+size; ++s) {
            if (m->charge[s] == 0) continue;
            for (int q = b0; q < b0+size; ++q) {
                if (m->charge[q] == 0) continue;
                double dr[3], r2 = 0;
                for (int d = 0; d < 3; ++d) { dr[d] = site[3*s+d] - site[3*q+d] - shift[d]; r2 += dr[d]*dr[d]; }
                if (r2 >= rcc2) continue;
                if (r2 < 1e-16 || !isfinite(r2)) {
                    PyErr_SetString(PyExc_ValueError, "Overlapping charges or unstable trajectory; check positions and timestep.");
                    return -1;
                }
                double coefficient;
                t->coulomb += coulomb_term(r2, m->qqr2e*m->charge[s]*m->charge[q], m->g, 0, &coefficient);
                t->virial += coefficient*r2;
                for (int d = 0; d < 3; ++d) { site_force[3*s+d] += coefficient*dr[d]; site_force[3*q+d] -= coefficient*dr[d]; }
            }
        }
    }
    if (m->coulomb) {
        /* Remove the intramolecular part that the mesh sum includes. */
        if (size == 3) for (int I = 0; I < m->molecules; ++I) for (int c = 0; c < 3; ++c) {
            int s = 3*I + PAIR[c][0], q = 3*I + PAIR[c][1];
            if (m->charge[s] == 0 || m->charge[q] == 0) continue;
            double dr[3], r2 = 0, coefficient;
            for (int d = 0; d < 3; ++d) { dr[d] = site[3*s+d] - site[3*q+d]; r2 += dr[d]*dr[d]; }
            t->coulomb += coulomb_term(r2, m->qqr2e*m->charge[s]*m->charge[q], m->g, 1, &coefficient);
            t->virial += coefficient*r2;
            for (int d = 0; d < 3; ++d) { site_force[3*s+d] += coefficient*dr[d]; site_force[3*q+d] -= coefficient*dr[d]; }
        }
        double kspace_virial, middle = seconds();
        clock_->pair += middle - start;
        if (pppm_compute(pppm, box, n, site, m->charge, m->qqr2e, site_force, &t->kspace, &kspace_virial) < 0) {
            PyErr_SetString(PyExc_FloatingPointError, "PPPM failed; check the box, charges and timestep.");
            return -1;
        }
        start = seconds();
        clock_->kspace += start - middle;
        t->virial += kspace_virial;
        if (size == 1) for (int k = 0; k < 3*n; ++k) f[k] += site_force[k];
        else for (int I = 0; I < m->molecules; ++I) for (int d = 0; d < 3; ++d) {
            int o = 3*I;
            double fm = site_force[3*o+d];
            f[3*o+d] += (1 - 2*m->alpha)*fm;
            f[3*(o+1)+d] += m->alpha*fm + site_force[3*(o+1)+d];
            f[3*(o+2)+d] += m->alpha*fm + site_force[3*(o+2)+d];
        }
    }
    t->molecular_virial = t->virial;
    clock_->pair -= start;
    if (size == 3) for (int I = 0; I < m->molecules; ++I) {
        int o = 3*I;
        double total = m->mass[o] + m->mass[o+1] + m->mass[o+2];
        for (int d = 0; d < 3; ++d) {
            double centre = (m->mass[o+1]*offset[3*(o+1)+d] + m->mass[o+2]*offset[3*(o+2)+d])/total;
            for (int a = o; a < o+3; ++a) t->molecular_virial -= (offset[3*a+d] - centre)*f[3*a+d];
        }
    }
    clock_->pair += seconds();
    if (!isfinite(t->lj + t->coulomb + t->kspace + t->virial + t->molecular_virial)) {
        PyErr_SetString(PyExc_FloatingPointError, "Non-finite energy; reduce timestep or relax the initial positions.");
        return -1;
    }
    return 0;
}

static double kinetic(const Model *m, const double *v) {
    double K = 0;
    for (int i = 0; i < m->n; ++i) for (int d = 0; d < 3; ++d) K += 0.5*m->mass[i]*v[3*i+d]*v[3*i+d];
    return K*m->mvv2e;
}

static double centre_kinetic(const Model *m, const double *v) {
    double K = 0;
    for (int I = 0; I < m->molecules; ++I) {
        double total = 0, p[3] = {0, 0, 0};
        for (int a = m->size*I; a < m->size*(I+1); ++a) {
            total += m->mass[a];
            for (int d = 0; d < 3; ++d) p[d] += m->mass[a]*v[3*a+d];
        }
        K += 0.5*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/total;
    }
    return K*m->mvv2e;
}

static void kick(const Model *m, double *v, const double *f, double dt) {
    for (int i = 0; i < m->n; ++i) for (int d = 0; d < 3; ++d) v[3*i+d] += dt*f[3*i+d]/(m->mass[i]*m->mvv2e);
}

static int drift(const Model *m, double *x, double *u, const double *v, const double *box, double duration) {
    for (int i = 0; i < m->n; ++i) for (int d = 0; d < 3; ++d) {
        int k = 3*i+d;
        double dx = duration*v[k];
        u[k] += dx; x[k] += dx;
        if (!isfinite(u[k]) || !isfinite(x[k])) {
            PyErr_SetString(PyExc_FloatingPointError, "Non-finite position; simulation is unstable.");
            return -1;
        }
        x[k] -= box[d]*floor(x[k]/box[d]);
    }
    return 0;
}

static void thermostat(const Model *m, double *v, double T, double gamma, double dt, uint64_t *rng) {
    double c = exp(-gamma*dt), momentum[3] = {0, 0, 0}, total = 0;
    for (int i = 0; i < m->n; ++i) {
        double sigma = sqrt(-expm1(-2*gamma*dt)*T/(m->mass[i]*m->mvv2e));
        total += m->mass[i];
        for (int d = 0; d < 3; ++d) {
            v[3*i+d] = c*v[3*i+d] + sigma*gaussian(rng);
            momentum[d] += m->mass[i]*v[3*i+d];
        }
    }
    for (int i = 0; i < m->n; ++i) for (int d = 0; d < 3; ++d) v[3*i+d] -= momentum[d]/total;
}

/* Solve a 3x3 system by Gaussian elimination with partial pivoting. */
static int solve3(double A[3][3], double b[3], double x[3]) {
    int order[3] = {0, 1, 2};
    for (int c = 0; c < 3; ++c) {
        int best = c;
        for (int r = c+1; r < 3; ++r) if (fabs(A[order[r]][c]) > fabs(A[order[best]][c])) best = r;
        int tmp = order[c]; order[c] = order[best]; order[best] = tmp;
        double pivot = A[order[c]][c];
        if (fabs(pivot) < 1e-300) return -1;
        for (int r = c+1; r < 3; ++r) {
            double factor = A[order[r]][c]/pivot;
            for (int k = c; k < 3; ++k) A[order[r]][k] -= factor*A[order[c]][k];
            b[order[r]] -= factor*b[order[c]];
        }
    }
    for (int c = 2; c >= 0; --c) {
        double s = b[order[c]];
        for (int k = c+1; k < 3; ++k) s -= A[order[c]][k]*x[k];
        x[c] = s/A[order[c]][c];
    }
    return 0;
}

/* g[c][a]: +1 if atom a is the first atom of constraint c, -1 if second. */
static double sign_of(int c, int a) { return PAIR[c][0] == a ? 1 : PAIR[c][1] == a ? -1 : 0; }

/* K[c][k] = g(k, first of c)/m_first - g(k, second of c)/m_second. */
static void coupling(const double *mass, double K[3][3]) {
    for (int c = 0; c < 3; ++c) for (int k = 0; k < 3; ++k)
        K[c][k] = sign_of(k, PAIR[c][0])/mass[PAIR[c][0]] - sign_of(k, PAIR[c][1])/mass[PAIR[c][1]];
}

static void bond_vectors(const double *x, const double *box, int o, double r[3][3]) {
    for (int c = 0; c < 3; ++c) for (int d = 0; d < 3; ++d)
        r[c][d] = minimum_image(x[3*(o+PAIR[c][0])+d] - x[3*(o+PAIR[c][1])+d], box[d]);
}

/* Remove relative velocity components along all three constraints. */
static int constrain_velocities(const Model *m, double *v, const double *x, const double *box) {
    for (int I = 0; I < m->molecules; ++I) {
        int o = 3*I;
        double r[3][3], K[3][3], A[3][3], b[3], mu[3];
        bond_vectors(x, box, o, r);
        coupling(m->mass + o, K);
        for (int c = 0; c < 3; ++c) {
            b[c] = 0;
            for (int d = 0; d < 3; ++d) b[c] -= r[c][d]*(v[3*(o+PAIR[c][0])+d] - v[3*(o+PAIR[c][1])+d]);
            for (int k = 0; k < 3; ++k) {
                double dot = 0;
                for (int d = 0; d < 3; ++d) dot += r[c][d]*r[k][d];
                A[c][k] = dot*K[c][k];
            }
        }
        if (solve3(A, b, mu) < 0) { PyErr_SetString(PyExc_ValueError, "Degenerate water geometry."); return -1; }
        for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k) {
            double s = sign_of(k, a)*mu[k]/m->mass[o+a];
            if (s != 0) for (int d = 0; d < 3; ++d) v[3*(o+a)+d] += s*r[k][d];
        }
    }
    return 0;
}

/* RATTLE positions: Newton iteration for multipliers so that drifting with the
   corrected velocities restores all three distances. */
static int constrained_drift(const Model *m, double *x, double *u, double *v, const double *box, double dt) {
    if (m->size == 3) for (int I = 0; I < m->molecules; ++I) {
        int o = 3*I;
        double old[3][3], trial[3][3], K[3][3], lambda[3] = {0, 0, 0};
        bond_vectors(x, box, o, old);
        coupling(m->mass + o, K);
        for (int c = 0; c < 3; ++c) for (int d = 0; d < 3; ++d)
            trial[c][d] = old[c][d] + dt*(v[3*(o+PAIR[c][0])+d] - v[3*(o+PAIR[c][1])+d]);
        int converged = 0;
        for (int iteration = 0; iteration < 100 && !converged; ++iteration) {
            double now[3][3], F[3], A[3][3], step[3];
            converged = 1;
            for (int c = 0; c < 3; ++c) {
                double length2 = 0;
                for (int d = 0; d < 3; ++d) {
                    now[c][d] = trial[c][d];
                    for (int k = 0; k < 3; ++k) now[c][d] += dt*K[c][k]*lambda[k]*old[k][d];
                    length2 += now[c][d]*now[c][d];
                }
                F[c] = length2 - m->constraint[c]*m->constraint[c];
                if (fabs(F[c]) > 1e-13*m->constraint[c]*m->constraint[c]) converged = 0;
            }
            if (converged) break;
            for (int c = 0; c < 3; ++c) {
                for (int k = 0; k < 3; ++k) {
                    double dot = 0;
                    for (int d = 0; d < 3; ++d) dot += now[c][d]*old[k][d];
                    A[c][k] = 2*dt*K[c][k]*dot;
                }
                F[c] = -F[c];
            }
            if (solve3(A, F, step) < 0) break;
            for (int k = 0; k < 3; ++k) lambda[k] += step[k];
        }
        if (!converged) { PyErr_SetString(PyExc_ValueError, "Rigid-water constraint failed; reduce timestep."); return -1; }
        for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k) {
            double s = sign_of(k, a)*lambda[k]/m->mass[o+a];
            if (s != 0) for (int d = 0; d < 3; ++d) v[3*(o+a)+d] += s*old[k][d];
        }
    }
    return drift(m, x, u, v, box, dt);
}

/* Move whole molecules affinely with their centres of mass (box lengths change by
   factor along each axis). Internal geometry is unchanged. */
static void scale_molecules(const Model *m, double *x, double *u, const double *box, const double factor[3]) {
    double *arrays[2] = {x, u};
    for (int I = 0; I < m->molecules; ++I) {
        int o = m->size*I;
        double total = 0;
        for (int a = o; a < o+m->size; ++a) total += m->mass[a];
        for (int which = 0; which < 2; ++which) {
            double *r = arrays[which];
            for (int d = 0; d < 3; ++d) {
                double origin = r[3*o+d], centre = 0, local[3];
                for (int a = 0; a < m->size; ++a) {
                    local[a] = a ? minimum_image(r[3*(o+a)+d] - origin, box[d]) : 0;
                    centre += m->mass[o+a]*local[a];
                }
                centre /= total;
                for (int a = 0; a < m->size; ++a) r[3*(o+a)+d] += (factor[d] - 1)*(r[3*(o+a)+d] - local[a] + centre);
            }
        }
    }
    for (int i = 0; i < m->n; ++i) for (int d = 0; d < 3; ++d) {
        double length = box[d]*factor[d];
        x[3*i+d] -= length*floor(x[3*i+d]/length);
    }
}

static void scale_centre_velocities(const Model *m, double *v, double factor) {
    for (int I = 0; I < m->molecules; ++I) {
        int o = m->size*I;
        double total = 0, p[3] = {0, 0, 0};
        for (int a = o; a < o+m->size; ++a) {
            total += m->mass[a];
            for (int d = 0; d < 3; ++d) p[d] += m->mass[a]*v[3*a+d];
        }
        for (int a = o; a < o+m->size; ++a) for (int d = 0; d < 3; ++d) v[3*a+d] += (factor - 1)*p[d]/total;
    }
}

/* state: mesh x,y,z, g, reference box x,y,z, estimated relative accuracy. */
static int prepare_pppm(Model *m, PPPM *pppm, int *ready, double *state, const double *box, double accuracy, int order) {
    int select = !(state[0] >= 1);
    for (int d = 0; d < 3 && !select; ++d) if (fabs(box[d]/state[4+d] - 1) > 0.05) select = 1;
    if (select) {
        double qsqsum = 0, g, estimate;
        int mesh[3];
        for (int i = 0; i < m->n; ++i) qsqsum += m->charge[i]*m->charge[i];
        int status = pppm_select(box, m->n, qsqsum, m->qqr2e, m->coulomb_cutoff, accuracy, order, mesh, &g, &estimate);
        if (status < 0) {
            PyErr_SetString(PyExc_ValueError, status == -2 ? "PPPM mesh would exceed 4096 points per side; increase coulomb_cutoff or kspace_accuracy." : "Could not choose PPPM parameters for this accuracy.");
            return -1;
        }
        for (int d = 0; d < 3; ++d) { state[d] = mesh[d]; state[4+d] = box[d]; }
        state[3] = g; state[7] = estimate;
    }
    int mesh[3] = {(int)state[0], (int)state[1], (int)state[2]};
    if (*ready && (pppm->n[0] != mesh[0] || pppm->n[1] != mesh[1] || pppm->n[2] != mesh[2] || pppm->g != state[3])) {
        pppm_free(pppm); *ready = 0;
    }
    if (!*ready) {
        if (pppm_init(pppm, mesh, state[3], order, m->n) < 0) { PyErr_SetString(PyExc_MemoryError, "Could not allocate the PPPM mesh."); return -1; }
        *ready = 1;
    }
    m->g = state[3];
    return 0;
}

static int get_buffer(PyObject *object, Py_buffer *view, int writable, Py_ssize_t itemsize, const char *format, int ndim, const Py_ssize_t *shape, const char *name) {
    if (PyObject_GetBuffer(object, view, PyBUF_C_CONTIGUOUS | PyBUF_FORMAT | (writable ? PyBUF_WRITABLE : 0)) < 0) return -1;
    int ok = view->itemsize == itemsize && view->format && strcmp(view->format, format) == 0 && view->ndim == ndim;
    for (int d = 0; ok && d < ndim; ++d) ok = shape[d] < 0 || view->shape[d] == shape[d];
    if (!ok) { PyErr_Format(PyExc_ValueError, "%s has the wrong dtype or shape.", name); return -1; }
    return 0;
}

PyObject *advance_molecular(PyObject *self, PyObject *args, PyObject *kwargs) {
    (void)self;
    static char *keywords[] = {"x", "v", "u", "f", "box", "mass", "types", "charges", "epsilon", "sigma", "pppm_state", "deform",
        "steps", "timestep", "cutoff", "coulomb_cutoff", "skin", "temperature", "friction", "seed", "heat_rate", "pressure",
        "barostat_mass", "rate", "mvv2e", "qqr2e", "molecule_size", "constraints", "alpha", "lj_shift", "coulomb",
        "accuracy", "order", NULL};
    enum { X, V, U, F, BOX, MASS, TYPES, CHARGES, EPSILON, SIGMA, STATE, BUFFERS };
    PyObject *objects[BUFFERS], *deform_object, *result = NULL;
    Py_buffer view[BUFFERS + 1];
    memset(view, 0, sizeof view);
    int requested, size, lj_shift, coulomb, order, done = 0, interrupted = 0, pppm_ready = 0;
    double dt, cutoff, coulomb_cutoff, skin, T, gamma, heat, targetP, baromass, rate, mvv2e, qqr2e, alpha, accuracy;
    double constraint[3];
    unsigned long long seed;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOOOOOiddddddKddddddi(ddd)dpidi", keywords,
            &objects[X], &objects[V], &objects[U], &objects[F], &objects[BOX], &objects[MASS], &objects[TYPES],
            &objects[CHARGES], &objects[EPSILON], &objects[SIGMA], &objects[STATE], &deform_object,
            &requested, &dt, &cutoff, &coulomb_cutoff, &skin, &T, &gamma, &seed, &heat, &targetP,
            &baromass, &rate, &mvv2e, &qqr2e, &size, &constraint[0], &constraint[1], &constraint[2],
            &alpha, &lj_shift, &coulomb, &accuracy, &order)) return NULL;
    Model m = {0};
    MoleculeList nb = {0};
    PPPM pppm;
    memset(&pppm, 0, sizeof pppm);
    double *site = NULL, *offset = NULL, *site_force = NULL;
    Terms terms = {0};
    Timings timing = {0};
    double began = seconds();
    Py_ssize_t any3[2] = {-1, 3}, square[2] = {-1, -1}, vector3[1] = {3}, state8[1] = {8};
    for (int k = X; k <= F; ++k) if (get_buffer(objects[k], &view[k], 1, sizeof(double), "d", 2, any3, "Positions, velocities and forces") < 0) goto cleanup;
    int n = (int)view[X].shape[0];
    Py_ssize_t atoms[1] = {n};
    for (int k = V; k <= F; ++k) if (view[k].shape[0] != n) { PyErr_SetString(PyExc_ValueError, "State shapes must agree."); goto cleanup; }
    if (get_buffer(objects[BOX], &view[BOX], 1, sizeof(double), "d", 1, vector3, "box") < 0 ||
        get_buffer(objects[MASS], &view[MASS], 0, sizeof(double), "d", 1, atoms, "mass") < 0 ||
        get_buffer(objects[TYPES], &view[TYPES], 0, sizeof(int), "i", 1, atoms, "types") < 0 ||
        get_buffer(objects[CHARGES], &view[CHARGES], 0, sizeof(double), "d", 1, atoms, "charges") < 0 ||
        get_buffer(objects[EPSILON], &view[EPSILON], 0, sizeof(double), "d", 2, square, "epsilon") < 0 ||
        get_buffer(objects[SIGMA], &view[SIGMA], 0, sizeof(double), "d", 2, square, "sigma") < 0 ||
        get_buffer(objects[STATE], &view[STATE], 1, sizeof(double), "d", 1, state8, "pppm_state") < 0) goto cleanup;
    int types = (int)view[EPSILON].shape[0];
    if (view[EPSILON].shape[1] != types || view[SIGMA].shape[0] != types || view[SIGMA].shape[1] != types) {
        PyErr_SetString(PyExc_ValueError, "epsilon and sigma tables must be square and equal."); goto cleanup;
    }
    double *deform = NULL;
    if (deform_object != Py_None) {
        Py_ssize_t six[1] = {6};
        if (get_buffer(deform_object, &view[BUFFERS], 0, sizeof(double), "d", 1, six, "deform") < 0) goto cleanup;
        deform = view[BUFFERS].buf;
    }
    for (int a = 0; a <= BUFFERS; ++a) for (int b = a+1; b <= BUFFERS; ++b) {
        if (!view[a].obj || !view[b].obj) continue;
        uintptr_t p = (uintptr_t)view[a].buf, q = (uintptr_t)view[b].buf;
        if (p < q + view[b].len && q < p + view[a].len) { PyErr_SetString(PyExc_ValueError, "Core array buffers must not overlap."); goto cleanup; }
    }
    double *x = view[X].buf, *v = view[V].buf, *u = view[U].buf, *f = view[F].buf, *box = view[BOX].buf, *state = view[STATE].buf;
    m.n = n; m.size = size; m.types = types; m.coulomb = coulomb;
    m.mass = view[MASS].buf; m.type = view[TYPES].buf; m.charge = view[CHARGES].buf;
    m.epsilon = view[EPSILON].buf; m.sigma = view[SIGMA].buf;
    m.cutoff = cutoff; m.coulomb_cutoff = coulomb ? coulomb_cutoff : 0; m.alpha = alpha;
    m.qqr2e = qqr2e; m.mvv2e = mvv2e;
    memcpy(m.constraint, constraint, sizeof constraint);
    if ((size != 1 && size != 3) || n < 2 || n > 20000 || n % size || requested < 0 || !(dt > 0) || !(cutoff > 0) ||
        !(skin > 0) || T < 0 || gamma < 0 || !(mvv2e > 0) || !(qqr2e > 0) || !isfinite(dt+cutoff+skin+T+gamma+heat+targetP+baromass+rate+mvv2e+qqr2e+alpha) ||
        (targetP >= 0 && !(baromass > 0)) || (coulomb && (!(coulomb_cutoff > 0) || !(accuracy > 0) || order < 2 || order > 7)) ||
        (size == 3 && (!(constraint[0] > 0) || !(constraint[1] > 0) || !(constraint[2] > 0) || alpha < 0 || alpha >= 0.5)) ||
        (deform && (targetP >= 0 || !(deform[3] > 0) || !(deform[5] >= 1) || deform[4] < 0))) {
        PyErr_SetString(PyExc_ValueError, "Invalid molecular integrator parameters."); goto cleanup;
    }
    m.molecules = n/size;
    for (int i = 0; i < n; ++i) {
        if (!(m.mass[i] > 0) || !isfinite(m.mass[i]) || m.type[i] < 0 || m.type[i] >= types || !isfinite(m.charge[i])) {
            PyErr_SetString(PyExc_ValueError, "Masses must be positive, types valid and charges finite."); goto cleanup;
        }
        for (int d = 0; d < 3; ++d) if (!isfinite(x[3*i+d]) || !isfinite(v[3*i+d]) || !isfinite(u[3*i+d])) {
            PyErr_SetString(PyExc_ValueError, "State arrays must be finite."); goto cleanup;
        }
    }
    m.shift = calloc((size_t)types*types, sizeof(double));
    if (!m.shift) { PyErr_NoMemory(); goto cleanup; }
    for (int k = 0; k < types*types; ++k) {
        if (!(m.epsilon[k] >= 0) || !(m.sigma[k] >= 0) || (m.epsilon[k] > 0 && !(m.sigma[k] > 0)) || !isfinite(m.epsilon[k] + m.sigma[k])) {
            PyErr_SetString(PyExc_ValueError, "Pair epsilon must be nonnegative and sigma positive where epsilon > 0."); goto cleanup;
        }
        double s6 = pow(m.sigma[k]/cutoff, 6);
        if (lj_shift) m.shift[k] = 4*m.epsilon[k]*s6*(s6 - 1);
    }
    if (size == 3) m.extent = fmax(constraint[0], constraint[1]);
    double reach = fmax(m.cutoff, m.coulomb_cutoff) + 2*m.extent;
    for (int d = 0; d < 3; ++d) if (!(box[d] > 2*reach) || !isfinite(box[d])) {
        PyErr_SetString(PyExc_ValueError, "Each box length must exceed twice (cutoff + molecular diameter)."); goto cleanup;
    }
    nb.radius = reach + skin;
    nb.reference = malloc((size_t)n*3*sizeof(double));
    site = malloc((size_t)n*3*sizeof(double));
    offset = malloc((size_t)n*3*sizeof(double));
    site_force = malloc((size_t)n*3*sizeof(double));
    if (!nb.reference || !site || !offset || !site_force) { PyErr_NoMemory(); goto cleanup; }
    if (coulomb && prepare_pppm(&m, &pppm, &pppm_ready, state, box, accuracy, order) < 0) goto cleanup;
    if (build_list(&nb, &m, x, u, box) < 0 || forces(&m, &nb, &pppm, x, box, site, offset, site_force, f, &terms, &timing) < 0) goto cleanup;
    uint64_t rng = seed ? seed : UINT64_C(88172645463325252);
    const double dof = 3.0*n - 3 - (size == 3 ? 3.0*m.molecules : 0);
    const double centre_dof = 3.0*m.molecules - 3;
    /* Molecular MTK: barostat couples to centre-of-mass motion and the molecular
       virial. Conserved with gamma=heat=0: K + U + P*V + baromass*rate^2/2. */
    while (done < requested) {
        if (targetP >= 0) {
            double K = centre_kinetic(&m, v), volume = box[0]*box[1]*box[2];
            rate += 0.5*dt*((2*K + terms.molecular_virial - 3*targetP*volume) + 6*K/centre_dof)/baromass;
            scale_centre_velocities(&m, v, exp(-0.5*dt*(1 + 3/centre_dof)*rate));
        }
        kick(&m, v, f, 0.5*dt);
        if (targetP >= 0) {
            double factor = exp(0.5*dt*rate), factors[3] = {factor, factor, factor};
            if (!isfinite(factor) || factor < 0.9 || factor > 1.1) { PyErr_SetString(PyExc_ValueError, "Unstable barostat; increase pressure_time or reduce timestep."); goto cleanup; }
            scale_molecules(&m, x, u, box, factors);
            for (int d = 0; d < 3; ++d) box[d] *= factor;
        }
        if (gamma > 0) {
            if (constrained_drift(&m, x, u, v, box, 0.5*dt) < 0) goto cleanup;
            thermostat(&m, v, T, gamma, dt, &rng);
            if (size == 3 && constrain_velocities(&m, v, x, box) < 0) goto cleanup;
            if (constrained_drift(&m, x, u, v, box, 0.5*dt) < 0) goto cleanup;
        } else if (constrained_drift(&m, x, u, v, box, dt) < 0) goto cleanup;
        if (targetP >= 0) {
            double factor = exp(0.5*dt*rate), factors[3] = {factor, factor, factor};
            scale_molecules(&m, x, u, box, factors);
            for (int d = 0; d < 3; ++d) box[d] *= factor;
        }
        if (deform) {
            double step = deform[4] + done + 1, factors[3];
            for (int d = 0; d < 3; ++d) factors[d] = deform[d]*(1 + (deform[3] - 1)*step/deform[5])/box[d];
            scale_molecules(&m, x, u, box, factors);
            for (int d = 0; d < 3; ++d) box[d] = deform[d]*(1 + (deform[3] - 1)*step/deform[5]);
        }
        for (int d = 0; d < 3; ++d) if (!(box[d] > 2*reach)) {
            PyErr_SetString(PyExc_ValueError, "Box contracted below twice (cutoff + molecular diameter)."); goto cleanup;
        }
        if (coulomb && prepare_pppm(&m, &pppm, &pppm_ready, state, box, accuracy, order) < 0) goto cleanup;
        double listed = seconds();
        if (needs_rebuild(&nb, n, u, box, skin) && build_list(&nb, &m, x, u, box) < 0) goto cleanup;
        timing.list += seconds() - listed;
        if (forces(&m, &nb, &pppm, x, box, site, offset, site_force, f, &terms, &timing) < 0) goto cleanup;
        kick(&m, v, f, 0.5*dt);
        if (size == 3 && constrain_velocities(&m, v, x, box) < 0) goto cleanup;
        if (targetP >= 0) {
            scale_centre_velocities(&m, v, exp(-0.5*dt*(1 + 3/centre_dof)*rate));
            double K = centre_kinetic(&m, v), volume = box[0]*box[1]*box[2];
            rate += 0.5*dt*((2*K + terms.molecular_virial - 3*targetP*volume) + 6*K/centre_dof)/baromass;
        }
        if (heat != 0) {
            double K = kinetic(&m, v);
            if (K <= 0 || K + heat*dt <= 0) { PyErr_SetString(PyExc_ValueError, "Heat removal would exhaust kinetic energy."); goto cleanup; }
            double factor = sqrt((K + heat*dt)/K);
            for (int k = 0; k < 3*n; ++k) v[k] *= factor;
        }
        ++done;
        if ((done % 16 == 0 || done == requested) && PyErr_CheckSignals() < 0) {
            if (!PyErr_ExceptionMatches(PyExc_KeyboardInterrupt)) goto cleanup;
            PyErr_Clear(); interrupted = 1; break;
        }
    }
    (void)dof;
    timing.total = seconds() - began;
    result = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:K,s:i,s:i,s:d,s:{s:d,s:d,s:d,s:d}}",
        "lj", terms.lj, "coulomb", terms.coulomb, "kspace", terms.kspace,
        "virial", terms.virial, "molecular_virial", terms.molecular_virial,
        "pressure_numerator", 2*centre_kinetic(&m, v) + terms.molecular_virial,
        "rng", (unsigned long long)rng, "done", done, "interrupted", interrupted, "rate", rate,
        "timings", "list", timing.list, "pair", timing.pair, "kspace", timing.kspace, "total", timing.total);
cleanup:
    for (int k = 0; k <= BUFFERS; ++k) if (view[k].obj) PyBuffer_Release(&view[k]);
    free(m.shift); free(nb.pairs); free(nb.reference); free(site); free(offset); free(site_force);
    if (pppm_ready) pppm_free(&pppm);
    return result;
}

PyObject *project_molecular_velocities(PyObject *self, PyObject *args) {
    (void)self;
    PyObject *objects[4];
    double constraint[3];
    if (!PyArg_ParseTuple(args, "OOOO(ddd)", &objects[0], &objects[1], &objects[2], &objects[3],
                          &constraint[0], &constraint[1], &constraint[2])) return NULL;
    Py_buffer view[4];
    memset(view, 0, sizeof view);
    PyObject *result = NULL;
    Py_ssize_t any3[2] = {-1, 3}, vector3[1] = {3}, any[1] = {-1};
    if (get_buffer(objects[0], &view[0], 0, sizeof(double), "d", 2, any3, "positions") < 0 ||
        get_buffer(objects[1], &view[1], 1, sizeof(double), "d", 2, any3, "velocities") < 0 ||
        get_buffer(objects[2], &view[2], 0, sizeof(double), "d", 1, vector3, "box") < 0 ||
        get_buffer(objects[3], &view[3], 0, sizeof(double), "d", 1, any, "mass") < 0) goto cleanup;
    int n = (int)view[0].shape[0];
    if (view[1].shape[0] != n || view[3].shape[0] != n || n % 3) { PyErr_SetString(PyExc_ValueError, "Water arrays must hold whole molecules."); goto cleanup; }
    Model m = {0};
    m.n = n; m.size = 3; m.molecules = n/3; m.mass = view[3].buf;
    memcpy(m.constraint, constraint, sizeof constraint);
    if (constrain_velocities(&m, view[1].buf, view[0].buf, view[2].buf) < 0) goto cleanup;
    result = Py_NewRef(Py_None);
cleanup:
    for (int k = 0; k < 4; ++k) if (view[k].obj) PyBuffer_Release(&view[k]);
    return result;
}

PyObject *pppm_parameters(PyObject *self, PyObject *args) {
    (void)self;
    double box[3], qsqsum, qqr2e, cutoff, accuracy, g, estimate;
    int atoms, order, mesh[3];
    if (!PyArg_ParseTuple(args, "(ddd)iddddi", &box[0], &box[1], &box[2], &atoms, &qsqsum, &qqr2e, &cutoff, &accuracy, &order)) return NULL;
    if (pppm_select(box, atoms, qsqsum, qqr2e, cutoff, accuracy, order, mesh, &g, &estimate) < 0) {
        PyErr_SetString(PyExc_ValueError, "Could not choose PPPM parameters."); return NULL;
    }
    return Py_BuildValue("((iii)dd)", mesh[0], mesh[1], mesh[2], g, estimate);
}

PyObject *fft_transform(PyObject *self, PyObject *args) {
    (void)self;
    PyObject *object;
    int sign;
    if (!PyArg_ParseTuple(args, "Oi", &object, &sign)) return NULL;
    Py_buffer view = {0};
    Py_ssize_t shape[3] = {-1, -1, -1};
    if (get_buffer(object, &view, 1, 2*sizeof(double), "Zd", 3, shape, "data") < 0) {
        if (view.obj) PyBuffer_Release(&view);
        return NULL;
    }
    FFT3D fft;
    int n[3] = {(int)view.shape[0], (int)view.shape[1], (int)view.shape[2]};
    if (fft3d_init(&fft, n) < 0) { PyBuffer_Release(&view); PyErr_SetString(PyExc_ValueError, "FFT sizes must have factors 2, 3 and 5 only."); return NULL; }
    fft3d(&fft, view.buf, sign < 0 ? -1 : 1);
    fft3d_free(&fft);
    PyBuffer_Release(&view);
    Py_RETURN_NONE;
}
