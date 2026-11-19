#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "molecular.h"

typedef struct { int i, j; } Pair;
/* Pairs within the list radius, ordered by i and then j whichever way they were
   found, so forces are summed in one fixed order. reference and box are the
   unwrapped positions and box at the build. The rest is work space for the
   cell build, kept between builds. */
typedef struct {
    Pair *pairs; size_t size, capacity; double *reference; double box[3];
    Pair *found; size_t found_capacity, *offset;
    double *sorted;
    int *cell, *members, *first, first_capacity;
} Neighbors;
/* What one run keeps between calls into the solver: the neighbour list, the
   forces (f, which the solver works in and copies out when it returns), and the
   positions, box and parameters exactly as the previous call left them. A call
   that finds them bit-for-bit unchanged continues with that list and those
   forces; any difference (positions, box, cutoff, atom count, pair or bond
   parameters) makes it start afresh. */
#define WORKSPACE_SCALARS 12
/* A direction without periodic images (a wall on its two faces, or the direction
   a flat system lacks) is given this period: so long that the nearest image of a
   neighbour is always the neighbour itself, and the pair loops need no branch. */
#define NO_IMAGE 1e30
/* One-body potentials that depend on one coordinate, one row each: kind, axis and
   three parameters. Every term reports the force the atoms exert on it. The kinds
   so far are the walls: kind 0 at coordinate 0 and kind 1 at box[axis], with
   U = stiffness/3 * depth^3 for an atom that has passed the wall by depth. The
   force, stiffness * depth^2, starts from zero without a kink, so a time step
   that straddles the wall costs no more accuracy than one between two atoms
   (a spring wall of the same hardness lost about ten times more energy accuracy).
   Parameters: stiffness, speed, target. A wall of kind 1 with speed > 0 is a
   piston: it moves towards target at that speed and stops there.
   For each term the caller keeps three numbers: the force on it now (counted
   outwards for a wall), the impulse it has received and the work it has done on
   the atoms by moving; the solver adds to the last two. */
#define EXTERNAL_MAX 8
#define EXTERNAL_WIDTH 5
#define EXTERNAL_SUMS 3
typedef struct {
    Neighbors nb;
    int n, current, externals;
    double *x, *u, *f, *tables, box[3], scalars[WORKSPACE_SCALARS], potential, virial;
    double external[EXTERNAL_MAX*EXTERNAL_WIDTH], external_force[EXTERNAL_MAX];
    size_t table_size[2];
} Workspace;
static const char WORKSPACE[] = "fys2160_md._core.workspace";

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
static int reserve_pairs(Pair **pairs, size_t *capacity, size_t needed) {
    if (needed <= *capacity) return 0;
    size_t cap = *capacity ? *capacity : 1024;
    while (cap < needed) cap *= 2;
    Pair *p = realloc(*pairs, cap * sizeof(Pair));
    if (!p) { PyErr_NoMemory(); return -1; }
    *pairs = p; *capacity = cap;
    return 0;
}
/* Cells at least radius/reach wide, so that every pair within radius lies at
   most reach cells apart. Returns the reach (1 or 2) with the cheapest build, or
   0 when checking all pairs is cheaper or the box holds too few cells. The
   costs are rough times per atom in nanoseconds, measured for this code: one
   distance per candidate, and a fixed part for each run of cells that is scanned.
   only = 1 or 2 asks for that reach whatever it costs. A flat system
   (dimension 2) has one layer of cells. */
static int choose_cells(int n, const double *box, int dimension, double radius, int grid[3], int only) {
    double best = only ? HUGE_VAL : 0.55 * n;
    int reach = 0;
    for (int k = 1; k <= 2; ++k) {
        if (only && only != k) continue;
        int g[3] = {1, 1, 1}, width = 2*k + 1, fits = 1, scanned = dimension == 3 ? width*width : width;
        int limit = (int)(dimension == 3 ? cbrt(4.0*k*k*k*n) : sqrt(4.0*k*k*n)) + 1;
        for (int d = 0; d < dimension; ++d) {
            /* The margin keeps rounding of cell indices from hiding a pair. */
            double count = floor(box[d] * k / (radius * (1 + 1e-9)));
            g[d] = count < limit ? (int)count : limit;
            if (g[d] < width) fits = 0;
        }
        if (!fits) continue;
        double cost = 8.0*scanned + 0.9*scanned*width*n/((double)g[0]*g[1]*g[2]);
        if (cost < best) { best = cost; reach = k; memcpy(grid, g, sizeof g); }
    }
    return reach;
}
/* box gives the extent the cells divide, period the image distances: the same
   numbers, except NO_IMAGE in a direction that is not periodic. There an atom
   beyond a wall is put in the outermost cell, which is at least as near to its
   neighbours as the atom is. */
static int build_with_cells(Neighbors *nb, int n, const double *x, const double *box, const double *period,
                            int dimension, const int *periodic, double rmax2, const int grid[3], int reach) {
    int cells = grid[0]*grid[1]*grid[2];
    if (cells + 1 > nb->first_capacity) {
        int *first = realloc(nb->first, ((size_t)cells + 1) * sizeof(int));
        if (!first) { PyErr_NoMemory(); return -1; }
        nb->first = first; nb->first_capacity = cells + 1;
    }
    /* Atoms sorted by cell, cells counted along x first: members[first[c] .. first[c+1])
       are the atoms of cell c and sorted holds their positions, so a run of cells
       along x is one stretch of memory. */
    int *first = nb->first, *members = nb->members, *cell = nb->cell;
    double *sorted = nb->sorted;
    memset(first, 0, ((size_t)cells + 1) * sizeof(int));
    for (int i = 0; i < n; ++i) for (int d = 0; d < 3; ++d) {
        if (d >= dimension) { cell[3*i+d] = 0; continue; }
        double t = x[3*i+d] / box[d];
        t = periodic[d] ? t - floor(t) : t < 0 ? 0 : t > 1 ? 1 : t;
        int c = (int)(t * grid[d]);
        cell[3*i+d] = c < grid[d] ? c : grid[d] - 1;
    }
    for (int i = 0; i < n; ++i) ++first[(cell[3*i+2]*grid[1] + cell[3*i+1])*grid[0] + cell[3*i] + 1];
    for (int c = 0; c < cells; ++c) first[c+1] += first[c];
    for (int i = 0; i < n; ++i) {
        int m = first[(cell[3*i+2]*grid[1] + cell[3*i+1])*grid[0] + cell[3*i]]++;
        members[m] = i;
        memcpy(sorted + 3*m, x + 3*i, 3*sizeof(double));
    }
    for (int c = cells; c > 0; --c) first[c] = first[c-1];
    first[0] = 0;
    /* The partners i < j of each j, taking j in rising order; ordering the pairs
       by i afterwards leaves the partners of each i in rising order too. The
       product with 1/box picks another image than the solver's division only for
       a separation of half a box, which is beyond the radius when cells are used. */
    const double inverse[3] = {1/period[0], 1/period[1], 1/period[2]};
    const int reach2 = dimension == 3 ? reach : 0;
    size_t found = 0;
    for (int j = 0; j < n; ++j) {
        const double xj[3] = {x[3*j], x[3*j+1], x[3*j+2]};
        /* The cells within reach along x, as one run or as two when it wraps. */
        int low = cell[3*j] - reach, high = cell[3*j] + reach, runs = 1, run[2][2] = {{low, high}, {0, 0}};
        if (!periodic[0]) { run[0][0] = low < 0 ? 0 : low; run[0][1] = high < grid[0] ? high : grid[0] - 1; }
        else if (low < 0) { runs = 2; run[0][0] = 0; run[1][0] = low + grid[0]; run[1][1] = grid[0] - 1; }
        else if (high >= grid[0]) { runs = 2; run[0][1] = grid[0] - 1; run[1][0] = 0; run[1][1] = high - grid[0]; }
        for (int a2 = -reach2; a2 <= reach2; ++a2) for (int a1 = -reach; a1 <= reach; ++a1) {
            int c2 = cell[3*j+2] + a2, c1 = cell[3*j+1] + a1;
            /* Beyond a wall there are no cells. */
            if ((!periodic[2] && (c2 < 0 || c2 >= grid[2])) || (!periodic[1] && (c1 < 0 || c1 >= grid[1]))) continue;
            c2 += c2 < 0 ? grid[2] : c2 >= grid[2] ? -grid[2] : 0;
            c1 += c1 < 0 ? grid[1] : c1 >= grid[1] ? -grid[1] : 0;
            int row = (c2*grid[1] + c1)*grid[0];
            for (int r = 0; r < runs; ++r) {
                int begin = first[row + run[r][0]], end = first[row + run[r][1] + 1];
                if (found + (size_t)(end - begin) > nb->found_capacity &&
                    reserve_pairs(&nb->found, &nb->found_capacity, found + (size_t)(end - begin)) < 0) return -1;
                Pair *out = nb->found;
                for (int m = begin; m < end; ++m) {
                    double r2 = 0;
                    for (int d = 0; d < 3; ++d) {
                        double dx = sorted[3*m+d] - xj[d];
                        dx -= period[d] * nearbyint(dx * inverse[d]);
                        r2 += dx * dx;
                    }
                    out[found] = (Pair){members[m], j};
                    found += (r2 < rmax2) & (members[m] < j);
                }
            }
        }
    }
    if (reserve_pairs(&nb->pairs, &nb->capacity, found) < 0) return -1;
    size_t *offset = nb->offset;
    memset(offset, 0, ((size_t)n + 1) * sizeof(size_t));
    for (size_t p = 0; p < found; ++p) ++offset[nb->found[p].i + 1];
    for (int i = 0; i < n; ++i) offset[i+1] += offset[i];
    for (size_t p = 0; p < found; ++p) nb->pairs[offset[nb->found[p].i]++] = nb->found[p];
    nb->size = found;
    return 0;
}
/* method: 0 checks all pairs, 1 or 2 uses cells with that reach, -1 takes the
   cheapest. All give the same list in the same order. */
static int build_neighbors(Neighbors *nb, int n, const double *x, const double *u, const double *box,
                           const double *period, int dimension, const int *periodic, double radius, int method) {
    double rmax2 = radius * radius;
    int grid[3], reach = method ? choose_cells(n, box, dimension, radius, grid, method > 0 ? method : 0) : 0;
    if (method > 0 && !reach) {
        PyErr_SetString(PyExc_ValueError, "The box is too small for cells of this reach.");
        return -1;
    }
    if (reach) {
        if (build_with_cells(nb, n, x, box, period, dimension, periodic, rmax2, grid, reach) < 0) return -1;
    } else {
        /* In a box this small a pair half a box apart may be listed under either
           image; the two are equally far, and the forces choose their own. */
        const double inverse[3] = {1/period[0], 1/period[1], 1/period[2]};
        size_t size = 0;
        for (int i = 0; i < n; ++i) {
            const double xi[3] = {x[3*i], x[3*i+1], x[3*i+2]};
            if (reserve_pairs(&nb->pairs, &nb->capacity, size + (size_t)(n - i)) < 0) return -1;
            Pair *out = nb->pairs;
            for (int j = i + 1; j < n; ++j) {
                double r2 = 0;
                for (int d = 0; d < 3; ++d) {
                    double dx = xi[d] - x[3*j+d];
                    dx -= period[d] * nearbyint(dx * inverse[d]);
                    r2 += dx * dx;
                }
                out[size] = (Pair){i, j};
                size += r2 < rmax2;
            }
        }
        nb->size = size;
    }
    memcpy(nb->reference, u, (size_t)n * 3 * sizeof(double));
    memcpy(nb->box, period, 3*sizeof(double));
    return 0;
}
static void release_workspace(Workspace *w) {
    free(w->nb.pairs); free(w->nb.reference); free(w->nb.found); free(w->nb.offset);
    free(w->nb.cell); free(w->nb.members); free(w->nb.first); free(w->nb.sorted);
    free(w->x); free(w->u); free(w->f); free(w->tables);
    memset(w, 0, sizeof *w);
}
/* Size the arrays for n atoms; a different atom count discards what was kept. */
static int prepare_workspace(Workspace *w, int n) {
    if (w->n == n) return 0;
    release_workspace(w);
    size_t count = (size_t)n * 3 * sizeof(double);
    w->nb.reference = malloc(count); w->nb.sorted = malloc(count); w->x = malloc(count); w->u = malloc(count); w->f = malloc(count);
    w->tables = malloc((size_t)n * 4 * sizeof(double));
    w->nb.offset = malloc(((size_t)n + 1) * sizeof(size_t));
    w->nb.cell = malloc((size_t)n * 3 * sizeof(int)); w->nb.members = malloc((size_t)n * sizeof(int));
    if (!w->nb.reference || !w->nb.sorted || !w->x || !w->u || !w->f || !w->tables || !w->nb.offset || !w->nb.cell || !w->nb.members) {
        release_workspace(w); PyErr_NoMemory(); return -1;
    }
    w->n = n;
    return 0;
}
static int same(const double *kept, const double *now, size_t count) {
    return !count || memcmp(kept, now, count * sizeof(double)) == 0;
}
/* table a and table b are constant inputs of the forces (at most 2n doubles each). */
static int unchanged(const Workspace *w, int n, const double *x, const double *u, const double *box,
                     const double *scalars, const double *a, size_t na, const double *b, size_t nb) {
    return w->current && w->n == n && same(w->scalars, scalars, WORKSPACE_SCALARS) && same(w->box, box, 3) &&
           same(w->x, x, (size_t)3*n) && (!u || same(w->u, u, (size_t)3*n)) &&
           w->table_size[0] == na && w->table_size[1] == nb && same(w->tables, a, na) && same(w->tables + 2*(size_t)n, b, nb);
}
/* tables is 0 when the call began with them unchanged. */
static void remember(Workspace *w, int n, const double *x, const double *u, const double *box,
                     const double *scalars, int tables, const double *a, size_t na, const double *b, size_t nb,
                     double potential, double virial) {
    size_t count = (size_t)n * 3 * sizeof(double);
    memcpy(w->x, x, count);
    if (u) memcpy(w->u, u, count);
    memcpy(w->box, box, sizeof w->box); memcpy(w->scalars, scalars, sizeof w->scalars);
    if (tables) {
        if (na) memcpy(w->tables, a, na * sizeof(double));
        if (nb) memcpy(w->tables + 2*(size_t)n, b, nb * sizeof(double));
        w->table_size[0] = na; w->table_size[1] = nb;
    }
    w->potential = potential; w->virial = virial; w->current = 1;
}
static void free_workspace(PyObject *capsule) {
    Workspace *w = PyCapsule_GetPointer(capsule, WORKSPACE);
    if (w) { release_workspace(w); free(w); }
}
static PyObject *workspace(PyObject *self, PyObject *args) {
    (void)self; (void)args;
    Workspace *w = calloc(1, sizeof *w);
    if (!w) return PyErr_NoMemory();
    PyObject *capsule = PyCapsule_New(w, WORKSPACE, free_workspace);
    if (!capsule) free(w);
    return capsule;
}
/* An optional trailing argument: None or absent gives a workspace for this call only. */
static Workspace *borrowed_workspace(PyObject *object, Workspace *local) {
    if (!object || object == Py_None) return local;
    return PyCapsule_GetPointer(object, WORKSPACE);
}
static int needs_rebuild(const Neighbors *nb, int n, const double *u,
                         const double *box, double cutoff, double skin) {
    double smallest_scale=box[0]/nb->box[0];
    for(int d=1;d<3;++d) smallest_scale=fmin(smallest_scale,box[d]/nb->box[d]);
    double allowance=0.5*(smallest_scale*(cutoff+skin)-cutoff);
    if(allowance<=0)return 1;
    for(int i=0;i<n;++i) {
        double displacement2=0;
        for(int d=0;d<3;++d) {
            double dx=u[3*i+d]-nb->reference[3*i+d]*box[d]/nb->box[d];
            displacement2+=dx*dx;
        }
        if(displacement2>allowance*allowance)return 1;
    }
    return 0;
}
/* Independent nonbonded LJ selection and bonded-pair exclusion.
   Flexible bonds use a polynomial (including the harmonic special case). */
static int forces(const Neighbors *nb, int n, const double *x, const double *box,
                  double cutoff, int model, int pair, double epsilon, double sigma, const double *parameters, int exclude, double length, double k2, double k3, double k4, const double *bond_parameters, double *f, double *potential, double *virial) {
    memset(f, 0, (size_t)n * 3 * sizeof(double));
    *potential = 0; *virial = 0;

    double rc2 = cutoff*cutoff;
    for (size_t p = 0; p < nb->size; ++p) {
        int i = nb->pairs[p].i, j = nb->pairs[p].j;
        if(model && exclude && i/2==j/2)continue; /* bonded intramolecular pair excluded */
        double dr[3], r2 = 0;
        for (int d = 0; d < 3; ++d) {
            dr[d] = minimum_image(x[3*i+d] - x[3*j+d], box[d]);
            r2 += dr[d] * dr[d];
        }
        if (r2 >= rc2) continue;
        if (r2 < 1.e-16 || !isfinite(r2)) {
            PyErr_SetString(PyExc_ValueError, "Overlapping atoms or unstable trajectory; check positions and timestep.");
            return -1;
        }
        if(parameters){epsilon=sqrt(parameters[2*i])*sqrt(parameters[2*j]);sigma=.5*parameters[2*i+1]+.5*parameters[2*j+1];}
        double invc2=sigma*sigma/(cutoff*cutoff),invc6=invc2*invc2*invc2;
        double shift=epsilon*(pair?2*invc6*sqrt(invc6)-3*invc6:4*invc6*(invc6-1));
        double inv2 = 1/r2, s2 = sigma*sigma*inv2, inv6 = s2*s2*s2;
        double coefficient;
        if(pair) {
            double inv9=inv6*sqrt(inv6);
            coefficient=epsilon*18*(inv9-inv6)*inv2;
            *potential+=epsilon*(2*inv9-3*inv6)-shift;
        } else {
            coefficient=epsilon*24*inv6*(2*inv6-1)*inv2;
            *potential+=epsilon*4*inv6*(inv6-1)-shift;
        }
        *virial += coefficient * r2;
        for (int d = 0; d < 3; ++d) {
            double force = coefficient * dr[d];
            f[3*i+d] += force; f[3*j+d] -= force;
        }
    }
    if(model==1)for(int i=0;i<n;i+=2) {
        if(bond_parameters){const double *b=bond_parameters+4*(i/2);length=b[0];k2=b[1];k3=b[2];k4=b[3];}
        double dr[3],r2=0;
        for(int d=0;d<3;++d){dr[d]=minimum_image(x[3*i+d]-x[3*(i+1)+d],box[d]);r2+=dr[d]*dr[d];}
        double r=sqrt(r2),q=r-length;
        if(r<1.e-8){PyErr_SetString(PyExc_ValueError,"Collapsed molecular bond.");return -1;}
        double coefficient=-(2*k2*q+3*k3*q*q+4*k4*q*q*q)/r;
        *potential+=k2*q*q+k3*q*q*q+k4*q*q*q*q;
        *virial+=coefficient*r2;
        for(int d=0;d<3;++d){f[3*i+d]+=coefficient*dr[d];f[3*(i+1)+d]-=coefficient*dr[d];}
    }
    if (!isfinite(*potential) || !isfinite(*virial)) {
        PyErr_SetString(PyExc_FloatingPointError, "Non-finite energy; reduce timestep or relax the initial positions.");
        return -1;
    }
    return 0;
}
/* Positions are wrapped into the box only along periodic directions. */
static int drift(int n, double *x, double *u, const double *v,
                 const double *box, const int *periodic, double duration) {
    for (int i = 0; i < n; ++i) for (int d = 0; d < 3; ++d) {
        int k = 3*i+d;
        double dx = duration * v[k];
        u[k] += dx; x[k] += dx;
        if (!isfinite(u[k]) || !isfinite(x[k])) {
            PyErr_SetString(PyExc_FloatingPointError, "Non-finite position; simulation is unstable.");
            return -1;
        }
        if (periodic[d]) x[k] -= box[d] * floor(x[k]/box[d]);
    }
    return 0;
}
static void kick(int n, double *v, const double *f, const double *mass, double dt) {
    for (int i = 0; i < n; ++i) for (int d = 0; d < 3; ++d)
        v[3*i+d] += dt * f[3*i+d] / mass[i];
}
/* The thermostat samples the canonical distribution, with total momentum zero
   along the directions that conserve it (the periodic ones; walls take up momentum). */
static void thermostat(int n, double *v, const double *mass, double T,
                       double gamma, double dt, uint64_t *rng, int dimension, const int *conserved) {
    double c = exp(-gamma * dt), spread = -expm1(-2*gamma*dt) * T, momentum[3] = {0,0,0}, totalmass = 0;
    for (int i = 0; i < n; ++i) {
        double sigma = sqrt(spread / mass[i]);
        totalmass += mass[i];
        for (int d = 0; d < dimension; ++d) {
            v[3*i+d] = c*v[3*i+d] + sigma*gaussian(rng);
            momentum[d] += mass[i]*v[3*i+d];
        }
    }
    for (int i = 0; i < n; ++i) for (int d = 0; d < 3; ++d)
        if (conserved[d]) v[3*i+d] -= momentum[d]/totalmass;
}
/* Adds the forces of the external terms to f and their energy to *potential.
   force[t] becomes the force of the atoms on term t: for a wall, the sum over the
   atoms that have passed it, counted positive outwards. */
static void external_forces(int terms, const double *external, int n, const double *x,
                            const double *box, double *f, double *potential, double *force) {
    for (int t = 0; t < terms; ++t) {
        const double *e = external + EXTERNAL_WIDTH*t;
        int upper = e[0] != 0, axis = (int)e[1];
        double stiffness = e[2], at = upper ? box[axis] : 0, side = upper ? 1 : -1, total = 0, energy = 0;
        for (int i = 0; i < n; ++i) {
            double depth = side * (x[3*i+axis] - at);
            if (depth > 0) {
                f[3*i+axis] -= side * stiffness * depth * depth;
                total += stiffness * depth * depth; energy += stiffness * depth * depth * depth / 3;
            }
        }
        force[t] = total; *potential += energy;
    }
}
/* One time step of every piston: box[axis] moves towards the target and stops on
   it. moved[t] is the displacement. The last step may be a thousandth longer, so
   that rounding in the sum of the steps does not leave one tiny step over. */
static void move_pistons(int terms, const double *external, double *box, double dt, double *moved) {
    for (int t = 0; t < terms; ++t) {
        const double *e = external + EXTERNAL_WIDTH*t;
        int axis = (int)e[1];
        moved[t] = 0;
        if (e[0] == 0 || !(e[3] > 0) || box[axis] == e[4]) continue;
        double left = e[4] - box[axis], step = e[3] * dt;
        if (fabs(left) <= 1.001 * step) { moved[t] = left; box[axis] = e[4]; }
        else { moved[t] = left > 0 ? step : -step; box[axis] += moved[t]; }
    }
}

static double kinetic(int n,const double *v,const double *mass) {
    double K=0;for(int i=0;i<n;++i)for(int d=0;d<3;++d)K+=0.5*mass[i]*v[3*i+d]*v[3*i+d];return K;
}
/* Velocity constraint, returning the trace of its impulse virial. */
static double constrain_velocities(int n,double *v,const double *x,const double *box,const double *mass) {
    double impulse_virial=0;
    for(int i=0;i<n;i+=2) {
        double r[3],r2=0,dot=0;
        for(int d=0;d<3;++d){r[d]=minimum_image(x[3*i+d]-x[3*(i+1)+d],box[d]);r2+=r[d]*r[d];dot+=r[d]*(v[3*i+d]-v[3*(i+1)+d]);}
        double inverse_mass=1/mass[i]+1/mass[i+1],lambda=-dot/(r2*inverse_mass);
        for(int d=0;d<3;++d){v[3*i+d]+=lambda*r[d]/mass[i];v[3*(i+1)+d]-=lambda*r[d]/mass[i+1];}
        impulse_virial+=lambda*r2;
    }
    return impulse_virial;
}
/* RATTLE position constraint: correction along the bond at the start of drift. */
static int rigid_drift(int n,double *x,double *u,double *v,const double *box,
                       const double *mass,double dt,double length,const double *bond_parameters,double *impulse_virial) {
    for(int i=0;i<n;i+=2) {
        if(bond_parameters)length=bond_parameters[4*(i/2)];
        double old[3],trial[3],old2=0,trial2=0,dot=0;
        for(int d=0;d<3;++d){
            old[d]=minimum_image(x[3*i+d]-x[3*(i+1)+d],box[d]);
            trial[d]=old[d]+dt*(v[3*i+d]-v[3*(i+1)+d]);
            old2+=old[d]*old[d];trial2+=trial[d]*trial[d];dot+=trial[d]*old[d];
        }
        double discriminant=dot*dot-old2*(trial2-length*length);
        if(discriminant<=0){PyErr_SetString(PyExc_ValueError,"Rigid-bond constraint failed; reduce timestep.");return -1;}
        /* Stable form of the root near zero. */
        double lambda=-(trial2-length*length)/(dot+sqrt(discriminant));
        double inverse_mass=1/mass[i]+1/mass[i+1];
        *impulse_virial+=lambda*old2/(inverse_mass*dt);
        for(int d=0;d<3;++d){
            v[3*i+d]+=lambda*old[d]/(mass[i]*inverse_mass*dt);
            v[3*(i+1)+d]-=lambda*old[d]/(mass[i+1]*inverse_mass*dt);
        }
    }
    static const int periodic[3]={1,1,1};
    return drift(n,x,u,v,box,periodic,dt);
}
/* Molecular pressure uses COM kinetic energy and intermolecular COM virial.
   It avoids noisy finite-difference constraint-force estimates for rigid bonds. */
static double molecular_pressure_numerator(const Neighbors *nb,int n,const double *x,
                 const double *v,const double *mass,const double *box,double cutoff,int pair,double epsilon,double sigma,const double *parameters) {
    double value=0;
    for(int i=0;i<n;i+=2)for(int d=0;d<3;++d){
        double momentum=mass[i]*v[3*i+d]+mass[i+1]*v[3*(i+1)+d];
        value+=momentum*momentum/(mass[i]+mass[i+1]);
    }
    for(size_t p=0;p<nb->size;++p){
        int i=nb->pairs[p].i,j=nb->pairs[p].j;
        if(i/2==j/2)continue;
        double dr[3],cm[3],r2=0;
        int a=2*(i/2),b=2*(j/2);
        for(int d=0;d<3;++d){
            dr[d]=minimum_image(x[3*i+d]-x[3*j+d],box[d]);r2+=dr[d]*dr[d];
            /* COM positions relative to interacting atom, with its same pair image. */
            int oi=(i==a?a+1:a),oj=(j==b?b+1:b);
            double offset_i=mass[oi]/(mass[i]+mass[oi])*minimum_image(x[3*oi+d]-x[3*i+d],box[d]);
            double offset_j=mass[oj]/(mass[j]+mass[oj])*minimum_image(x[3*oj+d]-x[3*j+d],box[d]);
            cm[d]=dr[d]+offset_i-offset_j;
        }
        if(r2>=cutoff*cutoff)continue;
        if(parameters){epsilon=sqrt(parameters[2*i])*sqrt(parameters[2*j]);sigma=.5*parameters[2*i+1]+.5*parameters[2*j+1];}
        double inv2=1/r2,s2=sigma*sigma*inv2,inv6=s2*s2*s2;
        double coef=epsilon*(pair ? 18*(inv6*sqrt(inv6)-inv6)*inv2 : 24*inv6*(2*inv6-1)*inv2);
        for(int d=0;d<3;++d)value+=cm[d]*coef*dr[d];
    }
    return value;
}

static PyObject *advance(PyObject *self, PyObject *args) {
    (void)self;
    PyObject *objects[6], *result = NULL, *parameter_object=NULL, *bond_object=NULL, *workspace_object=NULL;
    PyObject *external_object=NULL, *sums_object=NULL;
    Py_buffer parameter_buffer={0}, bond_buffer={0}, external_buffer={0}, sums_buffer={0};
    Workspace local={0}, *ws;
    const double *parameters=NULL, *bond_parameters=NULL;
    Py_buffer buffers[6] = {{0}};
    int requested, done = 0, interrupted = 0;
    double dt, cutoff, skin, T, gamma, heat, targetP, baromass, rate, potential=0, virial=0;
    int model, pair=-1, exclude=1, measure=1, dimension=3, terms=0;
    const double *external=NULL;
    double *sums=NULL, wall[EXTERNAL_MAX]={0}, before[EXTERNAL_MAX]={0}, moved[EXTERNAL_MAX]={0};
    double length=.7, k2=35., k3=-40., k4=59., epsilon=1., sigma=1.;
    unsigned long long seed;
    if (!PyArg_ParseTuple(args, "OOOOOOidddddKidddd|iiddddddOOOiiOO", &objects[0], &objects[1],
        &objects[2], &objects[3], &objects[4], &objects[5], &requested,
        &dt, &cutoff, &skin, &T, &gamma, &seed, &model, &heat, &targetP, &baromass, &rate, &pair, &exclude, &length, &k2, &k3, &k4, &epsilon, &sigma, &parameter_object, &bond_object, &workspace_object, &measure, &dimension, &external_object, &sums_object)) return NULL;
    if(!(ws=borrowed_workspace(workspace_object,&local)))return NULL;
    if(pair==-1)pair=(model!=0);
    for (int k=0; k<6; ++k) {
        int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT | (k<5 ? PyBUF_WRITABLE : 0);
        if (PyObject_GetBuffer(objects[k], &buffers[k], flags)<0) goto cleanup;
        if (buffers[k].itemsize != sizeof(double) || !buffers[k].format ||
            strcmp(buffers[k].format,"d") != 0) {
            PyErr_SetString(PyExc_TypeError, "Core arrays must be C-contiguous native float64 buffers.");
            goto cleanup;
        }
    }
    if (buffers[0].ndim != 2 || buffers[0].shape[1] != 3 ||
        buffers[0].shape[0] < 2 || buffers[0].shape[0] > 20000) {
        PyErr_SetString(PyExc_ValueError, "Positions must have shape (N,3), with 2 <= N <= 20000.");
        goto cleanup;
    }
    int n = (int)buffers[0].shape[0];
    if(parameter_object){
        if(PyObject_GetBuffer(parameter_object,&parameter_buffer,PyBUF_C_CONTIGUOUS|PyBUF_FORMAT)<0)goto cleanup;
        if(parameter_buffer.ndim!=2 || parameter_buffer.shape[0]!=n || parameter_buffer.shape[1]!=2 ||
           parameter_buffer.itemsize!=sizeof(double) || !parameter_buffer.format || strcmp(parameter_buffer.format,"d")!=0){
            PyErr_SetString(PyExc_ValueError,"Pair parameters must be a contiguous float64 (N,2) array.");goto cleanup;
        }
        parameters=parameter_buffer.buf;
        for(int i=0;i<5;++i){
            uintptr_t a=(uintptr_t)parameter_buffer.buf,b=(uintptr_t)buffers[i].buf;
            if(a<b+buffers[i].len && b<a+parameter_buffer.len){
                PyErr_SetString(PyExc_ValueError,"Pair parameters must not overlap mutable state.");goto cleanup;
            }
        }
    }
    if(bond_object){
        if(PyObject_GetBuffer(bond_object,&bond_buffer,PyBUF_C_CONTIGUOUS|PyBUF_FORMAT)<0)goto cleanup;
        int count=model?n/2:0;
        if(bond_buffer.ndim!=2 || bond_buffer.shape[0]!=count || bond_buffer.shape[1]!=4 ||
           bond_buffer.itemsize!=sizeof(double) || !bond_buffer.format || strcmp(bond_buffer.format,"d")!=0){
            PyErr_SetString(PyExc_ValueError,"Bond parameters must be a contiguous float64 (number of bonds,4) array.");goto cleanup;
        }
        bond_parameters=bond_buffer.buf;
        for(int i=0;i<5;++i){
            uintptr_t a=(uintptr_t)bond_buffer.buf,b=(uintptr_t)buffers[i].buf;
            if(a<b+buffers[i].len && b<a+bond_buffer.len){
                PyErr_SetString(PyExc_ValueError,"Bond parameters must not overlap mutable state.");goto cleanup;
            }
        }
    }
    for (int k=1;k<4;++k) if (buffers[k].ndim!=2 || buffers[k].shape[0]!=n || buffers[k].shape[1]!=3) {
        PyErr_SetString(PyExc_ValueError, "Position, velocity, unwrapped and force shapes must agree.");
        goto cleanup;
    }
    if (buffers[4].ndim!=1 || buffers[4].shape[0]!=3 || buffers[5].ndim!=1 || buffers[5].shape[0]!=n) {
        PyErr_SetString(PyExc_ValueError, "Box and mass shapes must be (3,) and (N,)."); goto cleanup;
    }
    for (int a=0;a<6;++a) for(int b=a+1;b<6;++b) {
        uintptr_t p=(uintptr_t)buffers[a].buf,q=(uintptr_t)buffers[b].buf;
        if (p<q+buffers[b].len && q<p+buffers[a].len) {
            PyErr_SetString(PyExc_ValueError, "Core array buffers must not overlap."); goto cleanup;
        }
    }
    if (!(epsilon>0) || !(sigma>0) || !isfinite(epsilon+sigma) || pair<0 || pair>1 || (exclude!=0 && exclude!=1) || !(length>0) || k2<0 || k4<0 || !isfinite(length+k2+k3+k4) || requested<0 || !(dt>0) || !(cutoff>0) || !(skin>0) || T<0 || gamma<0 ||
        !isfinite(dt+cutoff+skin+T+gamma+heat+targetP+baromass+rate) || model<0 || model>2 ||
        (model && n%2) || (targetP>=0 && (model!=0 || baromass<=0))) {
        PyErr_SetString(PyExc_ValueError, "Invalid integrator parameters."); goto cleanup;
    }
    double *x=buffers[0].buf,*v=buffers[1].buf,*u=buffers[2].buf;
    double *f=NULL,*box=buffers[4].buf,*mass=buffers[5].buf;
    if(external_object && external_object!=Py_None) {
        if(!sums_object || PyObject_GetBuffer(external_object,&external_buffer,PyBUF_C_CONTIGUOUS|PyBUF_FORMAT)<0 ||
           PyObject_GetBuffer(sums_object,&sums_buffer,PyBUF_C_CONTIGUOUS|PyBUF_FORMAT|PyBUF_WRITABLE)<0) {
            if(!PyErr_Occurred())PyErr_SetString(PyExc_ValueError,"External terms need an array for their sums.");
            goto cleanup;
        }
        if(external_buffer.ndim!=2 || external_buffer.shape[1]!=EXTERNAL_WIDTH || external_buffer.shape[0]<1 || external_buffer.shape[0]>EXTERNAL_MAX ||
           sums_buffer.ndim!=2 || sums_buffer.shape[0]!=external_buffer.shape[0] || sums_buffer.shape[1]!=EXTERNAL_SUMS ||
           external_buffer.itemsize!=sizeof(double) || !external_buffer.format || strcmp(external_buffer.format,"d")!=0 ||
           sums_buffer.itemsize!=sizeof(double) || !sums_buffer.format || strcmp(sums_buffer.format,"d")!=0) {
            PyErr_SetString(PyExc_ValueError,"External terms must be a float64 (terms,5) array with a (terms,3) array for their sums.");goto cleanup;
        }
        terms=(int)external_buffer.shape[0]; external=external_buffer.buf; sums=sums_buffer.buf;
    }
    /* A direction is periodic unless the system is flat there or has walls on it. */
    int periodic[3]={1,1,dimension==3}, everywhere=dimension==3;
    if(dimension!=2 && dimension!=3) {PyErr_SetString(PyExc_ValueError,"The dimension must be 2 or 3.");goto cleanup;}
    for(int t=0;t<terms;++t) {
        const double *e=external+EXTERNAL_WIDTH*t;
        if(!isfinite(e[0]+e[1]+e[2]+e[3]+e[4]) || (e[0]!=0 && e[0]!=1) || (e[1]!=0 && e[1]!=1 && e[1]!=2) || e[1]>=dimension ||
           !(e[2]>0) || e[3]<0 || (e[3]>0 && (e[0]==0 || !(e[4]>0))) || !isfinite(sums[EXTERNAL_SUMS*t+1]+sums[EXTERNAL_SUMS*t+2])) {
            PyErr_SetString(PyExc_ValueError,"Invalid wall or piston.");goto cleanup;
        }
        periodic[(int)e[1]]=0;
    }
    if((dimension==2 || terms) && (model!=0 || targetP>=0)) {
        PyErr_SetString(PyExc_ValueError,"A flat system and walls take unbonded atoms, without pressure control.");goto cleanup;
    }
    double fixed[3];
    for(int d=0;d<3;++d) {
        if(periodic[d] ? !(box[d]>2*cutoff) || !isfinite(box[d]) : d<dimension && (!(box[d]>0) || !isfinite(box[d]))) {
            PyErr_SetString(PyExc_ValueError, periodic[d] ? "Each box length must be greater than twice the cutoff." : "The distance between two walls must be positive.");goto cleanup;
        }
        fixed[d]=periodic[d]?box[d]:NO_IMAGE;
        if(!periodic[d])everywhere=0;
    }
    /* The image distances: the box itself when it is periodic all round (and may
       then change under a barostat), otherwise fixed for the call. */
    const double *period=everywhere?box:fixed;
    for(int i=0;i<n;++i) {
        if (!(mass[i]>0) || !isfinite(mass[i])) {PyErr_SetString(PyExc_ValueError,"Masses must be positive.");goto cleanup;}
        for(int d=0;d<3;++d) if(!isfinite(v[3*i+d])) {
            PyErr_SetString(PyExc_ValueError,"State arrays must be finite.");goto cleanup;
        }
    }
    const double scalars[WORKSPACE_SCALARS]={cutoff,skin,model,pair,exclude,epsilon,sigma,length,k2,k3,k4,dimension};
    size_t parameter_count=parameters?(size_t)2*n:0, bond_count=bond_parameters&&model?(size_t)2*n:0;
    Neighbors *nb=&ws->nb;
    const int continued=unchanged(ws,n,x,u,box,scalars,parameters,parameter_count,bond_parameters,bond_count) &&
                        ws->externals==terms && same(ws->external,external,(size_t)terms*EXTERNAL_WIDTH);
    ws->current=0;
    if(continued) {
        /* Continue the previous call: its list and forces still hold, and the
           positions and parameters are the ones it checked. */
        f=ws->f; potential=ws->potential; virial=ws->virial;
        memcpy(wall,ws->external_force,sizeof wall);
    } else {
        for(int k=0;k<3*n;++k) if(!isfinite(x[k])||!isfinite(u[k])) {
            PyErr_SetString(PyExc_ValueError,"State arrays must be finite.");goto cleanup;
        }
        if(parameters)for(int i=0;i<2*n;++i)if(!(parameters[i]>0)||!isfinite(parameters[i])){
            PyErr_SetString(PyExc_ValueError,"Pair epsilon and sigma must be finite and positive.");goto cleanup;
        }
        if(bond_parameters)for(int i=0;i<(model?n/2:0);++i){const double *b=bond_parameters+4*i;
            if(!(b[0]>0)||b[1]<0||b[3]<0||!isfinite(b[0]+b[1]+b[2]+b[3])){
                PyErr_SetString(PyExc_ValueError,"Invalid bond parameters.");goto cleanup;
            }
        }
        if(prepare_workspace(ws,n)<0)goto cleanup;
        f=ws->f;
        if(build_neighbors(nb,n,x,u,box,period,dimension,periodic,cutoff+skin,-1)<0 || forces(nb,n,x,period,cutoff,model,pair,epsilon,sigma,parameters,exclude,length,k2,k3,k4,bond_parameters,f,&potential,&virial)<0) goto cleanup;
        external_forces(terms,external,n,x,box,f,&potential,wall);
    }
    uint64_t rng = seed ? seed : UINT64_C(88172645463325252);
    double dof=3*n-3-(model==2?n/2:0), pressure_numerator=0;
    /* Symmetric isotropic pressure-control splitting. rate = d(log L)/dt.
       With gamma=heat=0 the conserved extended enthalpy is
       K + U + targetP*V + baromass*rate^2/2. The 3/dof terms account
       for the finite-size, centre-of-mass-removed kinetic degrees of freedom. */
    while(done<requested) {
        if(targetP>=0){
            double K=kinetic(n,v,mass),V=box[0]*box[1]*box[2];
            rate+=0.5*dt*((2*K+virial-3*targetP*V)+6*K/dof)/baromass;
            double factor=exp(-0.5*dt*(1+3/dof)*rate);
            for(int k=0;k<3*n;++k)v[k]*=factor;
        }
        kick(n,v,f,mass,0.5*dt);
        double constraint_impulse=0;
        if(targetP>=0){
            double factor=exp(0.5*dt*rate);
            if(!isfinite(factor)||factor<0.9||factor>1.1){PyErr_SetString(PyExc_ValueError,"Unstable barostat; increase pressure damping or reduce timestep.");goto cleanup;}
            for(int d=0;d<3;++d)box[d]*=factor;
            for(int k=0;k<3*n;++k){x[k]*=factor;u[k]*=factor;}
        }
        if(gamma>0) {
            if(model==2){if(rigid_drift(n,x,u,v,box,mass,0.5*dt,length,bond_parameters,&constraint_impulse)<0)goto cleanup;}
            else if(drift(n,x,u,v,box,periodic,0.5*dt)<0)goto cleanup;
            thermostat(n,v,mass,T,gamma,dt,&rng,dimension,periodic);
            if(model==2){constrain_velocities(n,v,x,box,mass);if(rigid_drift(n,x,u,v,box,mass,0.5*dt,length,bond_parameters,&constraint_impulse)<0)goto cleanup;}
            else if(drift(n,x,u,v,box,periodic,0.5*dt)<0)goto cleanup;
        } else if(model==2){if(rigid_drift(n,x,u,v,box,mass,dt,length,bond_parameters,&constraint_impulse)<0)goto cleanup;}
        else if(drift(n,x,u,v,box,periodic,dt)<0)goto cleanup;
        if(targetP>=0){
            double factor=exp(0.5*dt*rate);
            for(int d=0;d<3;++d){box[d]*=factor;if(!(box[d]>2*cutoff)){PyErr_SetString(PyExc_ValueError,"Box contracted below twice the cutoff.");goto cleanup;}}
            for(int k=0;k<3*n;++k){x[k]*=factor;u[k]*=factor;}
        }
        /* The pistons move with the drift, so the forces that follow are those at
           the new wall positions. Over the step each term receives the impulse
           dt*(force before + force after)/2, the negative of what the two half
           kicks give the atoms, and a wall that moved has done the work of minus
           that mean force over its displacement. */
        if(terms){memcpy(before,wall,sizeof wall);move_pistons(terms,external,box,dt,moved);}
        if(needs_rebuild(nb,n,u,period,cutoff,skin) && build_neighbors(nb,n,x,u,box,period,dimension,periodic,cutoff+skin,-1)<0)goto cleanup;
        if(forces(nb,n,x,period,cutoff,model,pair,epsilon,sigma,parameters,exclude,length,k2,k3,k4,bond_parameters,f,&potential,&virial)<0)goto cleanup;
        if(terms) {
            external_forces(terms,external,n,x,box,f,&potential,wall);
            for(int t=0;t<terms;++t) {
                double mean=0.5*(before[t]+wall[t]);
                sums[EXTERNAL_SUMS*t+1]+=dt*mean; sums[EXTERNAL_SUMS*t+2]-=mean*moved[t];
            }
        }
        kick(n,v,f,mass,0.5*dt);
        if(model==2)constrain_velocities(n,v,x,box,mass);
        if(targetP>=0){
            double factor=exp(-0.5*dt*(1+3/dof)*rate);
            for(int k=0;k<3*n;++k)v[k]*=factor;
            double K=kinetic(n,v,mass),V=box[0]*box[1]*box[2];
            rate+=0.5*dt*((2*K+virial-3*targetP*V)+6*K/dof)/baromass;
        }
        if(heat!=0){
            double K=kinetic(n,v,mass);
            if(K<=0 || K+heat*dt<=0){PyErr_SetString(PyExc_ValueError,"Heat removal would exhaust kinetic energy.");goto cleanup;}
            double factor=sqrt((K+heat*dt)/K);
            for(int k=0;k<3*n;++k)v[k]*=factor;
        }
        ++done;
        if ((done%64==0 || done==requested) && PyErr_CheckSignals()<0) {
            if(!PyErr_ExceptionMatches(PyExc_KeyboardInterrupt))goto cleanup;
            PyErr_Clear(); interrupted=1; break;
        }
    }
    /* The molecular pressure takes a pass over the pairs; measure=0 leaves it for a later call of zero steps. */
    pressure_numerator=!model?2*kinetic(n,v,mass)+virial:measure?molecular_pressure_numerator(nb,n,x,v,mass,box,cutoff,pair,epsilon,sigma,parameters):NAN;
    result=Py_BuildValue("(ddKiidd)",potential,virial,(unsigned long long)rng,done,interrupted,rate,pressure_numerator);
    if(result) {
        memcpy(buffers[3].buf,f,(size_t)n*3*sizeof(double));
        for(int t=0;t<terms;++t)sums[EXTERNAL_SUMS*t]=wall[t];
        if(ws!=&local) {
            remember(ws,n,x,u,box,scalars,!continued,parameters,parameter_count,bond_parameters,bond_count,potential,virial);
            ws->externals=terms;
            if(terms)memcpy(ws->external,external,(size_t)terms*EXTERNAL_WIDTH*sizeof(double));
            memcpy(ws->external_force,wall,sizeof wall);
        }
    }

cleanup:
    if(sums_buffer.obj)PyBuffer_Release(&sums_buffer);
    if(external_buffer.obj)PyBuffer_Release(&external_buffer);
    if(bond_buffer.obj)PyBuffer_Release(&bond_buffer);
    if(parameter_buffer.obj)PyBuffer_Release(&parameter_buffer);
    release_workspace(&local);
    for(int k=0;k<6;++k) if(buffers[k].obj) PyBuffer_Release(&buffers[k]);
    return result;
}
/* Isolated all-pairs Plummer-softened gravity. No cutoff, images or walls.
   U_ij=-G*m_i*m_j/sqrt(r^2+a^2); W=sum r_ij dot F_ij. */
static int gravity_forces(int n, const double *x, const double *mass,
                          double G, double softening, double *f, double *U, double *W) {
    memset(f,0,(size_t)3*n*sizeof(double));
    *U=0; *W=0;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
        double dr[3],r2=0;
        for(int d=0;d<3;++d){dr[d]=x[3*i+d]-x[3*j+d];r2+=dr[d]*dr[d];}
        if(r2 == 0 && softening == 0){PyErr_SetString(PyExc_ValueError,"Coincident particles make unsoftened gravity singular.");return -1;}
        double inv=1/sqrt(r2+softening*softening);
        double energy=-G*mass[i]*mass[j]*inv, coefficient=energy*inv*inv;
        *U+=energy; *W+=coefficient*r2;
        for(int d=0;d<3;++d){double fd=coefficient*dr[d];f[3*i+d]+=fd;f[3*j+d]-=fd;}
    }
    if(!isfinite(*U)||!isfinite(*W)){
        PyErr_SetString(PyExc_ValueError,"Nonfinite gravitational energy; reduce timestep or check parameters.");return -1;
    }
    for(int k=0;k<3*n;++k)if(!isfinite(f[k])){
        PyErr_SetString(PyExc_ValueError,"Nonfinite gravitational force.");return -1;
    }
    return 0;
}

static PyObject *advance_gravity(PyObject *self, PyObject *args) {
    (void)self;
    PyObject *objects[4], *result=NULL, *workspace_object=NULL;
    Py_buffer b[4]={{0}};
    Workspace local={0}, *ws;
    int requested,done=0,interrupted=0;
    double dt,G,a,heat,U,W;
    if(!PyArg_ParseTuple(args,"OOOOidddd|O",&objects[0],&objects[1],&objects[2],&objects[3],
                         &requested,&dt,&G,&a,&heat,&workspace_object))return NULL;
    if(!(ws=borrowed_workspace(workspace_object,&local)))return NULL;
    for(int k=0;k<4;++k){
        if(PyObject_GetBuffer(objects[k],&b[k],PyBUF_C_CONTIGUOUS|PyBUF_FORMAT|(k<3?PyBUF_WRITABLE:0))<0)goto cleanup;
        if(b[k].itemsize!=sizeof(double)||!b[k].format||strcmp(b[k].format,"d")){
            PyErr_SetString(PyExc_TypeError,"Gravity arrays must be contiguous native float64 buffers.");goto cleanup;
        }
    }
    if(b[0].ndim!=2||b[0].shape[1]!=3||b[0].shape[0]<2||b[0].shape[0]>4096){
        PyErr_SetString(PyExc_ValueError,"Direct gravity requires positions (N,3), 2 <= N <= 4096.");goto cleanup;
    }
    int n=(int)b[0].shape[0];
    for(int k=1;k<3;++k)if(b[k].ndim!=2||b[k].shape[0]!=n||b[k].shape[1]!=3){
        PyErr_SetString(PyExc_ValueError,"Gravity state shapes must match positions.");goto cleanup;
    }
    if(b[3].ndim!=1||b[3].shape[0]!=n){PyErr_SetString(PyExc_ValueError,"Mass shape must be (N,).");goto cleanup;}
    for(int i=0;i<4;++i)for(int j=i+1;j<4;++j){
        uintptr_t p=(uintptr_t)b[i].buf,q=(uintptr_t)b[j].buf;
        if(p<q+b[j].len&&q<p+b[i].len){PyErr_SetString(PyExc_ValueError,"Gravity buffers must not overlap.");goto cleanup;}
    }
    if(requested<0||!(dt>0)||!(G>0)||a<0||!isfinite(dt+G+a+heat)){
        PyErr_SetString(PyExc_ValueError,"Gravity requires positive finite timestep and G, nonnegative softening, and finite heat_rate.");goto cleanup;
    }
    double *x=b[0].buf,*v=b[1].buf,*f=NULL,*mass=b[3].buf;
    for(int i=0;i<n;++i){
        if(!(mass[i]>0)||!isfinite(mass[i])){PyErr_SetString(PyExc_ValueError,"Masses must be finite and positive.");goto cleanup;}
        for(int d=0;d<3;++d)if(!isfinite(x[3*i+d])||!isfinite(v[3*i+d])){
            PyErr_SetString(PyExc_ValueError,"Gravity state must be finite.");goto cleanup;
        }
    }
    /* There is no neighbour list here; the workspace keeps the forces between calls. */
    const double scalars[WORKSPACE_SCALARS]={G,a}, origin[3]={0,0,0};
    const int continued=unchanged(ws,n,x,NULL,origin,scalars,mass,(size_t)n,NULL,0);
    ws->current=0;
    if(continued){
        f=ws->f; U=ws->potential; W=ws->virial;
    } else {
        if(prepare_workspace(ws,n)<0)goto cleanup;
        f=ws->f;
        if(gravity_forces(n,x,mass,G,a,f,&U,&W)<0)goto cleanup;
    }
    while(done<requested){
        kick(n,v,f,mass,.5*dt);
        for(int k=0;k<3*n;++k){
            x[k]+=dt*v[k];
            if(!isfinite(x[k])){PyErr_SetString(PyExc_ValueError,"Unstable gravity step.");goto cleanup;}
        }
        if(gravity_forces(n,x,mass,G,a,f,&U,&W)<0)goto cleanup;
        kick(n,v,f,mass,.5*dt);
        if(heat!=0){
            double K=kinetic(n,v,mass);
            if(K<=0||K+heat*dt<=0){PyErr_SetString(PyExc_ValueError,"Heat removal would exhaust kinetic energy.");goto cleanup;}
            double factor=sqrt((K+heat*dt)/K);
            for(int k=0;k<3*n;++k)v[k]*=factor;
        }
        ++done;
        if((done%64==0||done==requested)&&PyErr_CheckSignals()<0){
            if(!PyErr_ExceptionMatches(PyExc_KeyboardInterrupt))goto cleanup;
            PyErr_Clear();interrupted=1;break;
        }
    }
    result=Py_BuildValue("ddii",U,W,done,interrupted);
    if(result){
        memcpy(b[2].buf,f,(size_t)n*3*sizeof(double));
        if(ws!=&local)remember(ws,n,x,NULL,origin,scalars,!continued,mass,(size_t)n,NULL,0,U,W);
    }
cleanup:
    release_workspace(&local);
    for(int k=0;k<4;++k)if(b[k].obj)PyBuffer_Release(&b[k]);
    return result;
}

/* Test helper: the (i, j) pairs within radius, by checking all pairs (reach 0) or with cells.
   Optionally the dimension (z is ignored in 2) and which directions are periodic. */
static PyObject *neighbor_pairs(PyObject *self, PyObject *args) {
    (void)self;
    PyObject *object, *result = NULL;
    Py_buffer view = {0};
    Workspace w = {0};
    double box[3], period[3], radius;
    int reach, dimension = 3, periodic[3] = {1, 1, 1};
    if (!PyArg_ParseTuple(args, "O(ddd)di|i(iii)", &object, &box[0], &box[1], &box[2], &radius, &reach,
                          &dimension, &periodic[0], &periodic[1], &periodic[2])) return NULL;
    if (PyObject_GetBuffer(object, &view, PyBUF_C_CONTIGUOUS | PyBUF_FORMAT) < 0) return NULL;
    if (view.itemsize != sizeof(double) || !view.format || strcmp(view.format, "d") != 0 || view.ndim != 2 ||
        view.shape[1] != 3 || view.shape[0] < 2 || view.shape[0] > 20000 || !(radius > 0) || reach < 0 || reach > 2 ||
        (dimension != 2 && dimension != 3) || !(box[0] > 0) || !(box[1] > 0) || (dimension == 3 && !(box[2] > 0))) {
        PyErr_SetString(PyExc_ValueError, "Positions must be float64 (N,3), with a positive box and radius, and reach 0, 1 or 2.");
        goto cleanup;
    }
    int n = (int)view.shape[0];
    for (int d = 0; d < 3; ++d) {
        if (d >= dimension) periodic[d] = 0;
        period[d] = periodic[d] ? box[d] : NO_IMAGE;
    }
    if (prepare_workspace(&w, n) < 0 || build_neighbors(&w.nb, n, view.buf, view.buf, box, period, dimension, periodic, radius, reach) < 0) goto cleanup;
    result = PyBytes_FromStringAndSize((const char *)w.nb.pairs, (Py_ssize_t)(w.nb.size * sizeof(Pair)));
cleanup:
    release_workspace(&w);
    PyBuffer_Release(&view);
    return result;
}

static PyMethodDef methods[] = {
    {"workspace", workspace, METH_NOARGS, "What one run keeps between calls: the neighbour list and the state the solver last returned."},
    {"neighbor_pairs", neighbor_pairs, METH_VARARGS, "Pairs within a radius as int (i, j) bytes, from all pairs or from cells of reach 1 or 2 (test helper)."},
    {"advance_gravity", advance_gravity, METH_VARARGS, "Advance isolated softened gravity in-place."},
    {"advance", advance, METH_VARARGS, "Advance LJ dynamics in-place, in three dimensions or two, with optional walls and pistons; return energy, virial, RNG, steps, interrupt flag, barostat rate and pressure numerator."},
    {"molecular_workspace", molecular_workspace, METH_NOARGS, "What one run of charged atoms or water keeps between calls: molecule list, PPPM mesh and forces."},
    {"advance_molecular", (PyCFunction)(void(*)(void))advance_molecular, METH_VARARGS|METH_KEYWORDS, "Advance charged atoms or rigid water in-place; return a dict of energies and state."},
    {"project_molecular_velocities", project_molecular_velocities, METH_VARARGS, "Remove velocity components along rigid-water constraints in-place."},
    {"pppm_parameters", pppm_parameters, METH_VARARGS, "Choose PPPM mesh and Ewald parameter from a relative accuracy."},
    {"fft_transform", fft_transform, METH_VARARGS, "In-place unnormalized 3D FFT of a complex128 array (test helper)."},
    {NULL,NULL,0,NULL}
};
static struct PyModuleDef module = {PyModuleDef_HEAD_INIT,"_core",NULL,-1,methods};
PyMODINIT_FUNC PyInit__core(void) {return PyModule_Create(&module);}
