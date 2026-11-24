/* PPPM written from Hockney & Eastwood (1988) and Deserno & Holm, J. Chem. Phys.
   109, 7678 and 7694 (1998). Parameter selection follows the LAMMPS kspace pppm
   documentation so meshes and splitting parameters agree for cross-checks. */
#include "pppm.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "tilt.h"

#define PI 3.14159265358979323846
#define MAX_ORDER 7

/* Deserno-Holm coefficients for the ik-differentiated RMS force error. */
static const double ERROR_COEFFICIENTS[MAX_ORDER+1][MAX_ORDER] = {
    {0},
    {2.0/3.0},
    {1.0/50.0, 5.0/294.0},
    {1.0/588.0, 7.0/1440.0, 21.0/3872.0},
    {1.0/4320.0, 3.0/1936.0, 7601.0/2271360.0, 143.0/28800.0},
    {1.0/23232.0, 7601.0/13628160.0, 143.0/69120.0, 517231.0/106536960.0, 106640677.0/11737571328.0},
    {691.0/68140800.0, 13.0/57600.0, 47021.0/35512320.0, 9694607.0/2095994880.0,
     733191589.0/59609088000.0, 326190917.0/11700633600.0},
    {1.0/345600.0, 3617.0/35512320.0, 745739.0/838397952.0, 56399353.0/12773376000.0,
     25091609.0/1560084480.0, 1755948832039.0/36229939200000.0, 4887769399.0/37838389248.0},
};

static double kspace_error(double h, double length, int atoms, double q2, double g, int order) {
    double sum = 0;
    for (int m = 0; m < order; ++m) sum += ERROR_COEFFICIENTS[order][m] * pow(h*g, 2.0*m);
    return q2 * pow(h*g, order) * sqrt(g*length*sqrt(2*PI)*sum/atoms) / (length*length);
}

static double kspace_error_mesh(const double box[3], const int mesh[3], int atoms, double q2, double g, int order) {
    double sum = 0;
    for (int d = 0; d < 3; ++d) {
        double e = kspace_error(box[d]/mesh[d], box[d], atoms, q2, g, order);
        sum += e*e;
    }
    return sqrt(sum/3);
}

/* lengths: lx, ly, lz, whose product is the volume of a tilted box too. */
static double rspace_error(const double lengths[3], int atoms, double q2, double g, double cutoff) {
    return 2*q2*exp(-g*g*cutoff*cutoff) / sqrt(atoms*cutoff*lengths[0]*lengths[1]*lengths[2]);
}

/* box holds the widths across each axis: the lengths of a rectangular box, the
   distances between opposite faces of a tilted one, for which the mesh spacing
   along each reciprocal direction is width/mesh. */
int pppm_select(const double box[3], const double lengths[3], int atoms, double qsqsum, double qqr2e, double cutoff,
                double accuracy, int order, int mesh[3], double *g_out, double *estimate) {
    if (order < 2 || order > MAX_ORDER || atoms < 1 || !(qsqsum > 0) || !(accuracy > 0)) return -1;
    double tolerance = accuracy * qqr2e, q2 = qsqsum * qqr2e;
    double volume = lengths[0]*lengths[1]*lengths[2];
    double g = tolerance * sqrt(atoms*cutoff*volume) / (2*q2);
    g = g >= 1 ? (1.35 - 0.15*log(tolerance))/cutoff : sqrt(-log(g))/cutoff;
    for (int d = 0; d < 3; ++d) {
        double h = 1/g;
        int n = (int)(box[d]/h) + 1;
        double error = kspace_error(h, box[d], atoms, q2, g, order);
        while (error > tolerance) {
            error = kspace_error(h, box[d], atoms, q2, g, order);
            ++n;
            h = box[d]/n;
            if (n > 4096) return -2;
        }
        while (!fft_factorable(n)) ++n;
        mesh[d] = n;
    }
    /* Balance real- and k-space errors by Newton iteration on g. */
    for (int iteration = 0; iteration < 10000; ++iteration) {
        double f = rspace_error(lengths, atoms, q2, g, cutoff) - kspace_error_mesh(box, mesh, atoms, q2, g, order);
        double step = 1e-6;
        double f2 = rspace_error(lengths, atoms, q2, g+step, cutoff) - kspace_error_mesh(box, mesh, atoms, q2, g+step, order);
        g -= f / ((f2 - f)/step);
        if (!(g > 0) || !isfinite(g)) return -3;
        f = rspace_error(lengths, atoms, q2, g, cutoff) - kspace_error_mesh(box, mesh, atoms, q2, g, order);
        if (fabs(f) < 1e-5) {
            double r = rspace_error(lengths, atoms, q2, g, cutoff), k = kspace_error_mesh(box, mesh, atoms, q2, g, order);
            *g_out = g;
            *estimate = sqrt(r*r + k*k) / qqr2e;
            return 0;
        }
    }
    return -3;
}

void pppm_free(PPPM *p) {
    free(p->greens); free(p->rho); free(p->work); free(p->weights); free(p->start);
    for (int d = 0; d < 3; ++d) free(p->k[d]);
    fft3d_free(&p->fft);
    memset(p, 0, sizeof *p);
}

int pppm_init(PPPM *p, const int mesh[3], double g, int order, int sites) {
    memset(p, 0, sizeof *p);
    if (order < 2 || order > MAX_ORDER || !(g > 0) || sites < 0) return -1;
    size_t points = 1;
    for (int d = 0; d < 3; ++d) {
        if (mesh[d] < 2 || !fft_factorable(mesh[d])) return -1;
        p->n[d] = mesh[d];
        /* The meshes hold 2*points doubles; size_t is 32-bit in WebAssembly. */
        if (points > SIZE_MAX/(2*sizeof(double))/(size_t)mesh[d]) { pppm_free(p); return -1; }
        points *= (size_t)mesh[d];
        p->k[d] = malloc((size_t)mesh[d]*sizeof(double));
        if (!p->k[d]) { pppm_free(p); return -1; }
    }
    p->order = order; p->g = g; p->sites = sites;
    p->greens = malloc(points*sizeof(double));
    p->rho = malloc(2*points*sizeof(double));
    p->work = malloc(2*points*sizeof(double));
    p->weights = malloc((size_t)(sites ? sites : 1)*3*order*sizeof(double));
    p->start = malloc((size_t)(sites ? sites : 1)*3*sizeof(int));
    if (!p->greens || !p->rho || !p->work || !p->weights || !p->start || fft3d_init(&p->fft, mesh) < 0) {
        pppm_free(p); return -1;
    }
    return 0;
}

/* sum over aliases m of (sin x/(x + pi m))^(2P), x = pi k/n. */
static double alias_sum(double x, int order) {
    if (fabs(sin(x)) < 1e-300) return 1.0;
    double total = 0, s = sin(x);
    for (int m = -50; m <= 50; ++m) total += pow(s/(x + PI*m), 2*order);
    return total;
}

static double sinc_power(double x, int power) {
    return x == 0 ? 1.0 : pow(sin(x)/x, power);
}

static int wave_index(int i, int n) { return i - n*((2*i)/n); }

/* Optimal ik influence function G(k) = 4 pi/k^2 sum_m (k.k_m/k_m^2) exp(-k_m^2/4g^2) U^2(k_m)
   / (sum_m U^2(k_m))^2, with the LAMMPS alias truncation for the numerator. */
static int update_greens(PPPM *p, const double box[3]) {
    const int *n = p->n, order = p->order;
    const double g = p->g;
    int aliases[3];
    double unit[3], *wave[3] = {0}, *factor[3] = {0}, *denominator[3] = {0};
    int status = -1;
    for (int d = 0; d < 3; ++d) {
        unit[d] = 2*PI/box[d];
        aliases[d] = (int)((g*box[d]/(PI*n[d])) * pow(-log(1e-7), 0.25));
        int width = 2*aliases[d] + 1;
        wave[d] = malloc((size_t)n[d]*width*sizeof(double));
        factor[d] = malloc((size_t)n[d]*width*sizeof(double));
        denominator[d] = malloc((size_t)n[d]*sizeof(double));
        if (!wave[d] || !factor[d] || !denominator[d]) goto done;
        for (int i = 0; i < n[d]; ++i) {
            int w = wave_index(i, n[d]);
            p->k[d][i] = unit[d]*w;
            denominator[d][i] = alias_sum(PI*w/n[d], order);
            for (int a = 0; a < width; ++a) {
                double q = unit[d]*(w + n[d]*(a - aliases[d]));
                wave[d][i*width+a] = q;
                factor[d][i*width+a] = exp(-0.25*(q/g)*(q/g)) * sinc_power(0.5*q*box[d]/n[d], 2*order);
            }
        }
    }
    const int wx = 2*aliases[0]+1, wy = 2*aliases[1]+1, wz = 2*aliases[2]+1;
    size_t index = 0;
    for (int i = 0; i < n[0]; ++i) for (int j = 0; j < n[1]; ++j) for (int l = 0; l < n[2]; ++l, ++index) {
        double k[3] = {p->k[0][i], p->k[1][j], p->k[2][l]};
        double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
        if (k2 == 0) { p->greens[index] = 0; continue; }
        double sum = 0;
        for (int a = 0; a < wx; ++a) {
            double qx = wave[0][i*wx+a], fx = factor[0][i*wx+a];
            for (int b = 0; b < wy; ++b) {
                double qy = wave[1][j*wy+b], fxy = fx*factor[1][j*wy+b];
                for (int c = 0; c < wz; ++c) {
                    double qz = wave[2][l*wz+c];
                    double dot = k[0]*qx + k[1]*qy + k[2]*qz, q2 = qx*qx + qy*qy + qz*qz;
                    sum += dot/q2 * fxy*factor[2][l*wz+c];
                }
            }
        }
        double s = denominator[0][i]*denominator[1][j]*denominator[2][l];
        p->greens[index] = 12.5663706/k2 * sum / (s*s);
    }
    memcpy(p->box, box, 3*sizeof(double));
    memset(p->tilt, 0, sizeof p->tilt);
    status = 0;
done:
    for (int d = 0; d < 3; ++d) { free(wave[d]); free(factor[d]); free(denominator[d]); }
    return status;
}

/* The reciprocal edges of a tilted box: k = 2 pi R w for integer w, where R is the
   inverse transpose of the edge matrix, lower triangular, with rows (R00, 0, 0),
   (R10, R11, 0), (R20, R21, R22). */
static void reciprocal(const double box[3], const double tilt[3], double R[6]) {
    double lx = box[0], ly = box[1], lz = box[2], xy = tilt[0], xz = tilt[1], yz = tilt[2];
    R[0] = 1/lx;
    R[1] = -xy/(lx*ly); R[2] = 1/ly;
    R[3] = (xy*yz - ly*xz)/(lx*ly*lz); R[4] = -yz/(ly*lz); R[5] = 1/lz;
}
static void wave_vector(const double R[6], double w0, double w1, double w2, double k[3]) {
    k[0] = 2*PI*R[0]*w0;
    k[1] = 2*PI*(R[1]*w0 + R[2]*w1);
    k[2] = 2*PI*(R[3]*w0 + R[4]*w1 + R[5]*w2);
}

/* The same influence function in a tilted box. Charges are spread in fractional
   coordinates, so the assignment factors U stay separable in the mesh indices, but
   the Gaussian exp(-q^2/4g^2) of each alias q = 2 pi R (w + n m) is not. */
static int update_greens_tilted(PPPM *p, const double box[3], const double tilt[3]) {
    const int *n = p->n, order = p->order;
    const double g = p->g;
    int aliases[3];
    double R[6], width[3], *factor[3] = {0}, *denominator[3] = {0};
    int status = -1;
    reciprocal(box, tilt, R);
    face_distances(box, tilt, width);
    for (int d = 0; d < 3; ++d) {
        aliases[d] = (int)((g*width[d]/(PI*n[d])) * pow(-log(1e-7), 0.25));
        int span = 2*aliases[d] + 1;
        factor[d] = malloc((size_t)n[d]*span*sizeof(double));
        denominator[d] = malloc((size_t)n[d]*sizeof(double));
        if (!factor[d] || !denominator[d]) goto done;
        for (int i = 0; i < n[d]; ++i) {
            int w = wave_index(i, n[d]);
            denominator[d][i] = alias_sum(PI*w/n[d], order);
            for (int a = 0; a < span; ++a) factor[d][i*span+a] = sinc_power(PI*(w + n[d]*(a - aliases[d]))/n[d], 2*order);
        }
    }
    const int sx = 2*aliases[0]+1, sy = 2*aliases[1]+1, sz = 2*aliases[2]+1;
    size_t index = 0;
    for (int i = 0; i < n[0]; ++i) for (int j = 0; j < n[1]; ++j) for (int l = 0; l < n[2]; ++l, ++index) {
        int w0 = wave_index(i, n[0]), w1 = wave_index(j, n[1]), w2 = wave_index(l, n[2]);
        double k[3];
        wave_vector(R, w0, w1, w2, k);
        double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
        if (k2 == 0) { p->greens[index] = 0; continue; }
        double sum = 0;
        for (int a = 0; a < sx; ++a) {
            double fx = factor[0][i*sx+a];
            for (int b = 0; b < sy; ++b) {
                double fxy = fx*factor[1][j*sy+b];
                for (int c = 0; c < sz; ++c) {
                    double q[3];
                    wave_vector(R, w0 + n[0]*(a - aliases[0]), w1 + n[1]*(b - aliases[1]), w2 + n[2]*(c - aliases[2]), q);
                    double q2 = q[0]*q[0] + q[1]*q[1] + q[2]*q[2], dot = k[0]*q[0] + k[1]*q[1] + k[2]*q[2];
                    sum += dot/q2 * exp(-0.25*q2/(g*g)) * fxy*factor[2][l*sz+c];
                }
            }
        }
        double s = denominator[0][i]*denominator[1][j]*denominator[2][l];
        p->greens[index] = 12.5663706/k2 * sum / (s*s);
    }
    memcpy(p->box, box, 3*sizeof(double));
    memcpy(p->tilt, tilt, 3*sizeof(double));
    status = 0;
done:
    for (int d = 0; d < 3; ++d) { free(factor[d]); free(denominator[d]); }
    return status;
}

/* Centred cardinal B-spline weights M_P(y + j), y in (0, 1], j = 0..P-1. */
static void bspline(double y, int order, double *w) {
    w[0] = 1;
    for (int p = 2; p <= order; ++p) {
        double previous = 0;
        for (int j = 0; j < p; ++j) {
            double current = j < p-1 ? w[j] : 0;
            w[j] = ((y + j)*current + (p - y - j)*previous) / (p - 1);
            previous = current;
        }
    }
}

int pppm_compute(PPPM *p, const double box[3], const double *tilt, int sites, const double *x, const double *q,
                 double qqr2e, double *f, double *energy, double *virial, double *tensor) {
    if (sites != p->sites) return -1;
    const int skew = tilted(tilt);
    if (skew) {
        if ((memcmp(box, p->box, 3*sizeof(double)) != 0 || memcmp(tilt, p->tilt, 3*sizeof(double)) != 0) &&
            update_greens_tilted(p, box, tilt) < 0) return -1;
    } else if ((memcmp(box, p->box, 3*sizeof(double)) != 0 || tilted(p->tilt)) && update_greens(p, box) < 0) return -1;
    const int *n = p->n, order = p->order;
    const size_t points = (size_t)n[0]*n[1]*n[2];
    const double volume = box[0]*box[1]*box[2];
    const double cell = volume / points;
    double qsum = 0, qsqsum = 0, R[6];
    if (skew) reciprocal(box, tilt, R);
    memset(p->rho, 0, 2*points*sizeof(double));
    for (int s = 0; s < sites; ++s) {
        qsum += q[s]; qsqsum += q[s]*q[s];
        double fraction[3];
        if (skew) fractional(x + 3*s, box, tilt, fraction);
        for (int d = 0; d < 3; ++d) {
            double u = (skew ? fraction[d] * n[d] : x[3*s+d] * n[d] / box[d]) - 0.5*order;
            double first = floor(u) + 1;
            bspline(first - u, order, p->weights + (3*s+d)*order);
            int start = (int)fmod(first, (double)n[d]);
            p->start[3*s+d] = start < 0 ? start + n[d] : start;
        }
        if (q[s] == 0) continue;
        const double *wx = p->weights + 3*s*order, *wy = wx + order, *wz = wy + order;
        double scale = q[s] / cell;
        for (int a = 0; a < order; ++a) {
            int i = (p->start[3*s] + a) % n[0];
            for (int b = 0; b < order; ++b) {
                int j = (p->start[3*s+1] + b) % n[1];
                double wxy = scale * wx[a] * wy[b];
                double *row = p->rho + 2*((size_t)i*n[1] + j)*n[2];
                for (int c = 0; c < order; ++c) {
                    int l = (p->start[3*s+2] + c) % n[2];
                    row[2*l] += wxy * wz[c];
                }
            }
        }
    }
    fft3d(&p->fft, p->rho, -1);
    double e = 0, w = 0, inverse = 1.0 / points, sums[6] = {0, 0, 0, 0, 0, 0};
    size_t index = 0;
    for (int i = 0; i < n[0]; ++i) for (int j = 0; j < n[1]; ++j) for (int l = 0; l < n[2]; ++l, ++index) {
        double re = p->rho[2*index], im = p->rho[2*index+1], G = p->greens[index];
        double k[3], kx, ky, kz;
        if (skew) {
            /* The field uses the wave vector without its Nyquist components, which keeps
               it real: the Nyquist index is its own negative. */
            int w0 = wave_index(i, n[0]), w1 = wave_index(j, n[1]), w2 = wave_index(l, n[2]), field[3];
            wave_vector(R, w0, w1, w2, k);
            field[0] = 2*i == n[0] ? 0 : w0; field[1] = 2*j == n[1] ? 0 : w1; field[2] = 2*l == n[2] ? 0 : w2;
            double kf[3];
            wave_vector(R, field[0], field[1], field[2], kf);
            kx = kf[0]; ky = kf[1]; kz = kf[2];
        } else {
            k[0] = p->k[0][i]; k[1] = p->k[1][j]; k[2] = p->k[2][l];
            kx = 2*i == n[0] ? 0 : k[0]; ky = 2*j == n[1] ? 0 : k[1]; kz = 2*l == n[2] ? 0 : k[2];
        }
        double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
        double term = G * (re*re + im*im) * inverse * inverse;
        e += term;
        w += term * (1 - 0.5*k2/(p->g*p->g));
        if (tensor && k2 > 0) {
            /* -dE/d(strain ab) of this wave: E (delta_ab - 2 k_a k_b (1/k^2 + 1/(4 g^2))). */
            double c = -2*term*(1/k2 + 0.25/(p->g*p->g));
            sums[0] += term + c*k[0]*k[0]; sums[1] += term + c*k[1]*k[1]; sums[2] += term + c*k[2]*k[2];
            sums[3] += c*k[0]*k[1]; sums[4] += c*k[0]*k[2]; sums[5] += c*k[1]*k[2];
        }
        /* Field of (d/dx, d/dy) packed as one complex transform: i kx phi - ky phi. */
        double phi_re = G*re*inverse, phi_im = G*im*inverse;
        p->work[2*index] = -kx*phi_im - ky*phi_re;
        p->work[2*index+1] = kx*phi_re - ky*phi_im;
        p->rho[2*index] = -kz*phi_im;
        p->rho[2*index+1] = kz*phi_re;
    }
    *energy = qqr2e * (0.5*volume*e - p->g/sqrt(PI)*qsqsum - 0.5*PI*qsum*qsum/(p->g*p->g*volume));
    *virial = qqr2e * 0.5*volume*w;
    if (tensor) for (int a = 0; a < 6; ++a) tensor[a] = qqr2e * 0.5*volume*sums[a];
    fft3d(&p->fft, p->work, 1);
    fft3d(&p->fft, p->rho, 1);
    for (int s = 0; s < sites; ++s) {
        if (q[s] == 0) continue;
        const double *wx = p->weights + 3*s*order, *wy = wx + order, *wz = wy + order;
        double gradient[3] = {0, 0, 0};
        for (int a = 0; a < order; ++a) {
            int i = (p->start[3*s] + a) % n[0];
            for (int b = 0; b < order; ++b) {
                int j = (p->start[3*s+1] + b) % n[1];
                double wxy = wx[a] * wy[b];
                size_t row = ((size_t)i*n[1] + j)*n[2];
                for (int c = 0; c < order; ++c) {
                    size_t m = row + (p->start[3*s+2] + c) % n[2];
                    double weight = wxy * wz[c];
                    gradient[0] += weight * p->work[2*m];
                    gradient[1] += weight * p->work[2*m+1];
                    gradient[2] += weight * p->rho[2*m];
                }
            }
        }
        for (int d = 0; d < 3; ++d) f[3*s+d] -= qqr2e * q[s] * gradient[d];
    }
    return isfinite(*energy) && isfinite(*virial) ? 0 : -1;
}
