/* Recursive decimation-in-time FFT for lengths with factors 2, 3 and 5. */
#include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int fft_factorable(int n) {
    if (n < 1) return 0;
    while (n % 2 == 0) n /= 2;
    while (n % 3 == 0) n /= 3;
    while (n % 5 == 0) n /= 5;
    return n == 1;
}

static void plan_free(FFTPlan *p) {
    free(p->twiddle); free(p->scratch);
    p->twiddle = p->scratch = NULL;
}

static int plan_init(FFTPlan *p, int n) {
    memset(p, 0, sizeof *p);
    if (!fft_factorable(n)) return -1;
    p->n = n;
    p->twiddle = malloc(2 * (size_t)n * sizeof(double));
    p->scratch = malloc(2 * (size_t)n * sizeof(double));
    if (!p->twiddle || !p->scratch) { plan_free(p); return -1; }
    for (int j = 0; j < n; ++j) {
        double angle = 6.283185307179586476925 * j / n;
        p->twiddle[2*j] = cos(angle); p->twiddle[2*j+1] = sin(angle);
    }
    int k = 0, m = n;
    while (m > 1) {
        int radix = m % 2 == 0 ? 2 : m % 3 == 0 ? 3 : 5;
        m /= radix;
        p->factors[k++] = radix; p->factors[k++] = m;
    }
    return 0;
}

/* out holds radix sub-transforms of length m: out[q*m+u] = Y_q[u]. */
static void butterfly(double *out, int fstride, const FFTPlan *p, int m, int radix, int sign) {
    const double *tw = p->twiddle;
    const int n = p->n;
    if (radix == 2) {
        for (int u = 0; u < m; ++u) {
            int t = fstride * u;
            double c = tw[2*t], s = sign * tw[2*t+1];
            double *a = out + 2*u, *b = out + 2*(u+m);
            double br = b[0]*c - b[1]*s, bi = b[0]*s + b[1]*c;
            b[0] = a[0] - br; b[1] = a[1] - bi;
            a[0] += br; a[1] += bi;
        }
        return;
    }
    double y[10];
    for (int u = 0; u < m; ++u) {
        for (int q = 0; q < radix; ++q) { y[2*q] = out[2*(u+q*m)]; y[2*q+1] = out[2*(u+q*m)+1]; }
        for (int q1 = 0; q1 < radix; ++q1) {
            int k = u + q1 * m;
            double re = y[0], im = y[1];
            for (int q = 1; q < radix; ++q) {
                int t = (int)(((long long)fstride * q * k) % n);
                double c = tw[2*t], s = sign * tw[2*t+1];
                re += y[2*q]*c - y[2*q+1]*s;
                im += y[2*q]*s + y[2*q+1]*c;
            }
            out[2*k] = re; out[2*k+1] = im;
        }
    }
}

static void transform(double *out, const double *in, int fstride, const int *factors, const FFTPlan *p, int sign) {
    int radix = factors[0], m = factors[1];
    if (m == 1) {
        for (int q = 0; q < radix; ++q) {
            out[2*q] = in[2*q*fstride]; out[2*q+1] = in[2*q*fstride+1];
        }
    } else {
        for (int q = 0; q < radix; ++q)
            transform(out + 2*q*m, in + 2*q*fstride, fstride*radix, factors + 2, p, sign);
    }
    butterfly(out, fstride, p, m, radix, sign);
}

int fft3d_init(FFT3D *fft, const int n[3]) {
    memset(fft, 0, sizeof *fft);
    int largest = 1;
    for (int d = 0; d < 3; ++d) {
        fft->n[d] = n[d];
        if (plan_init(&fft->plan[d], n[d]) < 0) { fft3d_free(fft); return -1; }
        if (n[d] > largest) largest = n[d];
    }
    fft->line = malloc(2 * (size_t)largest * sizeof(double));
    if (!fft->line) { fft3d_free(fft); return -1; }
    return 0;
}

void fft3d_free(FFT3D *fft) {
    for (int d = 0; d < 3; ++d) plan_free(&fft->plan[d]);
    free(fft->line); fft->line = NULL;
}

static void lines(FFT3D *fft, double *data, int axis, int sign) {
    const int nx = fft->n[0], ny = fft->n[1], nz = fft->n[2];
    const FFTPlan *p = &fft->plan[axis];
    if (p->n == 1) return;
    int count, stride;
    if (axis == 2) { count = nx*ny; stride = 1; }
    else if (axis == 1) { count = nx*nz; stride = nz; }
    else { count = ny*nz; stride = ny*nz; }
    for (int c = 0; c < count; ++c) {
        size_t base = axis == 2 ? (size_t)c*nz : axis == 1 ? (size_t)(c/nz)*ny*nz + c%nz : (size_t)c;
        for (int j = 0; j < p->n; ++j) {
            fft->line[2*j] = data[2*(base + (size_t)j*stride)];
            fft->line[2*j+1] = data[2*(base + (size_t)j*stride)+1];
        }
        transform(p->scratch, fft->line, 1, p->factors, p, sign);
        for (int j = 0; j < p->n; ++j) {
            data[2*(base + (size_t)j*stride)] = p->scratch[2*j];
            data[2*(base + (size_t)j*stride)+1] = p->scratch[2*j+1];
        }
    }
}

void fft3d(FFT3D *fft, double *data, int sign) {
    lines(fft, data, 2, sign);
    lines(fft, data, 1, sign);
    lines(fft, data, 0, sign);
}
