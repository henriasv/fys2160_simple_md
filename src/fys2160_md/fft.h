/* Small mixed-radix (2, 3, 5) complex FFT for PPPM meshes. Data are
   interleaved (re, im) doubles; transforms are unnormalized. */
#ifndef FYS2160_FFT_H
#define FYS2160_FFT_H

typedef struct {
    int n;
    int factors[64];      /* pairs (p, m) for each recursion level */
    double *twiddle;      /* exp(-2*pi*i*j/n), j = 0..n-1 */
    double *scratch;      /* 2*n doubles for one line */
} FFTPlan;

typedef struct {
    int n[3];
    FFTPlan plan[3];
    double *line;         /* 2*max(n) doubles, gathered strided lines */
} FFT3D;

int fft_factorable(int n);
int fft3d_init(FFT3D *fft, const int n[3]);
void fft3d_free(FFT3D *fft);
/* sign=-1: forward exp(-ikx); sign=+1: inverse exp(+ikx), not normalized. */
void fft3d(FFT3D *fft, double *data, int sign);

#endif
