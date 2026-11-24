/* A tilted (triclinic) box, as in LAMMPS: edges a = (lx, 0, 0), b = (xy, ly, 0) and
   c = (xz, yz, lz), with tilt = {xy, xz, yz}. A NULL or zero tilt is a rectangular
   box, for which callers keep their per-axis code. Shared by the pair solver, the
   charged-atom solver and PPPM. */
#ifndef FYS2160_TILT_H
#define FYS2160_TILT_H
#include <math.h>

static const double NO_TILT[3] = {0, 0, 0};
/* The nearest image takes whole edges off z, then y, then x: it rounds the
   fractional coordinates, which finds the nearest image of any separation shorter
   than half the smallest distance between opposite faces. */
static inline void tilted_image(double dr[3], const double *box, const double *tilt) {
    double k = nearbyint(dr[2]/box[2]);
    dr[2] -= k*box[2]; dr[1] -= k*tilt[2]; dr[0] -= k*tilt[1];
    k = nearbyint(dr[1]/box[1]);
    dr[1] -= k*box[1]; dr[0] -= k*tilt[0];
    dr[0] -= box[0]*nearbyint(dr[0]/box[0]);
}
/* Fractional coordinates s, with x = s0 a + s1 b + s2 c. */
static inline void fractional(const double *x, const double *box, const double *tilt, double *s) {
    s[2] = x[2]/box[2];
    s[1] = (x[1] - tilt[2]*s[2])/box[1];
    s[0] = (x[0] - tilt[0]*s[1] - tilt[1]*s[2])/box[0];
}
/* Moves a position by whole edges into the box, the parallelepiped 0 <= s < 1. */
static inline void tilted_wrap(double *x, const double *box, const double *tilt) {
    double s[3];
    fractional(x, box, tilt, s);
    double k0 = floor(s[0]), k1 = floor(s[1]), k2 = floor(s[2]);
    x[0] -= k0*box[0] + k1*tilt[0] + k2*tilt[1];
    x[1] -= k1*box[1] + k2*tilt[2];
    x[2] -= k2*box[2];
}
/* The distance between the two faces across each axis: volume over face area. */
static inline void face_distances(const double *box, const double *tilt, double *w) {
    double lx = box[0], ly = box[1], lz = box[2], xy = tilt[0], xz = tilt[1], yz = tilt[2], V = lx*ly*lz;
    w[0] = V/sqrt(ly*lz*ly*lz + xy*lz*xy*lz + (xy*yz - ly*xz)*(xy*yz - ly*xz));
    w[1] = V/(lx*sqrt(lz*lz + yz*yz));
    w[2] = lz;
}
static inline int tilted(const double *tilt) { return tilt && (tilt[0] != 0 || tilt[1] != 0 || tilt[2] != 0); }

#endif
