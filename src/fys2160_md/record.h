/* Per-step copies for measurements. After step k of a call the solver can copy
   the positions, velocities, unwrapped positions and box into row k of caller
   arrays: (steps, N, 3) for the first three and (steps, 3) for the box. Without
   them nothing is copied and the step loop is unchanged. */
#ifndef FYS2160_RECORD_H
#define FYS2160_RECORD_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <string.h>

typedef struct { Py_buffer view[4]; double *data[4]; } Record;

/* record is NULL, None, or a tuple (positions, velocities, unwrapped, box), each
   None or a writable C-contiguous float64 array with at least `steps` rows. */
static int record_open(Record *r, PyObject *record, int n, int steps) {
    if (!record || record == Py_None) return 0;
    if (!PyTuple_Check(record) || PyTuple_GET_SIZE(record) != 4) {
        PyErr_SetString(PyExc_ValueError, "record must be None or a tuple of four arrays (or None).");
        return -1;
    }
    for (int k = 0; k < 4; ++k) {
        PyObject *object = PyTuple_GET_ITEM(record, k);
        if (object == Py_None) continue;
        Py_buffer *view = &r->view[k];
        if (PyObject_GetBuffer(object, view, PyBUF_C_CONTIGUOUS | PyBUF_FORMAT | PyBUF_WRITABLE) < 0) return -1;
        int ndim = k < 3 ? 3 : 2;
        if (view->itemsize != sizeof(double) || !view->format || strcmp(view->format, "d") != 0 ||
            view->ndim != ndim || view->shape[0] < steps || view->shape[ndim-1] != 3 || (k < 3 && view->shape[1] != n)) {
            PyErr_SetString(PyExc_ValueError, "Record arrays must be float64 (steps, N, 3), and (steps, 3) for the box.");
            return -1;
        }
        r->data[k] = view->buf;
    }
    return 0;
}

static void record_step(Record *r, int k, int n, const double *x, const double *v, const double *u, const double *box) {
    const double *from[4] = {x, v, u, box};
    for (int a = 0; a < 4; ++a) if (r->data[a] && from[a]) {
        size_t count = a < 3 ? (size_t)3*(size_t)n : 3;
        memcpy(r->data[a] + (size_t)k*count, from[a], count*sizeof(double));
    }
}

static void record_close(Record *r) {
    for (int k = 0; k < 4; ++k) if (r->view[k].obj) PyBuffer_Release(&r->view[k]);
}
#endif
