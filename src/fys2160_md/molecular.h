/* Charged atoms and rigid three-site water: Python entry points. */
#ifndef FYS2160_MOLECULAR_H
#define FYS2160_MOLECULAR_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>

PyObject *advance_molecular(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject *project_molecular_velocities(PyObject *self, PyObject *args);
PyObject *pppm_parameters(PyObject *self, PyObject *args);
PyObject *fft_transform(PyObject *self, PyObject *args);

#endif
