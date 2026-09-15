#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int i, j; } Pair;
typedef struct { Pair *pairs; size_t size, capacity; double *reference; double box[3]; } Neighbors;

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
static int build_neighbors(Neighbors *nb, int n, const double *x,
                           const double *u, const double *box, double radius) {
    nb->size = 0;
    double rmax2 = radius * radius;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
        double r2 = 0;
        for (int d = 0; d < 3; ++d) {
            double dx = minimum_image(x[3*i+d] - x[3*j+d], box[d]);
            r2 += dx * dx;
        }
        if (r2 < rmax2) {
            if (nb->size == nb->capacity) {
                size_t cap = nb->capacity ? nb->capacity * 2 : 1024;
                Pair *p = realloc(nb->pairs, cap * sizeof(Pair));
                if (!p) { PyErr_NoMemory(); return -1; }
                nb->pairs = p; nb->capacity = cap;
            }
            nb->pairs[nb->size++] = (Pair){i, j};
        }
    }
    memcpy(nb->reference, u, (size_t)n * 3 * sizeof(double));
    memcpy(nb->box, box, 3*sizeof(double));
    return 0;
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
static int drift(int n, double *x, double *u, const double *v,
                 const double *box, double duration) {
    for (int i = 0; i < n; ++i) for (int d = 0; d < 3; ++d) {
        int k = 3*i+d;
        double dx = duration * v[k];
        u[k] += dx; x[k] += dx;
        if (!isfinite(u[k]) || !isfinite(x[k])) {
            PyErr_SetString(PyExc_FloatingPointError, "Non-finite position; simulation is unstable.");
            return -1;
        }
        x[k] -= box[d] * floor(x[k]/box[d]);
    }
    return 0;
}
static void kick(int n, double *v, const double *f, const double *mass, double dt) {
    for (int i = 0; i < n; ++i) for (int d = 0; d < 3; ++d)
        v[3*i+d] += dt * f[3*i+d] / mass[i];
}
/* The thermostat samples the canonical distribution with total momentum zero. */
static void thermostat(int n, double *v, const double *mass, double T,
                       double gamma, double dt, uint64_t *rng) {
    double c = exp(-gamma * dt), momentum[3] = {0,0,0}, totalmass = 0;
    for (int i = 0; i < n; ++i) {
        double sigma = sqrt(-expm1(-2*gamma*dt) * T / mass[i]);
        totalmass += mass[i];
        for (int d = 0; d < 3; ++d) {
            v[3*i+d] = c*v[3*i+d] + sigma*gaussian(rng);
            momentum[d] += mass[i]*v[3*i+d];
        }
    }
    for (int i = 0; i < n; ++i) for (int d = 0; d < 3; ++d)
        v[3*i+d] -= momentum[d]/totalmass;
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
    return drift(n,x,u,v,box,dt);
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
    PyObject *objects[6], *result = NULL, *parameter_object=NULL, *bond_object=NULL;
    Py_buffer parameter_buffer={0}, bond_buffer={0};
    const double *parameters=NULL, *bond_parameters=NULL;
    Py_buffer buffers[6] = {{0}};
    int requested, done = 0, interrupted = 0;
    double dt, cutoff, skin, T, gamma, heat, targetP, baromass, rate, potential=0, virial=0;
    int model, pair=-1, exclude=1;
    double length=.7, k2=35., k3=-40., k4=59., epsilon=1., sigma=1.;
    unsigned long long seed;
    Neighbors nb = {0};
    if (!PyArg_ParseTuple(args, "OOOOOOidddddKidddd|iiddddddOO", &objects[0], &objects[1],
        &objects[2], &objects[3], &objects[4], &objects[5], &requested,
        &dt, &cutoff, &skin, &T, &gamma, &seed, &model, &heat, &targetP, &baromass, &rate, &pair, &exclude, &length, &k2, &k3, &k4, &epsilon, &sigma, &parameter_object, &bond_object)) return NULL;
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
        for(int i=0;i<2*n;++i)if(!(parameters[i]>0)||!isfinite(parameters[i])){
            PyErr_SetString(PyExc_ValueError,"Pair epsilon and sigma must be finite and positive.");goto cleanup;
        }
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
        for(int i=0;i<count;++i){const double *b=bond_parameters+4*i;
            if(!(b[0]>0)||b[1]<0||b[3]<0||!isfinite(b[0]+b[1]+b[2]+b[3])){
                PyErr_SetString(PyExc_ValueError,"Invalid bond parameters.");goto cleanup;
            }
        }
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
    double *f=buffers[3].buf,*box=buffers[4].buf,*mass=buffers[5].buf;
    for(int d=0;d<3;++d) if (!(box[d]>2*cutoff) || !isfinite(box[d])) {
        PyErr_SetString(PyExc_ValueError, "Each box length must be greater than twice the cutoff.");goto cleanup;
    }
    for(int i=0;i<n;++i) {
        if (!(mass[i]>0) || !isfinite(mass[i])) {PyErr_SetString(PyExc_ValueError,"Masses must be positive.");goto cleanup;}
        for(int d=0;d<3;++d) if(!isfinite(x[3*i+d])||!isfinite(v[3*i+d])||!isfinite(u[3*i+d])) {
            PyErr_SetString(PyExc_ValueError,"State arrays must be finite.");goto cleanup;
        }
    }
    nb.reference = malloc((size_t)n*3*sizeof(double));
    if(!nb.reference) {PyErr_NoMemory();goto cleanup;}
    if(build_neighbors(&nb,n,x,u,box,cutoff+skin)<0 || forces(&nb,n,x,box,cutoff,model,pair,epsilon,sigma,parameters,exclude,length,k2,k3,k4,bond_parameters,f,&potential,&virial)<0) goto cleanup;
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
            else if(drift(n,x,u,v,box,0.5*dt)<0)goto cleanup;
            thermostat(n,v,mass,T,gamma,dt,&rng);
            if(model==2){constrain_velocities(n,v,x,box,mass);if(rigid_drift(n,x,u,v,box,mass,0.5*dt,length,bond_parameters,&constraint_impulse)<0)goto cleanup;}
            else if(drift(n,x,u,v,box,0.5*dt)<0)goto cleanup;
        } else if(model==2){if(rigid_drift(n,x,u,v,box,mass,dt,length,bond_parameters,&constraint_impulse)<0)goto cleanup;}
        else if(drift(n,x,u,v,box,dt)<0)goto cleanup;
        if(targetP>=0){
            double factor=exp(0.5*dt*rate);
            for(int d=0;d<3;++d){box[d]*=factor;if(!(box[d]>2*cutoff)){PyErr_SetString(PyExc_ValueError,"Box contracted below twice the cutoff.");goto cleanup;}}
            for(int k=0;k<3*n;++k){x[k]*=factor;u[k]*=factor;}
        }
        if(needs_rebuild(&nb,n,u,box,cutoff,skin) && build_neighbors(&nb,n,x,u,box,cutoff+skin)<0)goto cleanup;
        if(forces(&nb,n,x,box,cutoff,model,pair,epsilon,sigma,parameters,exclude,length,k2,k3,k4,bond_parameters,f,&potential,&virial)<0)goto cleanup;
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
    pressure_numerator=model?molecular_pressure_numerator(&nb,n,x,v,mass,box,cutoff,pair,epsilon,sigma,parameters):2*kinetic(n,v,mass)+virial;
    result=Py_BuildValue("(ddKiidd)",potential,virial,(unsigned long long)rng,done,interrupted,rate,pressure_numerator);

cleanup:
    if(bond_buffer.obj)PyBuffer_Release(&bond_buffer);
    if(parameter_buffer.obj)PyBuffer_Release(&parameter_buffer);
    free(nb.pairs);free(nb.reference);
    for(int k=0;k<6;++k) if(buffers[k].obj) PyBuffer_Release(&buffers[k]);
    return result;
}
static PyMethodDef methods[] = {
    {"advance", advance, METH_VARARGS, "Advance LJ dynamics in-place; return energy, virial, RNG, steps, interrupt flag, barostat rate and pressure numerator."},
    {NULL,NULL,0,NULL}
};
static struct PyModuleDef module = {PyModuleDef_HEAD_INIT,"_core",NULL,-1,methods};
PyMODINIT_FUNC PyInit__core(void) {return PyModule_Create(&module);}
