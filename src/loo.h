/*
* Copyright (C) 2015 Lulu Li
* <lululi@mit.edu>
*
* This file is part of OpenMC
*
* OpenMC is distributed under the MIT/X license. Permission is hereby
* granted, free of charge, to any person obtaining a copy of this
* software and associated documentation files (the “Software”), to
* deal in the Software without restriction, including without
* limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject
* to the following conditions:

* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
* BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
* ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/
#ifndef LOO_H_
#define LOO_H_

#include <cstdio>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* P0 is sin theta from Tabuchi-Yamamoto set for 1 polar angle */
#define P0 0.798184
#define SIN_THETA_45 0.70710678118654746
#define FOUR_PI 12.566370614359172
#define TWO_PI 6.283185307179586
#define ONE_OVER_FOUR_PI 0.07957747154594767
#define PI 3.141592653589793

class meshElement {
protected:
    int _ng, _nx, _ny, _nz;
    double *_value;

public:
    meshElement(int ng, int nx, int ny, int nz, void *p);
    meshElement(int ng, int nx, int ny, int nz);
    virtual ~meshElement();
    double getValue(int g, int i, int j, int k);
    void setValue(int g, int i, int j, int k, double value);
    void zero();
};

class energyElement {
protected:
    int _ng, _nx, _ny, _nz;
    double *_value;

public:
    energyElement(int ng, int nx, int ny, int nz, void *p);
    virtual ~energyElement();
    double getValue(int g1, int g2, int i, int j, int k);
};


class surfaceElement {
protected:
    int _ns, _ng, _nx, _ny, _nz;
    double *_value;

public:
    surfaceElement(int ns, int ng, int nx, int ny, int nz);
    surfaceElement(int ns, int ng, int nx, int ny, int nz, void *p);
    virtual ~surfaceElement();
    double getNs();
    double getValue(int s, int g, int i, int j, int k);
    void setValue(int s, int g, int i, int j, int k, double value);
    void zero();
};

class Loo {
private:
    int _nx;
    int _ny;
    int _nz;
    int _ng;
    int _nt;
    int _ns_2d;
    int _ns_3d;
    int _num_loop;
    int _num_track;
    int *_i_array;
    int *_j_array;
    int *_t_array;
    int *_t_arrayb;
    double _k;
    meshElement _track_length;
    meshElement _volume;
    meshElement _scalar_flux;
    meshElement _total_xs;
    meshElement _sum_quad_flux;
    meshElement _fission_source;
    energyElement _nfiss_xs;
    energyElement _scatt_xs;
    surfaceElement _length;
    surfaceElement _current;
    surfaceElement _quad_current;
    surfaceElement _quad_flux;
    surfaceElement _old_quad_flux;
    surfaceElement _quad_src;

public:
    Loo(int *indices, double *k, void *phxyz, void *pflx, void *ptxs,
        void *pfxs, void *psxs, void *pcur, void *pqcur);
    virtual ~Loo();

    // main methods
    /* computes _track_length for each mesh cell based on mesh cell
     * length _length. Notice the generated length has already been
     * projected into one polar angle (using TY 1 polar angle
     * set). Also computes volume of each mesh cell. */
    void computeTrackLengthVolume();

    /* generate track laydown */
    void generate2dTrack();

    /* computes quad fluxes from quad currents, also save an old copy */
    void computeQuadFlux();

    /* computes source term for every track from neutron balance using
     * quad fluxes and mesh-cell averaged cross-sections tallied
     * during MC. */
    void computeQuadSrc();

    /* iteratively solve the low-order problem using MOC (LOO) */
    void executeLoo();

    /* compute mesh cell energy-integrated fission source */
    void computeFissionSource();

    /* return surface length that this track is crossing with its
       start point (e = 0) or end point (e = 1) */
    double returnSurfaceLength(int i, int j, int k, int t, int e);

    // helper methods
    void printElement(meshElement element, std::string string);
    void printElement(energyElement element, std::string string);
    void printElement(surfaceElement element, std::string string);
    void verifyPartialCurrent(surfaceElement element1, surfaceElement element2);

};

extern "C" {
    Loo* new_loo(int *indices, double *k, void *phxyz, void *pflx, void *ptxs,
                 void *pfxs, void *psxs, void *pcur, void *pqcur);
}
#endif /* LOO_H_ */
