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
    void incrementValue(int g, int i, int j, int k, double value);
    void normalize(double ratio);
    void printElement(std::string string);
    void zero();
    double sum();
};

class energyElement {
protected:
    int _ng, _nx, _ny, _nz;
    double *_value;

public:
    energyElement(int ng, int nx, int ny, int nz, void *p);
    virtual ~energyElement();
    double getValue(int g1, int g2, int i, int j, int k);
    void printElement(std::string string);
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
    void normalize(double value);
    void printElement(std::string string);
    void zero();
};

class Loo {
private:
    int _nx;
    int _ny;
    int _nz;
    int _ng;
    int _nt;
    int _num_dimension;
    int _ns_2d;
    int _ns_3d;
    int _num_loop;
    int _num_track;
    int *_i_array;
    int *_j_array;
    int *_t_array;
    int *_t_arrayb;
    double _k;
    double _leakage;
    double *_albedo;

    /* track lengthes in LOO: first calculated in 2D then projected
     * into 3D using 1 polar angle*/
    meshElement _track_length;

    /* 3D volume of mesh cells */
    meshElement _volume;
    meshElement _old_scalar_flux;
    meshElement _scalar_flux;
    meshElement _total_xs;
    meshElement _abs_xs;
    meshElement _sum_quad_flux;
    meshElement _fission_source;
    /* total source at the end of last batch of MC, corresponding to
       the source before this batch of MC started) */
    meshElement _old_total_source;
    meshElement _total_source;
    energyElement _nfiss_xs;
    energyElement _scatt_xs;
    surfaceElement _length;
    surfaceElement _area;
    surfaceElement _current;
    surfaceElement _quad_current;
    surfaceElement _quad_flux;
    surfaceElement _old_quad_flux;
    surfaceElement _quad_src;

public:
    Loo(int *indices, double *k, double *albedo,
        void *phxyz, void *pflx, void *ptso, void *ptsn,
        void *ptxs, void *pfxs, void *psxs, void *pcur, void *pqcur);
    virtual ~Loo();

    // main methods
    /* compute the surface areas and volume for each mesh cell */
    void computeAreaVolume();

    /* computes _track_length for each mesh cell based on mesh cell
     * length _length. Notice the generated length has already been
     * projected into one polar angle (using TY 1 polar angle
     * set). */
    void computeTrackLength();

    /* generate track laydown */
    void generate2dTrack();

    /* process _scalar_flux and _current: the openmc generated
     * _scalar_flux and c_current are volume-integrated and
     * area-integrated respectively */
    void processFluxCurrent();

    /* compute absorption xs */
    void processXs();

    /* computes quad fluxes from quad currents, also save an old copy */
    void computeQuadFlux();

    /* computes source term for every track from neutron balance using
     * quad fluxes and mesh-cell averaged cross-sections tallied
     * during MC. */
    void computeQuadSourceFromClosure();

    /* iteratively solve the low-order problem using MOC (LOO) */
    void executeLoo();

    /* compute mesh cell energy-integrated fission source, and
     * return normalization factor such that the average is 1.0 */
    double computeFissionSource(double old_avg);

    /* compute mesh cell energy-dependent total source (fission +
     * scattering) and update the source term passed in by
     * reference */
    void computeTotalSource(meshElement& source);

    /* compute quad_src using the current total_source and class member
     * _old_total_source and _quad_src */
    void computeQuadSource(surfaceElement& quad_src, meshElement& total_source);

    /* the main sweeping routine, updating _quad_flux, sum_quad_flux,
     * net_current */
    void sweep(meshElement& sum_quad_flux, meshElement& net_current,
               surfaceElement quad_src);

    /* the sweeping routine for sweeping through one track (nt) in one
     * direction (0 is forward, 1 is backward) for one energy group (g),
     * updating _quad_flux, sum_quad_flux, net_current, return updated
     * psi */
    double sweepOneTrack(meshElement& sum_quad_flux,
                         meshElement& net_current,
                         surfaceElement quad_src,
                         double psi, int g, int nt, int direction);

    /* return bool representing whether a track starts from the mesh
     * geometry's specific surface, in the direction specified */
    bool startFromBoundary(int t, int i, int j, int k, int s, int dir);

    /* return bool representing whether a track starts from a vacuum
     * geometry boundary, in the direction specified */
    bool startFromAnyVacuumBoundary(int t, int i, int j, int k, int dir);

    /* return area of the surface that a track t crosses with its
       start point (e = 0) or end point (e = 1) */
    double getSurfaceArea(int t, int i, int j, int k, int e);

    /* compute new mesh-cell averaged scalar flux _scalar_flux using LOO1 */
    void computeScalarFlux(meshElement sum_quad_flux, meshElement net_current);

    /* normalize scalar flux, quad flux and leakage */
    void normalization(double old_avg);

    double computeL2Norm(meshElement fission_source);

    void computeK();

    // helper routines
    void verifyPartialCurrent(surfaceElement element1, surfaceElement element2);

};

extern "C" {
    Loo* new_loo(int *indices, double *k, double *albedo,
                 void *phxyz, void *pflx, void *ptso, void *ptsn,
                 void *ptxs, void *pfxs, void *psxs, void *pcur, void *pqcur);
}
#endif /* LOO_H_ */
