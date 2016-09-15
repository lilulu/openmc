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
#include <time.h>       /* time_t, time, ctime */
#include <vector>
#include <fenv.h>       /* catch floating point exceptions */

#define FOUR_PI 12.566370614359172
#define TWO_PI 6.283185307179586
#define ONE_OVER_FOUR_PI 0.07957747154594767
#define PI 3.141592653589793
#define SMALL_THRESHOLD 1e-5
#define TINY_THRESHOLD 1e-15

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
    void normalize(int g, int i, int j, int k, double ratio);
    void printElement(std::string string, FILE* pfile);
    void zero();
    void one();
    void copyTo(meshElement& destination);
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
    void setValue(int g1, int g2, int i, int j, int k, double value);
    void printElement(std::string string, FILE* pfile);
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
    void normalize(int s, int g, int i, int j, int k, double value);
    void printElement(std::string string, FILE* pfile);
    void zero();
    void initialize(double);
};

class Loo {
private:
    bool _mapped;
    int _nx;
    int _ny;
    int _nz;
    int _ng;
    /* number of low-order tracks in a mesh cell */
    int _nt; 
    int _num_dimension;
    int _ns_2d;
    int _ns_3d;
    /* number of loops possible */
    int _num_loop;
    /* max number of trakcs in any loop; _num_tracks[loop_id] holds
     * the actual number of tracks in a given loop in the case where
     * coremap is specified */
    int _num_track;
    int _loo_iter;

    /* the mesh's x index in [0, _nx) given a global track id in [0,
     * _num_loop * _num_track) */
    int *_i_array;
    /* the mesh's y index in [0, _ny) given a global track id in [0,
     * _num_loop * _num_track) */    
    int *_j_array;
    /* the forward track id in [0, _nt) given a global track id in [0,
     * _num_loop * _num_track) */
    int *_t_array;
    /* the backward track id in [0, _nt) given a global track id in [0,
     * _num_loop * _num_track) */
    int *_t_arrayb;
    
    /* coremap passed in from main routine, defined in cmfd.xml */
    int *_map;
    int *_num_loops;
    int *_num_tracks;
    int *_boundary_cells;
    double _k;
    double _leakage;
    /* the global albedo passed in from main routine, set in cmfd.xml */
    double *_albedo;

    /* track lengthes in LOO: first calculated in 2D then projected
     * into 3D using 1 polar angle*/
    meshElement _track_length;

    /* 3D volume of mesh cells, stored in (0, i, j, k) */
    meshElement _volume;
    meshElement _old_scalar_flux;
    meshElement _scalar_flux;
    meshElement _total_xs;
    meshElement _abs_xs;

    /* m+1/2, not modified during LOO iterations */
    meshElement _sum_quad_flux;
    
    /* old fission source corresponding to the source before this
       batch of MC started) */
    meshElement _previous_fission_source;
    
    /* intend to pass back into Fortran: $\Qbar_{F,g}^{l}$ fission
     * source (including vol) generated by LOO */
    meshElement _fission_source;

    /* energy-integrated fission source */
    meshElement _energy_integrated_fission_source;

    /* energy dependent total source */
    meshElement _total_source;
    
    energyElement _nfiss_xs;
    energyElement _scatt_xs;
    meshElement _p1_scatt_xs;
    surfaceElement _length;
    surfaceElement _area;
    surfaceElement _quad_current;
    surfaceElement _quad_flux;
    surfaceElement _old_quad_flux;
    
    /* $\qhat^{m+1/2}$ divided by $\Qbar_S^{m+1/2} + \Qbar_F^{m}/k$ */
    surfaceElement _quad_src_form_factor;

    /* albedos (computed only at boundary cells), given track id,
     * energy group, mesh indexes */
    surfaceElement _albedos;
    
    FILE* _pfile;

public:
    Loo(int *indices, double *k, double *albedo,
        void *phxyz, void *pflx, void *ptso,
        void *ptxs, void *pfxs, void *psxs, void *pp1sxs,
        void *pqcur, void *pfs, int *pcor);
    virtual ~Loo();

    /* open log file for printing */
    void openLogFile();

    // getters, per say
    /* return true for active mesh */
    bool getMapValue(int, int);
    bool getMapped();
    double getK();    
    
    /* compute the surface areas and volume for each mesh cell */
    void computeAreaVolume();

    /* computes _track_length for each mesh cell based on mesh cell
     * length _length. Notice the generated length has already been
     * projected into one polar angle (using TY 1 polar angle
     * set). */
    void computeTrackLength();

    /* generate track laydown */
    void generate2dTrack();
    void generate2dTracknxn();
    void generate2dTracknxnMapped();
    void generate2dTracknx1();

    /* generate albedo conditions */
    void generateAlbedos();
    
    /* process _scalar_flux and _quad_current: the openmc generated
     * _scalar_flux and _quad_current are volume-integrated and
     * area-integrated respectively */
    void processFluxCurrent();

    /* compute absorption xs */
    void processXs();

    /* normalize tallies to remove the effect that
     * _previous_fission_source is accumulated from one less batch
     * compared with other tallies */
    void normalizeTallies();

    /* computes quad fluxes from quad currents, also save an old copy */
    void computeQuadFlux();

    /* computes the scattering quad source form factor associated with
     * each track and the sum of the eight quadrature fluxes in each
     * mesh cell */ 
    void computeQuadSourceFormFactor();

    /* iteratively solve the low-order problem using MOC (LOO) */
    void executeLoo();

    /* compute and store mesh cell _energy_integrated_fission_source:
     * nu_sigma_f * flux * vol, and return the average of the mesh
     * cell _energy_integrated_fission_source */
    double computeEnergyIntegratedFissionSource();

    /* compute loo source for outputing to external routines */
    void computeSource();

    /* compute quad_src */
    void computeQuadSource(surfaceElement& quad_src);

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

    bool cellOnBoundary(int, int, int, int);
    
    /* return bool representing whether a track starts from the mesh
     * geometry's specific surface, in the direction specified */
    bool startFromBoundary(int t, int i, int j, int k, int s, int dir);

    /* return bool representing whether a track starts from a vacuum
     * geometry boundary, in the direction specified */
    bool startFromAnyVacuumBoundary(int t, int i, int j, int k, int dir);

    bool startFromAnyBoundary(int t, int i, int j, int k, int dir);

    /* return area of the surface that a track t crosses with its
       start point (e = 0) or end point (e = 1) */
    double getSurfaceArea(int t, int i, int j, int k, int e);

    /* compute new mesh-cell averaged scalar flux _scalar_flux using
     * LOO1, LOO2, respectively */
    void computeScalarFlux(meshElement sum_quad_flux);
    void computeScalarFlux2(meshElement net_current);

    /* normalize fission source, scalar flux, quad current (if
     * initialization is set to be true) and leakage such that the
     * average of mesh-cell energy-integrated fission source is
     * old_avg */
    void normalizationByEnergyIntegratedFissionSourceAvg(double old_avg,
                                                         bool initialization);

    /* compute the L2 norm of relative change between the passed in
     * energy_integrated_fission_source variable (old fs from last
     * iteration) and the stored _energy_integrated_fission_source
     * (current fs from this iteration) */
    double computeL2Norm(meshElement fission_source);

    /* compute keff using updated scalar fluxes and leakage */
    void computeK();

    // helper routines
    void readInReferenceParameters();
    void checkBalance();
    /* debug: print out the coremap passed in from main routine */
    void printCoreMap();
};

extern "C" {
    double new_loo(int *indices, double *k, double *albedo,
                   void *phxyz, void *pflx, void *ptso,
                   void *ptxs, void *pfxs, void *psxs, void *pp1sxs,
                   void *pqcur, void *pfs, int *pcor);
}
#endif /* LOO_H_ */
