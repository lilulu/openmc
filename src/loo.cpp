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
#include "loo.h"

Loo* new_loo(int *indices, double *k, void *phxyz, void *pflx, void *ptso,
             void *ptxs, void *pfxs, void *psxs, void *pcur, void *pqcur)
{
    /* set up loo object */
    Loo* loo = new Loo(indices, k, phxyz, pflx, ptso,
                       ptxs, pfxs, psxs, pcur, pqcur);

    /* computes _quad_flux from _quad_current */
    loo->computeQuadFlux();

    /* computes _quad_src from _quad_flux and total xs */
    loo->computeQuadSrc();

    /* execute loo iterative solver */
    loo->executeLoo();

    return loo;
}

meshElement::meshElement(int ng, int nx, int ny, int nz, void *p)
    : _ng(ng), _nx(nx), _ny(ny), _nz(nz), _value((double*)p) { }
meshElement::meshElement(int ng, int nx, int ny, int nz)
    : _ng(ng), _nx(nx), _ny(ny), _nz(nz) {
    _value = new double[_ng * _nx * _ny * _nz];}

meshElement::~meshElement(){ }

double meshElement::getValue(int g, int i, int j, int k) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g < _ng);
    int index_f = g + i * _ng + j * _ng * _nx + k * _ng * _nx * _ny;
    return _value[index_f];
}

void meshElement::setValue(int g, int i, int j, int k, double value){
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g < _ng);
    int index_f = g + i * _ng + j * _ng * _nx + k * _ng * _nx * _ny;
    _value[index_f] = value;
    return;
}

void meshElement::zero(){
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    setValue(g, i, j, k, 0.0);
                }}}}
}



/* applies: scattering & fission cross-sections */
energyElement::energyElement(int ng, int nx, int ny, int nz, void *p)
    : _ng(ng), _nx(nx), _ny(ny), _nz(nz), _value((double*)p) { }

energyElement::~energyElement(){ }

double energyElement::getValue(int g2, int g1, int i, int j, int k) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g1 < _ng);
    assert(g2 < _ng);
    int index_f = g2 + g1 * _ng + i * _ng * _ng + j * _ng *  _ng * _nx
        + k * _ng * _ng * _nx * _ny;
    return _value[index_f];
}

/* applies: currents & quad_currents */
surfaceElement::surfaceElement(int ns, int ng, int nx, int ny, int nz, void *p)
    : _ns(ns), _ng(ng), _nx(nx), _ny(ny), _nz(nz), _value((double*)p) {}

surfaceElement::surfaceElement(int ns, int ng, int nx, int ny, int nz)
    : _ns(ns), _ng(ng), _nx(nx), _ny(ny), _nz(nz) {
    _value = new double[_ns * _ng * _nx * _ny * _nz];}

surfaceElement::~surfaceElement(){ }

double surfaceElement::getNs(){
    return _ns;
}

double surfaceElement::getValue(int s, int g, int i, int j, int k) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g < _ng);
    assert(s < _ns);
    int index_f = s + g * _ns + i * _ns * _ng + j * _ns *  _ng * _nx
        + k * _ns * _ng * _nx * _ny;
    return _value[index_f];
}

void surfaceElement::setValue(int s, int g, int i, int j, int k, double value) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g < _ng);
    assert(s < _ns);
    int index_f = s + g * _ns + i * _ns * _ng + j * _ns *  _ng * _nx
        + k * _ns * _ng * _nx * _ny;
    _value[index_f] = value;
    return;
}

void surfaceElement::zero() {
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    for (int s = 0; s < _ns; s++) {
                        setValue(s, g, i, j, k, 0.0);
                    }}}}}
    return;
}

/**
 * Constructor
 * @param indices parameters that contain #cells x,y,z and #energy groups
 */
Loo::Loo(int *indices, double* k, void *phxyz, void *pflx, void *ptso,
         void *ptxs, void *pfxs, void *psxs, void *pcur, void *pqcur)
    : _nx(indices[0]),
      _ny(indices[1]),
      _nz(indices[2]),
      _ng(indices[3]),
      // FIXME: CMFD: ns = 12; LOO: ns = 16 in 2D, ns = 48 in 3D
      _nt(8),
      _ns_2d(16),
      _ns_3d(12),
      _num_loop(_nx),
      _num_track(4*_nx),
      _i_array(new int[_num_loop * _num_track]),
      _j_array(new int[_num_loop * _num_track]),
      _t_array(new int[_num_loop * _num_track]),
      _t_arrayb(new int[_num_loop *_num_track]),
      _k(k[0]),
      _track_length(1, _nx, _ny, _nz),
      _volume(1, _nx, _ny, _nz),
      _scalar_flux(_ng, _nx, _ny, _nz, pflx),
      _total_xs(_ng, _nx, _ny, _nz, ptxs),
      _sum_quad_flux(_ng, _nx, _ny, _nz),
      _fission_source(1, _nx, _ny, _nz),
      _old_total_source(_ng, _nx, _ny, _nz, ptso),
      _nfiss_xs(_ng, _nx, _ny, _nz, pfxs),
      _scatt_xs(_ng, _nx, _ny, _nz, psxs),
      _length(3, 1, _nx, _ny, _nz, phxyz),
      _current(_ns_3d, _ng, _nx, _ny, _nz, pcur),
      _quad_current(_ns_2d, _ng, _nx, _ny, _nz, pqcur),
      _quad_flux(_ns_2d, _ng, _nx, _ny, _nz),
      _old_quad_flux(_ns_2d, _ng, _nx, _ny, _nz),
      _quad_src(_nt, _ng, _nx, _ny, _nz)
{
    printf("k=%f\n", _k);
    computeTrackLengthVolume();
    generate2dTrack();
}

/**
 *  Destructor: clear memory
 */
Loo::~Loo(){
}

/* Computes _track_length and _volume based on _length */
void Loo::computeTrackLengthVolume() {
    double x, y, l;

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                x = _length.getValue(0, 0, i, j, k);
                y = _length.getValue(1, 0, i, j, k);
                _volume.setValue(0, i, j, k, x * y);
                l = 0.5 * sqrt(x * x + y * y) / P0;
                _track_length.setValue(0, i, j, k, l);
            }}}

    return;
}

void Loo::generate2dTrack() {
    /* nl = index of the loop (a loop is when tracks form a complete cycle) */
    /* nt = local index of track within a loop. */
    /* i = global index of track, taking into account nl and nt  */
    /* i_index, j_index = x, y index of cell that the current track is in */
    int nl, nt, i, i_index, j_index;

    /* len1, len2 are the number of tracks along the sides of the
     * current loop. len1 represents the tracks 45 degree from the
     * x-axis, and len2 represents the tracks 135 degree from the
     * x-axis. They are initialized below for the 1st loop, and will
     * be updated accordinly.*/
    int len1 = _num_track / 2 - 1;
    int len2 = 1;

    for (nl = 0; nl < _num_loop; nl++) {
        /* when we start the nl-th loop, initialize the cell
         * indexes. We assume that the loop starts from the bottom
         * boundary in this plane. */
        i_index = nl;
        j_index = _ny - 1;

        for (nt = 0; nt < _num_track; nt++) {
            /* update global counter for track index */
            i = nl * _num_track + nt;

            /* First side: start with 0 index, goes like 0,5,0,5,... */
            if (nt < len1) {
                if (nt % 2 == 1)
                    i_index += 1;
                else if (nt != 0)
                    j_index -= 1;

                if (nt % 2 == 0) {
                    _t_array[i] = 0;
                    _t_arrayb[i] = 1;
                }
                else {
                    _t_array[i] = 5;
                    _t_arrayb[i] = 4;
                }
            }
            /* 2nd side: start with odd index, goes like 2,7,... */
            else if (nt < len1 + len2) {
                if (nt % 2 == 0)
                    j_index -= 1;
                else if (nt != len1)
                    i_index -= 1;

                if (nt % 2 == 1) {
                    _t_array[i] = 2;
                    _t_arrayb[i] = 3;
                }
                else {
                    _t_array[i] = 7;
                    _t_arrayb[i] = 6;
                }
            }
            /* 3rd side: start with even index, goes like 4,1,...*/
            else if (nt < len1 + len2 + len1) {
                if (nt % 2 == 1)
                    i_index -= 1;
                else if (nt != len1 + len2)
                    j_index += 1;

                if (nt % 2 == 0) {
                    _t_array[i] = 4;
                    _t_arrayb[i] = 5;
                }
                else {
                    _t_array[i] = 1;
                    _t_arrayb[i] = 0;
                }
            }
            /* last side */
            else {
                if (nt % 2 == 0)
                    j_index += 1;
                else if (nt != len1 + len2 + len1)
                    i_index += 1;

                if (nt % 2 == 1) {
                    _t_array[i] = 6;
                    _t_arrayb[i] = 7;
                }
                else {
                    _t_array[i] = 3;
                    _t_arrayb[i] = 2;
                }
            }

            /* store the x, y cell indexes into the corresponding arrays */
            assert(i_index >= 0);
            assert(i_index < _nx);
            assert(j_index >= 0);
            assert(j_index < _ny);
            _i_array[i] = i_index;
            _j_array[i] = j_index;
        }

        /* update the side lengths for the next loop */
        len1 -= 2;
        len2 += 2;
    }
    return;
}

void Loo::computeQuadFlux(){
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    for (int s = 0; s < _quad_current.getNs(); s++) {
                        /* compute _quad_flux based on _quad_current */
                        _quad_flux.setValue(
                            s, g, i, j, k,
                            _quad_current.getValue(s, g, i, j, k) / P0 /
                            SIN_THETA_45);

                        /* store _quad_flux into _old_quad_flux. This
                         * is because we need a copy of the _quad_flux
                         * generated from MC */
                        _old_quad_flux.setValue(
                            s, g, i, j, k,
                            _quad_flux.getValue(s, g, i, j, k));
                    }}}}}
    //printElement(_quad_current, "quad_current");
    //printElement(_quad_flux, "quad_flux");
}

void Loo::computeQuadSrc(){
    double xs, l, ex, src, sum_quad_flux, out, in;

    int in_index[] = {13, 4, 10, 3, 12, 2, 11, 5};
    int out_index[] = {6, 8, 1, 15, 0, 9, 7, 14};

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    sum_quad_flux = 0;
                    xs = _total_xs.getValue(g, i, j, k);
                    l = _track_length.getValue(0, i, j, k);
                    ex = exp(-xs * l);

                    for (int t = 0; t < _nt; t++) {
                        in = _quad_flux.getValue(in_index[t], g, i, j, k);
                        out = _quad_flux.getValue(out_index[t], g, i, j, k);
                        src = xs * (out - ex * in) / (1.0 - ex);
                        /* Debug: print out if negative quad src is *
                         generated. -1e-5 is used as the cutoff
                         because it looks like at least early on there
                         were multiple slightly negative values. */
                        if (src < -1e-5){
                            printf("A negative quad src %e is generated for "
                                   "(%d %d %d) group %d track %d\n",
                                   src, i, j, k, g, t);
                        }
                        _quad_src.setValue(t, g, i, j, k, src);
                        sum_quad_flux += src / xs + (in - out) / (xs * l);
                    }
                    _sum_quad_flux.setValue(g, i, j, k, sum_quad_flux);
                }}}}
}

/* iteratively solve the low-order problem using MOC (LOO) */
void Loo::executeLoo(){
    int loo_iter, min_loo_iter, max_loo_iter;

     /* memory allocation for data structure internal to this routine */
    meshElement fission_source (_ng, _nx, _ny, _nz);
    meshElement mesh_src (_ng, _nx, _ny, _nz);
    meshElement net_current (_ng, _nx, _ny, _nz);
    meshElement sum_quad_flux (_ng, _nx, _ny, _nz);
    surfaceElement quad_src (_nt, _ng, _nx, _ny, _nz);

    /* loop control variables: min, max number of loo sweeps to be performed */
    min_loo_iter = 100;
    max_loo_iter = 1;

    /* initialization */
    mesh_src.zero();
    net_current.zero();
    sum_quad_flux.zero();
    quad_src.zero();

    /* compute total fission source generated by MC */
    computeFissionSource();

    /* iteratively solve the LOO problem */
    for (loo_iter = 0; loo_iter < max_loo_iter; loo_iter++) {
        /* reset net current and summmation of quad fluxes */
        net_current.zero();
        sum_quad_flux.zero();

        /* compute mesh_src: scattering and fission source for every
         * mesh every energy group */
        mesh_src.zero();
        computeMeshSource(mesh_src);

        
    }

    /* cleans up memory */
    return;
}

/* computes energy-integrated fission source: nu_sigma_f * flux * vol (?) */
void Loo::computeFissionSource() {
    double fission_source;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {

                /* initialize fission source for this cell to be zero */
                fission_source = 0;

                /* loop through energies, incrementing fission source counter */
                for (int g1 = 0; g1 < _ng; g1++) {
                    for (int g2 = 0; g2 < _ng; g2++) {
                        fission_source += _nfiss_xs.getValue(g2, g1, i, j, k)
                            * _scalar_flux.getValue(g2, i, j, k)
                            * _volume.getValue(0, i, j, k);
                    }}

                /* setter */
                _fission_source.setValue(0, i, j, k, fission_source);
            }}}
    return;
}

/* compute mesh cell energy-independent total source (fission +
 * scattering) and update the source term passed in by reference */
void Loo::computeMeshSource(meshElement& source) {
    double src;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {

                for (int g1 = 0; g1 < _ng; g1++) {
                    /* initialize source for this mesh this energy to be zero */
                    src = 0;

                    for (int g2 = 0; g2 < _ng; g2++) {
                        src += (_scatt_xs.getValue(g2, g1, i, j, k)
                                + _nfiss_xs.getValue(g2, g1, i, j, k) / _k)
                            * _scalar_flux.getValue(g2, i, j, k);
                    }
                    src *= _volume.getValue(0, i, j, k);
                    /* setter */
                    source.setValue(0, i, j, k, src);
                }}}}
    return;
}

/* return the surface length that track t is crossing with its start
 * point (e = 0) or end point (e = 1) in cell (i,j,k) */
double Loo::returnSurfaceLength(int i, int j, int k, int t, int e) {
    /* e represents start point (e = 0) or end point (e = 1) */
    assert(e > -1);
    assert(e < 2);

    /* t represents track id in [0, 7] */
    assert(t > -1);
    assert(t < _nt);

    /* index determines whether the cell side length in x, y or z
     * direction will be used */
    int index = 1;

    /* track 0, 2, 4, 6's entering points and track 1, 3, 5, 7's
     * exiting points are on surfaces along the x-direction */
    if (((e == 0) && (t % 2 == 0)) || ((e == 1) && (t % 2 == 1)))
        index = 0;

    /* _length is a surfaceElement of dimension 3 x 1 x nx x ny x nz*/
    return _length.getValue(index, 0, i, j, k);
}

void Loo::printElement(meshElement element, std::string string){
    printf("%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    printf("(%d %d %d) g = %d: %e \n",
                           i, j, k, g, element.getValue(g, i, j, k));
                        }}}}
    return;

}

void Loo::printElement(energyElement element, std::string string){
    printf("%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g1 = 0; g1 < _ng; g1++) {
                    for (int g2 = 0; g2 < _ng; g2++) {
                        printf("(%d %d %d) g1 = %d -> g2 = %d: %f \n",
                               i, j, k, g1, g2,
                               element.getValue(g1, g2, i, j, k));
                    }}}}}
    return;

}

void Loo::printElement(surfaceElement element, std::string string){
    printf("%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                // FIXME: temporary, should be _ng
                for (int g = 0; g < 1; g++) {
                    for (int s = 0; s < element.getNs(); s++) {
                        printf("(%d %d %d) g = %d, s = %d: %f \n",
                               i, j, k, g, s, element.getValue(s, g, i, j, k));
                    }}}}}
    return;
}

/* element1 = quad_current, element2 = current*/
void Loo::verifyPartialCurrent(surfaceElement element1,
                               surfaceElement element2){
    double delta;
    printf("delta \n");
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    // FIXME: temporary, should be _ns
                    for (int s = 0; s < 12; s++) {
                        delta = element1.getValue(s, g, i, j, k) -
                            element2.getValue(s, g, i, j, k);
                        //if (fabs(delta) > 0) {
                            printf("(%d %d %d) g = %d, s = %d: %f \n",
                                   i, j, k, g, s, delta);
                            //}
                    }}}}}
    return;
}
