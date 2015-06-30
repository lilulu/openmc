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

Loo* new_loo(int *indices, double *k, double *albedo,
             void *phxyz, void *pflx, void *ptso,
             void *ptxs, void *pfxs, void *psxs, void *pcur, void *pqcur)
{
    /* set up loo object */
    Loo* loo = new Loo(indices, k, albedo, phxyz, pflx, ptso,
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

void meshElement::incrementValue(int g, int i, int j, int k, double value){
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g < _ng);
    int index_f = g + i * _ng + j * _ng * _nx + k * _ng * _nx * _ny;
    _value[index_f] += value;
    return;
}

void meshElement::zero(){
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    setValue(g, i, j, k, 0.0);
                }}}}
    return;
}

double meshElement::sum() {
    double sum = 0;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    sum += getValue(g, i, j, k);
                }}}}
    return sum;
}

void meshElement::normalize(double ratio) {
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    setValue(g, i, j, k, getValue(g, i, j, k) * ratio);
                }}}}
    return;
}

/* print out every values in this data structure */
void meshElement::printElement(std::string string){
    printf("%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    printf("(%d %d %d) g = %d: %e \n",
                           i, j, k, g, getValue(g, i, j, k));
                        }}}}
    return;

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

void energyElement::printElement(std::string string){
    printf("%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g1 = 0; g1 < _ng; g1++) {
                    for (int g2 = 0; g2 < _ng; g2++) {
                        printf("(%d %d %d) g1 = %d -> g2 = %d: %f \n",
                               i, j, k, g1, g2,
                               getValue(g1, g2, i, j, k));
                    }}}}}
    return;

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

void surfaceElement::normalize(double ratio) {
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    for (int s = 0; s < _ns; s++) {
                        setValue(s, g, i, j, k,
                                 getValue(s, g, i, j, k) * ratio);
                    }}}}}
    return;
}


void surfaceElement::printElement(std::string string){
    printf("%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    for (int s = 0; s < _ns; s++) {
                        printf("(%d %d %d) g = %d, s = %d: %f \n",
                               i, j, k, g, s, getValue(s, g, i, j, k));
                    }}}}}
    return;
}


/**
 * Constructor
 * @param indices parameters that contain #cells x,y,z and #energy groups
 */
Loo::Loo(int *indices, double* k, double* albedo,
         void *phxyz, void *pflx, void *ptso,
         void *ptxs, void *pfxs, void *psxs, void *pcur, void *pqcur)
    : _nx(indices[0]),
      _ny(indices[1]),
      _nz(indices[2]),
      _ng(indices[3]),
      _nt(8),
      // FIXME: CMFD: ns = 12; LOO: ns = 16 in 2D
      // FIEMX: for 3D, _num_dimension 2-> 3, _ns 16-> 48
      _num_dimension(2),
      _ns_2d(16),
      _ns_3d(12),
      _num_loop(_nx),
      _num_track(4*_nx),
      _i_array(new int[_num_loop * _num_track]),
      _j_array(new int[_num_loop * _num_track]),
      _t_array(new int[_num_loop * _num_track]),
      _t_arrayb(new int[_num_loop *_num_track]),
      _k(k[0]),
      _leakage(0.0),
      _albedo(albedo),
      _track_length(1, _nx, _ny, _nz),
      _volume(1, _nx, _ny, _nz),
      _old_scalar_flux(_ng, _nx, _ny, _nz, pflx),
      _scalar_flux(_ng, _nx, _ny, _nz),
      _total_xs(_ng, _nx, _ny, _nz, ptxs),
      _abs_xs(_ng, _nx, _ny, _nz),
      _sum_quad_flux(_ng, _nx, _ny, _nz),
      _fission_source(1, _nx, _ny, _nz),
      _old_total_source(_ng, _nx, _ny, _nz, ptso),
      _nfiss_xs(_ng, _nx, _ny, _nz, pfxs),
      _scatt_xs(_ng, _nx, _ny, _nz, psxs),
      _length(3, 1, _nx, _ny, _nz, phxyz),
      _area(3, 1, _nx, _ny, _nz),
      _current(_ns_3d, _ng, _nx, _ny, _nz, pcur),
      _quad_current(_ns_2d, _ng, _nx, _ny, _nz, pqcur),
      _quad_flux(_ns_2d, _ng, _nx, _ny, _nz),
      _old_quad_flux(_ns_2d, _ng, _nx, _ny, _nz),
      _quad_src(_nt, _ng, _nx, _ny, _nz)
{
    computeAreaVolume();
    computeTrackLength();
    generate2dTrack();
    processFluxCurrent();
    processXs();
}

/**
 *  Destructor: clear memory
 */
Loo::~Loo(){
}

/* compute _area and _volume using _length */
void Loo::computeAreaVolume() {
    double volume;

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                volume = _length.getValue(0, 0, i, j, k)
                * _length.getValue(1, 0, i, j, k)
                * _length.getValue(2, 0, i, j, k);
                _volume.setValue(0, i, j, k, volume);

                for (int s = 0; s < 3; s++) {
                    _area.setValue(s, 0, i, j, k,
                                   volume / _length.getValue(s, 0, i, j, k));
                }}}}
    return;
}

/* Computes _track_length and _volume using _length */
void Loo::computeTrackLength() {
    double x, y, z, l;

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                x = _length.getValue(0, 0, i, j, k);
                y = _length.getValue(1, 0, i, j, k);
                z = _length.getValue(2, 0, i, j, k);
                _volume.setValue(0, i, j, k, x * y * z);
                l = 0.5 * sqrt(x * x + y * y) / P0;
                _track_length.setValue(0, i, j, k, l);
            }}}

    return;
}


/* Track generation in 2D. In generating tracks, we use the following convention
 * for numbering the 8 tracks inside of each mesh cell.
 *           _
 *    //\   \\
 *  4//5    3\\2
 * \//       _\\
 *   _
 *  \\        //\
 *  6\\7    1//0
 *   _\\   \//
 *
 */
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

/* process _scalar_flux and _current: the openmc generated
 * _scalar_flux and c_current are volume-integrated and
 * area-integrated respectively */
void Loo::processFluxCurrent() {
    double scalar_flux, quad_current, volume, area;

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                volume = _volume.getValue(0, i, j, k);
                for (int g = 0; g < _ng; g++) {

                    /* the scalar flux passed in from openmc has
                       volume in it. This step divides it by volume so
                       _scalar_flux is the real scalar flux.  */
                    scalar_flux = _old_scalar_flux.getValue(g, i, j, k)
                        / volume;
                    _old_scalar_flux.setValue(g, i, j, k, scalar_flux);
                    _scalar_flux.setValue(g, i, j, k, scalar_flux);

                    /* similarly, we need to divide the current by
                     * surface area. */
                    for (int s = 0; s < _quad_current.getNs(); s++) {
                        quad_current = _quad_current.getValue(s, g, i, j, k);

                        /* there are 8 quad current associated with
                         * each direction, hence we find the
                         * corresponding area by s / 8 */
                        area = _area.getValue(s / 8, 0, i, j, k);

                        _quad_current.setValue(s, g, i, j, k, quad_current
                                               / area);
                    }}}}}
    return;
}

/* compute absorption xs */
void Loo::processXs() {
    double abs_xs, scatt_xs;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g1 = 0; g1 < _ng; g1++) {
                    scatt_xs = 0;
                    for (int g2 = 0; g2 < _ng; g2++) {
                        scatt_xs += _scatt_xs.getValue(g1, g2, i, j, k);
                    }
                    abs_xs = _total_xs.getValue(g1, i, j, k) - scatt_xs;
                    _abs_xs.setValue(g1, i, j, k, abs_xs);
                }}}}
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
                            s, g, i, j, k, _quad_flux.getValue(s, g, i, j, k));
                    }}}}}
    return;
}

/* compute the quadrature source associated with each track and the
 * sum of the eight quadrature fluxes in each mesh cell */
void Loo::computeQuadSrc(){
    double xs, l, ex, src, sum_quad_flux, out, in;

    int in_index[] = {13, 5, 4, 11, 10, 2, 3, 12};
    int out_index[] = {6, 14, 8, 7, 1, 9, 15, 0};

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
    /* loo iteration control */
    int loo_iter, min_loo_iter, max_loo_iter;

     /* memory allocation for data structure internal to this routine */
    meshElement total_source (_ng, _nx, _ny, _nz);
    meshElement net_current (_ng, _nx, _ny, _nz);
    meshElement sum_quad_flux (_ng, _nx, _ny, _nz);
    surfaceElement quad_src (_nt, _ng, _nx, _ny, _nz);

    /* loop control variables: min, max number of loo sweeps to be performed */
    min_loo_iter = 100;
    max_loo_iter = 1;

    /* initialization */
    total_source.zero();
    quad_src.zero();

    /* compute total fission source generated by MC */
    computeFissionSource();

    /* iteratively solve the LOO problem */
    for (loo_iter = 0; loo_iter < max_loo_iter; loo_iter++) {
        /* reset net current, summmation of quad fluxes, leakage */
        net_current.zero();
        sum_quad_flux.zero();
        _leakage = 0;

        /* compute total_source: scattering and fission source for
         * every mesh every energy group */
        total_source.zero();
        computeTotalSource(total_source);

        /* update quad_src using total_source, _old_total_source,
         * _quad_src */
        computeQuadSource(quad_src, total_source);

        /* sweep through geometry, updating sum_quad_flux,
         * net_current, and _leakage */
        sweep(sum_quad_flux, net_current);

        /* compute new mesh-cell averaged scalar flux _scalar_flux
         * using LOO1 */
        computeScalarFlux(sum_quad_flux, net_current);

        /* normalize scalar fluxes, quad fluxes, and leakage, by
         * calling computeFissionSource */
        normalization();

        computeK();
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
    _fission_source.printElement("fission source");
    return;
}

/* compute mesh cell energy-independent total source (fission +
 * scattering) and update the source term passed in by reference */
void Loo::computeTotalSource(meshElement& source) {
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

/* update quad_src using the current total_source and class member
 * _old_total_source and _quad_src */
void Loo::computeQuadSource(surfaceElement& quad_src,
                            meshElement& total_source) {
    double src_ratio;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    /* computes the form factor applied to all the
                     * tracks in this cell */
                    src_ratio = total_source.getValue(g, i, j, k) /
                        _old_total_source.getValue(g, i, j, k);

                    /* loop through each track updating */
                    for (int t = 0; t < _nt; t++) {
                        quad_src.setValue(t, g, i, j, k,
                                          _quad_src.getValue(t, g, i, j, k)
                                          * src_ratio);
                    }}}}}
    return;
}

/* the main sweeping routine, updating _quad_flux, sum_quad_flux,
 * net_current */
void Loo::sweep(meshElement& sum_quad_flux, meshElement& net_current) {
    /* current angular flux */
    double psi;

    for (int g = 0; g < _ng; g++) {

        /* forward */
        for (int nl = 0; nl < _num_loop; nl++) {
            /* get initial psi for this loop of tracks. We always
             * start our track from the bottom boundary. That is, the
             * first track is mesh cell (nl, _ny - 1, 0)'s quad flux
             * 13 */
            psi = _quad_flux.getValue(13, g, nl, _ny - 1, 0) * _albedo[3];

            /* sweeping through tracks in the forward order */
            for (int nt = _num_track * nl; nt < _num_track * (nl + 1); nt++) {
                psi = sweepOneTrack(sum_quad_flux, net_current, psi, g, nt, 0);
            }

            /* handle exiting psi: store psi (if reflective) or tally
             * leakage (if vacuum)*/
            _quad_flux.setValue(13, g, nl, _ny - 1, 0, psi * _albedo[3]);
            _leakage += psi * getSurfaceArea(0, nl, _ny - 1, 0, 0)
                * (1 - _albedo[3]);
        }

        /* backward */
        for (int nl = 0; nl < _num_loop; nl++) {
            /* similar to forward loop, except track 12 instead of 13 */
            psi = _quad_flux.getValue(12, g, nl, _ny - 1, 0) * _albedo[3];

            /* sweeping through tracks in the forward order */
            for (int nt = _num_track * (nl + 1) - 1;
                 nt > _num_track * nl - 1; nt--) {
                psi = sweepOneTrack(sum_quad_flux,net_current, psi, g, nt, 1);
            }

            /* handle exiting psi: store psi (if reflective) or tally
             * leakage (if vacuum)*/
            _quad_flux.setValue(12, g, nl, _ny - 1, 0, psi * _albedo[3]);
            _leakage += psi * getSurfaceArea(7, nl, _ny - 1, 0, 0)
                * (1 - _albedo[3]);
        }
    }

    return;
}

/* the sweeping routine for sweeping through one track (nt) in one
 * direction (0 is forward, 1 is backward) for one energy group (g),
 * updating _quad_flux, sum_quad_flux, net_current, return updated
 * psi */
double Loo::sweepOneTrack(meshElement& sum_quad_flux, meshElement& net_current,
                          double psi, int g, int nt, int direction) {
    int i, j, k, t;
    double delta, xs, l, src;

    i = _i_array[nt];
    j = _j_array[nt];

    /* for 2D problem, we only have k = 0 */
    k = 0;

    if (direction == 0)
        t = _t_array[nt];
    else
        t = _t_arrayb[nt];

    /* if we are on a vaccuum boundary, tally leakage and reset psi */
    if (startFromAnyVacuumBoundary(t, i, j, k, direction)) {
        _leakage += psi * getSurfaceArea(t, i, j, k, 0);
        psi = 0;
    }

    /* compute delta */
    xs = _total_xs.getValue(g, i, j, k);
    l = _track_length.getValue(0, i, j, k);
    src = _quad_src.getValue(t, g, i, j, k);
    delta = (psi - src / xs) * (1.0 - exp(-xs * l));

    /* update sum_quad_flux */
    sum_quad_flux.incrementValue(g, i, j, k, delta / (xs * l) + src / xs);

    /* update flux and net current */
    net_current.incrementValue(g, i, j, k,-psi * getSurfaceArea(t, i, j, k, 0));
    psi -= delta;
    net_current.incrementValue(g, i, j, k, psi * getSurfaceArea(t, i, j, k, 1));

    return psi;
}

/* return bool representing whether a track starts from the mesh
 * geometry's specific surface, in the direction specified */
bool Loo::startFromBoundary(int t, int i, int j, int k, int s, int dir) {

    /* dir = 0 means forward:  t = 6, 2, 4, 0;
     * dir = 1 means backward: t = 5, 1, 3, 7 */

    /* left geometry boundary in 2D */
    if ((s == 0) && (t == 6 - dir) && (j == 0))
        return true;

    /* right geometry boundary in 2D */
    if ((s == 1) && (t == 2 - dir) && (i == _nx - 1))
        return true;

    /* top surface in 2D */
    if ((s == 2) && (t == 4 - dir) && (i == 0))
        return true;

    /* bottom surface in 2D */
    if ((s == 3) && (t == dir * 7) && (j == _ny - 1))
        return true;

    return false;
}

/* return bool representing whether a track starts from a vacuum
 * geometry boundary, in the direction specified */
bool Loo::startFromAnyVacuumBoundary(int t, int i, int j, int k, int dir) {
    /* loop over the 2 * _num_dimension geometry boundary */
    for (int s = 0; s < 2 * _num_dimension; s++) {
        if ((_albedo[s] < 1e-10) && (startFromBoundary(t, i, j, k, s, dir)))
            return true;
    }
    return false;
}

/* return area of the surface that a track t crosses with its
   start point (e = 0) or end point (e = 1) */
double Loo::getSurfaceArea(int t, int i, int j, int k, int e) {
    /* e represents whether we are talking about track t's start point
     * (e = 0) or end point (e = 1) */
    assert(e > -1);
    assert(e < 2);

    /* t represents track id in [0, 7] */
    assert(t > -1);
    assert(t < _nt);

    /* index determines whether the cell side length in x, y or z
     * direction will be used */
    int index = 1;

    /* track 1, 2, 5, 6's start points and track 0, 3, 4, 7's
     * ending points are on surfaces perpendicular to the x-axis */
    if ((t == 1) || (t == 2) || (t == 5) || (t == 6)) {
        if (e == 0) index = 0;
    }
    else {
        if (e == 1) index = 0;
    }

    /* _area is a surfaceElement of dimension 3 x 1 x nx x ny x nz */
    return _area.getValue(index, 0, i, j, k);
}

/* compute new mesh-cell averaged scalar flux _scalar_flux using LOO1 */
void Loo::computeScalarFlux(meshElement sum_quad_flux, meshElement net_current){
    double phi_ratio, phi;

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    phi_ratio = sum_quad_flux.getValue(g, i, j, k) /
                        _sum_quad_flux.getValue(g, i, j, k);

                    phi = _old_scalar_flux.getValue(g, i, j, k) * phi_ratio;
                    _scalar_flux.setValue(g, i, j, k, phi);
                }}}}

    return;
}

/* normalize scalar flux, quad flux and leakage */
void Loo::normalization() {
    double old_sum, new_sum, ratio;

    old_sum = _fission_source.sum();
    computeFissionSource();
    new_sum = _fission_source.sum();
    ratio = old_sum / new_sum;

    _scalar_flux.normalize(ratio);
    _quad_flux.normalize(ratio);
    _leakage *= ratio;
}

void Loo::computeK(){
    double fission_total, absorption_total;

    fission_total = _fission_source.sum();
    absorption_total = 0;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                        absorption_total += _abs_xs.getValue(g, i, j, k)
                            * _scalar_flux.getValue(g, i, j, k)
                            * _volume.getValue(0, i, j, k);
                }}}}

    _k = fission_total / (absorption_total + _leakage);

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
