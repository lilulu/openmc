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
#define WT 4 * 3.141592653589793

double new_loo(int *indices, double *k, double *albedo,
               void *phxyz, void *pflx, void *ptso, void *ptxs, void *pfxs,
               void *psxs, void *pcur, void *pqcur, void *pfs)
{
    /* set up loo object */
    Loo loo = Loo(indices, k, albedo, phxyz, pflx, ptso,
                  ptxs, pfxs, psxs, pcur, pqcur, pfs);

    /* computes _quad_flux from _quad_current */
    loo.computeQuadFlux();

    /* normalizes _previous_fission_source to have energy-integrated
     * FS average to be 1.0. similarly for the current terms:
     * _scalar_flux, _quad_flux, _leakage, _energy_integrated_fission_source */
    loo.normalizeTallies();

    /* computes _quad_src from _quad_flux and total xs */
    loo.computeQuadSourceFormFactor();

    /* execute loo iterative solver */
    loo.executeLoo();

    return loo.getK();
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

void meshElement::copyTo(meshElement& destination) {
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    destination.setValue(g, i, j, k, getValue(g, i, j, k));
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
void meshElement::printElement(std::string string, FILE* pfile){
    fprintf(pfile, "%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    fprintf(pfile, "(%d %d %d) g = %d: %13.8e \n",
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

void energyElement::printElement(std::string string, FILE* pfile){
    fprintf(pfile, "%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g1 = 0; g1 < _ng; g1++) {
                    for (int g2 = 0; g2 < _ng; g2++) {
                        fprintf(pfile, "(%d %d %d) g1 = %d -> g2 = %d: %f \n",
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


void surfaceElement::printElement(std::string string, FILE* pfile){
    fprintf(pfile, "%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    for (int s = 0; s < _ns; s++) {
                        fprintf(pfile, "(%d %d %d) g = %d, s = %d: %f \n",
                                i, j, k, g, s, getValue(s, g, i, j, k));
                    }}}}}
    return;
}


/**
 * Constructor
 * @param indices parameters that contain #cells x,y,z and #energy groups
 */
Loo::Loo(int *indices, double *k, double* albedo,
         void *phxyz, void *pflx, void *ptso,
         void *ptxs, void *pfxs, void *psxs, void *pcur, void *pqcur, void *pfs)
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
      // FIXME:_num_loop, _num_track may not general for all mesh layouts
      _num_loop(_ny),
      _num_track(4 * _nx),
      _loo_iter(0),
      _i_array(new int[_num_loop * _num_track]),
      _j_array(new int[_num_loop * _num_track]),
      _t_array(new int[_num_loop * _num_track]),
      _t_arrayb(new int[_num_loop *_num_track]),
      /* double */
      _k(k[0]),
      _leakage(0.0),
      _albedo(albedo),
      _rms(0.0),
      /* meshElement */
      _track_length(1, _nx, _ny, _nz),
      _volume(1, _nx, _ny, _nz),
      _old_scalar_flux(_ng, _nx, _ny, _nz),
      _scalar_flux(_ng, _nx, _ny, _nz, pflx),
      _total_xs(_ng, _nx, _ny, _nz, ptxs),
      _abs_xs(_ng, _nx, _ny, _nz),
      _sum_quad_flux(_ng, _nx, _ny, _nz),
      _energy_integrated_fission_source(1, _nx, _ny, _nz),
      _previous_fission_source(_ng, _nx, _ny, _nz, ptso),
      _fission_source(_ng, _nx, _ny, _nz, pfs),
      /* energyElement */
      _nfiss_xs(_ng, _nx, _ny, _nz, pfxs),
      _scatt_xs(_ng, _nx, _ny, _nz, psxs),
      /* surfaceElement */
      _length(3, 1, _nx, _ny, _nz, phxyz),
      _area(3, 1, _nx, _ny, _nz),
      _current(_ns_3d, _ng, _nx, _ny, _nz, pcur),
      _quad_current(_ns_2d, _ng, _nx, _ny, _nz, pqcur),
      _quad_flux(_ns_2d, _ng, _nx, _ny, _nz),
      _old_quad_flux(_ns_2d, _ng, _nx, _ny, _nz),
      _quad_src_form_factor(_nt, _ng, _nx, _ny, _nz),
      _pfile(NULL)
{
    //printf("albedos: %d %d %d %d %d %d\n", _albedo[0], _albedo[1], _albedo[2],
    //     _albedo[3], _albedo[4], _albedo[5]);
    // FIXME: the following routine does not make any difference?
    _albedo[0] = 0.0;
    _albedo[1] = 0.0;
    _albedo[2] = 1.0;
    _albedo[3] = 1.0;
    _albedo[4] = 1.0;
    _albedo[5] = 1.0;
    openLogFile();
    computeAreaVolume();
    computeTrackLength();
    generate2dTrack();
    processFluxCurrent();
    processXs();
}

/**
 *  Destructor: clear memory
 *  FIXME: actually destroy stuff
 */
Loo::~Loo(){
}

/* open log file for printing */
void Loo::openLogFile() {
    time_t rawtime;
    char buffer [255];

    time (&rawtime);
    //sprintf(buffer,"log_%s.txt",ctime(&rawtime) );
    sprintf(buffer, "log.txt");

    /* Lets convert space to _ in */
    char *p = buffer;
    for (; *p; ++p) {
        if (*p == ' ')
            *p = '_';
    }

    /* open log file */
    _pfile = fopen(buffer, "a");

    return;
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
                l = 0.5 * sqrt(x * x + y * y);
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
    if (_nx == _ny) generate2dTracknxn();
    else if (_ny == 1) generate2dTracknx1();
    else printf("mesh of %d x %d might not be supported.\n", _nx, _ny);
    return;
}

/* this routine handles track generation for nx1 mesh cells. In this
 * case there is only one loop (regardless of whether the number of
 * cells is even or odd) */
void Loo::generate2dTracknx1() {
    /* nl = index of the loop (a loop is when tracks form a complete cycle) */
    /* nt = local index of track within a loop. */
    /* i = global index of track, taking into account nl and nt  */
    /* i_index, j_index = x, y index of cell that the current track is in */
    int nl, nt, i, i_index, j_index;

    /* len1 is the number of tracks from start point (left-most cell's
     * track 0) to when it first hits the right most cell's
     * boundary. This formula is the same for nxn and nx1 cases. */
    int len1 = _num_track / 2 - 1;

    for (nl = 0; nl < _num_loop; nl++) {
        /* when we start the nl-th loop, initialize the cell
         * indexes. We assume that the loop starts from the bottom
         * boundary in this plane. */
        i_index = nl;
        j_index = _ny - 1;

        for (nt = 0; nt < _num_track; nt++) {
            /* update global counter for track index */
            i = nl * _num_track + nt;

            /* First side: has len1 number of tracks, start with 0
             * index, goes like 0,5,3,6,... */
            if (nt < len1) {
                /* every odd track enters the cell to the right */
                if (nt % 2 == 1)
                    i_index += 1;

                /* we know the track index goes like 0,5,3,6 */
                if (nt % 4 == 0) {
                    _t_array[i] = 0;
                    _t_arrayb[i] = 1;
                }
                else if (nt % 4 == 1) {
                    _t_array[i] = 5;
                    _t_arrayb[i] = 4;
                }
                else if (nt % 4 == 2) {
                    _t_array[i] = 3;
                    _t_arrayb[i] = 2;
                }
                else {
                    _t_array[i] = 6;
                    _t_arrayb[i] = 7;
                }
            }
            /* 2nd side: has len1 + 1 number of tracks (two tracks per
             * cell), start with 2 or 1, goes like 2,4,1,7,2,4,... */
            else if (nt <= 2 * len1) {
                /* the first track has index len1 which is an odd
                 * number; except this first track, every odd track
                 * enters the cell to the left */
                if ((nt % 2 == 1) && (nt > len1))
                    i_index -= 1;

                if (nt % 4 == 0) {
                    _t_array[i] = 7;
                    _t_arrayb[i] = 6;
                }
                else if (nt % 4 == 1) {
                    _t_array[i] = 2;
                    _t_arrayb[i] = 3;
                }
                else if (nt % 4 == 2) {
                    _t_array[i] = 4;
                    _t_arrayb[i] = 5;
                }
                else {
                    _t_array[i] = 1;
                    _t_arrayb[i] = 0;
                }
            }
            /* last index: 6*/
            else {
                _t_array[i] = 6;
                _t_arrayb[i] = 7;
            }

            /* store the x, y cell indexes into the corresponding arrays */
            assert(i_index >= 0);
            assert(i_index < _nx);
            assert(j_index >= 0);
            assert(j_index < _ny);
            _i_array[i] = i_index;
            _j_array[i] = j_index;
        }
    }
    return;
}

void Loo::generate2dTracknxn() {
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
 * _scalar_flux and _current are volume-integrated and
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
                    scalar_flux = _scalar_flux.getValue(g, i, j, k)
                        / volume;
                    _old_scalar_flux.setValue(g, i, j, k, scalar_flux);
                    _scalar_flux.setValue(g, i, j, k, scalar_flux);

                    /* similarly, we need to divide the current by
                     * surface area. */
                    for (int s = 0; s < _quad_current.getNs(); s++) {
                        quad_current = _quad_current.getValue(s, g, i, j, k);

                        /* there are 8 quad current associated with
                         * each direction, and the index happens to
                         * work out that s/8 gives us the surface index
                         * corresponding to surface s */
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

/* normalize tallies to remove the effect that
 * _previous_fission_source is accumulated from one less batch
 * compared with other tallies */
void Loo::normalizeTallies() {
    double sum, factor;

    /* normalize the (m)-th order fission source to have average of
     * 1.0 (i.e., the sum should add up to be number of mesh
     * cells) */
    sum = _previous_fission_source.sum();
    factor = (double) _nx * _ny * _nz / sum;
    _previous_fission_source.normalize(factor);

    /* compute fission source generated by MC. */
    normalizationByEnergyIntegratedFissionSourceAvg(1.0);

    _previous_fission_source.printElement("m-th FS", _pfile);
    _energy_integrated_fission_source.printElement("m+1/2-th FS", _pfile);
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
                            /* FIXME: change sin 45 to reflect
                             * rectangular cell */
                            _quad_current.getValue(s, g, i, j, k) /
                            SIN_THETA_45);

                        /* store _quad_flux into _old_quad_flux. This
                         * is because we need a copy of the _quad_flux
                         * generated from MC */
                        _old_quad_flux.setValue(
                            s, g, i, j, k, _quad_flux.getValue(s, g, i, j, k));
                    }}}}}
    return;
}

/* compute the scattering quad source form factor associated with each
 * track and the sum of the eight quadrature fluxes in each mesh
 * cell */
void Loo::computeQuadSourceFormFactor(){
    double xs, l, ex, src, src_form_factor, sum_quad_flux, out, in, fs;
    double scattering_source;
    int in_index[] = {13, 5, 4, 11, 10, 2, 3, 12};
    int out_index[] = {6, 14, 8, 7, 1, 9, 15, 0};

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    sum_quad_flux = 0;
                    l = _track_length.getValue(0, i, j, k);
                    fs = _previous_fission_source.getValue(g, i, j, k);
                    xs = _total_xs.getValue(g, i, j, k);
                    /* Debug */
                    if (xs < 1e-5) {
                        printf("(%d %d %d) g = %d has tiny xs %f,"
                               "this would cause infinite src.\n",
                               i, j, k, g, xs);}
                    ex = exp(-xs * l);

                    /* computes m+1/2 scattering source */
                    scattering_source = 0;
                    for (int g2 = 0; g2 < _ng; g2++) {
                        scattering_source +=
                            _scatt_xs.getValue(g2, g, i, j, k)
                            * _scalar_flux.getValue(g2, i, j, k);
                    }
                    scattering_source *= _volume.getValue(0, i, j, k);

                    for (int t = 0; t < _nt; t++) {
                        in = _quad_flux.getValue(in_index[t], g, i, j, k);
                        out = _quad_flux.getValue(out_index[t], g, i, j, k);
                        src = xs * (out - ex * in) / (1.0 - ex) - fs / (WT * _k);
                        /* Debug: print out if negative quad src is *
                         generated. -1e-5 is used as the cutoff
                         because it looks like at least early on there
                         were multiple slightly negative values. */
                        if (src < -1e-5){
                            printf("A negative scattering quad src %e is "
                                   "generated for (%d %d %d) group %d "
                                   "track %d, out = %f, total src = %f, fs = %f\n",
                                   src, i, j, k, g, t, out,
                                   xs * (out - ex * in) / (1.0 - ex),
                                   fs / (WT * _k));
                        }
                        src_form_factor = src / scattering_source;
                        _quad_src_form_factor.setValue(t, g, i, j, k,
                                                       src_form_factor);

                        /* sum of quad flux uses lagged fission source */
                        sum_quad_flux += (src + fs / (WT * _k)) / xs
                            + (in - out) / (xs * l);
                    }
                    _sum_quad_flux.setValue(g, i, j, k, sum_quad_flux);
                }}}}

    _quad_src_form_factor.printElement("m+1/2 qs form factor", _pfile);
    _sum_quad_flux.printElement("m+1/2 sum quad flux", _pfile);
}

/* iteratively solve the low-order problem using MOC (LOO) */
void Loo::executeLoo(){
    fprintf(_pfile, "data passed into LOO from openmc:\n");
    //_quad_current.printElement("quad current", _pfile);
    _scalar_flux.printElement("scalar flux", _pfile);
    //_total_xs.printElement("total xs", _pfile);
    //_scatt_xs.printElement("scat xs", _pfile);
    //_nfiss_xs.printElement("fission xs", _pfile);

    /* loo iteration control */
    int min_loo_iter, max_loo_iter;
    double openmc_fs_avg, eps, eps_converged;

    /* memory allocation for data structure internal to this routine */
    meshElement fission_source (1, _nx, _ny, _nz);
    meshElement net_current (_ng, _nx, _ny, _nz);
    meshElement sum_quad_flux (_ng, _nx, _ny, _nz);
    surfaceElement quad_src (_nt, _ng, _nx, _ny, _nz);

    /* loop control variables: min, max number of loo sweeps to be performed */
    eps_converged = 1e-8;
    min_loo_iter = 10;
    max_loo_iter = 1000;

    /* save _fission_source into fission_source */
    _energy_integrated_fission_source.copyTo(fission_source);

    /* iteratively solve the LOO problem */
    for (_loo_iter = 0; _loo_iter < max_loo_iter; _loo_iter++) {

        /* reset net current, summmation of quad fluxes, leakage */
        net_current.zero();
        sum_quad_flux.zero();
        _leakage = 0;

        /* update quad_src using FIXME */
        computeQuadSource(quad_src);

        /* sweep through geometry, updating sum_quad_flux,
         * net_current, and _leakage */
        sweep(sum_quad_flux, net_current, quad_src);

        /* compute new mesh-cell averaged scalar flux _scalar_flux
         * using LOO1 */
        computeScalarFlux(sum_quad_flux, net_current);

        /* normalize scalar fluxes, quad fluxes, and leakage, by
         * calling computeEnergyIntegratedFissionSource */
        normalizationByEnergyIntegratedFissionSourceAvg(1.0);

        /* compute L2 norm of relative difference between successive
         * iterations, and compute k */
        eps = computeL2Norm(fission_source);
        computeK();

        fprintf(_pfile, "k = %.7f, eps = %f \n", _k, eps);

        /* save _energy_integrated_fission_source into fission_source */
        _energy_integrated_fission_source.copyTo(fission_source);

        /* check on convergence criteria */
        if ((eps < eps_converged) && (_loo_iter > min_loo_iter))
            break;
    }

    /* re-normalize _fission_source such that the average is 1.0, so
     * we can examine the solution generated by LOO */
    double normalize_factor = 1.0 / computeEnergyIntegratedFissionSource();
    _energy_integrated_fission_source.normalize(normalize_factor);
    _energy_integrated_fission_source.printElement(
        "converged energy integrated fission source", _pfile);
    fprintf(_pfile, "**********************************\n");
    fclose(_pfile);

    /* FIXME: cleans up memory */

    return;
}

/* compute and store mesh cell _energy_integrated_fission_source:
 * nu_sigma_f * flux * vol, compute and set _rms which is the computed
 * fission sources' deviation from a flat source distribution, and
 * return the average of the mesh cell
 * _energy_integrated_fission_source */
double Loo::computeEnergyIntegratedFissionSource() {
    double fission_source, sum, avg, rms;
    int counter;

    counter = 0;
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

                if (fission_source > 0) counter++;

                /* store the newly computed fission source into
                 * _fission_source */
                _energy_integrated_fission_source.setValue(0, i, j, k, fission_source);
            }}}

    /* compute average of the current fission source distribution */
    avg = _energy_integrated_fission_source.sum() / ((double) counter);

    /* compute rms relative deviation from flat fission source
     * distribution */
    sum = 0;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                sum += pow(_energy_integrated_fission_source.getValue
                           (0, i, j, k) / avg - 1, 2.0);
            }}}
    rms = sqrt(sum / ((double) counter));
    _rms = rms;

    /* finally return the normalization factor to make the current
     * fission source's average to come out to match old_avg */
    return avg;
}

/* compute quad_src for low-order iteration $l$ */
void Loo::computeQuadSource(surfaceElement& quad_src) {
    double scattering_source, fission_source, total_source;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g1 = 0; g1 < _ng; g1++) {
                    /* initialize source for this mesh this energy to be zero */
                    scattering_source = 0;
                    fission_source = 0;

                    for (int g2 = 0; g2 < _ng; g2++) {
                        scattering_source +=
                            _scatt_xs.getValue(g2, g1, i, j, k)
                            * _scalar_flux.getValue(g2, i, j, k);
                        fission_source +=
                            _nfiss_xs.getValue(g2, g1, i, j, k)
                            * _scalar_flux.getValue(g2, i, j, k);
                    }

                    scattering_source *= _volume.getValue(0, i, j, k);
                    fission_source *= _volume.getValue(0, i, j, k);

                    _fission_source.setValue(g1, i, j, k, fission_source);
                    for (int t = 0; t < _nt; t++) {
                        total_source = fission_source / (WT * _k)
                            + scattering_source *
                            _quad_src_form_factor.getValue(t, g1, i, j, k);
                        quad_src.setValue(t, g1, i, j, k, total_source);
                    }}}}}
    return;
}

/* the main sweeping routine, updating _quad_flux, sum_quad_flux,
 * net_current */
void Loo::sweep(meshElement& sum_quad_flux, meshElement& net_current,
                surfaceElement quad_src) {
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
                psi = sweepOneTrack(sum_quad_flux, net_current,
                                    quad_src, psi, g, nt, 0);
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
                psi = sweepOneTrack(sum_quad_flux, net_current,
                                    quad_src, psi, g, nt, 1);
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
double Loo::sweepOneTrack(meshElement& sum_quad_flux,
                          meshElement& net_current,
                          surfaceElement quad_src,
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
    src = quad_src.getValue(t, g, i, j, k);
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

    if (_nx == _ny) {
        /* dir = 0 means forward:  t = 6, 2, 4, 0;
         * dir = 1 means backward: t = 5, 1, 3, 7 */

        /* left geometry boundary in 2D */
        if ((s == 0) && (t == 6 - dir) && (i == 0))
            return true;

        /* right geometry boundary in 2D */
        if ((s == 1) && (t == 2 - dir) && (i == _nx - 1))
            return true;

        /* top geometry boundary in 2D */
        if ((s == 2) && (t == 4 - dir) && (j == 0))
            return true;

        /* bottom geometry boundary in 2D */
        if ((s == 3) && (t == dir * 7) && (j == _ny - 1))
            return true;
        return false;
    }

    if (_ny == 1) {
        /* dir = 0 means forward:  t = 6, 1, 4, 0;
         * dir = 1 means backward: t = 5, 2, 3, 7 */

        /* left geometry boundary in 2D */
        if ((s == 0) && (t == 6 - dir) && (i == 0))
            return true;

        /* right geometry boundary in 2D */
        if ((s == 1) && (t == dir + 1) && (i == _nx - 1))
            return true;

        // FIXME: top and bottom boundaries do not quite work for
        // non-reflective BCs
        /* top geometry boundary in 2D */
        if ((s == 2) && (t == 4 - dir) && (j == 0))
            return true;

        /* bottom geometry boundary in 2D */
        if ((s == 3) && (t == dir * 7) && (j == _ny - 1))
            return true;
        return false;
    }
}

/* return bool representing whether a track starts from a vacuum
 * geometry boundary, in the direction specified */
bool Loo::startFromAnyVacuumBoundary(int t, int i, int j, int k, int dir) {
    /* loop over the 2 * _num_dimension geometry boundary */
    for (int s = 0; s < 2 * _num_dimension; s++) {
        if ((_albedo[s] < 1e-8) && (startFromBoundary(t, i, j, k, s, dir))) {
            //printf("track %d direction %d in (%d %d) starts at vac\n",
            //       t, dir, i, j);
            return true;
        }
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

/* normalize fission source, scalar flux, quad flux and leakage such
 * that the average of mesh-cell energy-integrated fission source is
 * old_avg */
void Loo::normalizationByEnergyIntegratedFissionSourceAvg(double avg) {
    double ratio;

    /* compute normalization factor such that the energy-integrated
     * fission source average is avg */
    ratio = avg / computeEnergyIntegratedFissionSource();

    /* normalize fission source, scalar flux, quad flux, and leakage */
    _energy_integrated_fission_source.normalize(ratio);
    _scalar_flux.normalize(ratio);
    //_quad_flux.normalize(ratio);
    _leakage *= ratio;
}

/* compute the L2 norm of relative change between the passed in
 * fission_source variable (old fs from last iteration) and the stored
 * _fission_source (current fs from this iteration) */
double Loo::computeL2Norm(meshElement fission_source) {
    int counter;
    double eps, sum;
    counter = 0;
    sum = 0;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                if (fission_source.getValue(0, i, j, k) > 0) {
                    counter++;
                    sum += pow(_energy_integrated_fission_source.getValue
                               (0, i, j, k) /
                               fission_source.getValue(0, i, j, k) - 1
                               , 2);
                }}}}
    sum /= (double) (counter);
    eps = sqrt(sum);
    return eps;
}

/* compute keff using updated scalar fluxes and leakage */
void Loo::computeK(){
    double fission_total, absorption_total;

    fission_total = _energy_integrated_fission_source.sum();
    absorption_total = 0;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                        absorption_total += _abs_xs.getValue(g, i, j, k)
                            * _scalar_flux.getValue(g, i, j, k)
                            * _volume.getValue(0, i, j, k);
                }}}}

    _leakage *= SIN_THETA_45;
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

double Loo::getRms() {
    return _rms;
}

double Loo::getK() {
    return _k;
}
