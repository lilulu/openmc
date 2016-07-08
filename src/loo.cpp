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
#define SIN_THETA_45 0.70710678118654746//1.0//0.9996937981813316//
#define P0 0.798184
#define WT_Q 8.0
#define REFERENCE 0          // control readInReferenceParameters
#define REFERENCE_TRANS_XS 0 // control whether to force right trans xs in 2G
#define DIVIDE 0             // divide by number of realizations in MC tallies

double new_loo(int *indices, double *k, double *albedo,
               void *phxyz, void *pflx, void *ptso, void *ptxs, void *pfxs,
               void *psxs, void *pp1sxs, void *pcur, void *pqcur, void *pfs)
{
    /* set up loo object */
    Loo loo = Loo(indices, k, albedo, phxyz, pflx, ptso,
                  ptxs, pfxs, psxs, pp1sxs, pcur, pqcur, pfs);

    /* divides scalar fluxes by volume, and currents by area */
    loo.processFluxCurrent();

    /* normalizes all the tallied terms such that the
     * energy-integrated FS average to be 1.0 */
    // DEBUG: when divided by number of realizations, do not need to
    // call normalizeTallies() but should normalize old FS when it's lagged MC
    if (DIVIDE) loo.normalizationByEnergyIntegratedFissionSourceAvg(1.0, true);
    else loo.normalizeTallies();

    /* computes _abs_xs from _total_xs and _scatt_xs */
    loo.processXs();

    /* Debug option: we can run the reference solution by reading in
     * parameters */
    //if (REFERENCE) loo.readInReferenceParameters();
    if (false) loo.readInReferenceParameters();

    /* computes _quad_flux from _quad_current */
    loo.computeQuadFlux();

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

void meshElement::one(){
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    setValue(g, i, j, k, 1.0);
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
                    if (getValue(g, i, j, k) > 1e-15) {
                        setValue(g, i, j, k, getValue(g, i, j, k) * ratio);
                    }}}}}
    return;
}

void meshElement::normalize(int g, int i, int j, int k, double ratio) {
    setValue(g, i, j, k, getValue(g, i, j, k) * ratio);
    return;
}

/* print out every values in this data structure */
void meshElement::printElement(std::string string, FILE* pfile){
    fprintf(pfile, "%s, printed k, j, i, g \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    fprintf(pfile, "%.8f, ", getValue(g, i, j, k));
                }}}}
    
    fprintf(pfile, "\n");
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

void energyElement::setValue(int g2, int g1, int i, int j, int k, double value) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g1 < _ng);
    assert(g2 < _ng);
    int index_f = g2 + g1 * _ng + i * _ng * _ng + j * _ng *  _ng * _nx
        + k * _ng * _ng * _nx * _ny;
    _value[index_f] = value;
    return;
}

void energyElement::printElement(std::string string, FILE* pfile){
    fprintf(pfile, "%s, printed in order of k, j, k, g1, g2\n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g1 = 0; g1 < _ng; g1++) {
                    for (int g2 = 0; g2 < _ng; g2++) {
                        fprintf(pfile, "%f, ", getValue(g1, g2, i, j, k));
                    }}}}}
    fprintf(pfile, "\n");

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

void surfaceElement::normalize(double ratio) {
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    for (int s = 0; s < _ns; s++) {
                        if (getValue(s, g, i, j, k) > 1e-15) {
                            setValue(s, g, i, j, k, getValue(s, g, i, j, k) * ratio);
                        }}}}}}
    return;
}

void surfaceElement::normalize(int s, int g, int i, int j, int k, double ratio) {
    setValue(s, g, i, j, k, getValue(s, g, i, j, k) * ratio);
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


void surfaceElement::one(){
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    for (int s = 0; s < _ns; s++) {
                        setValue(s, g, i, j, k, 1.0);
                    }}}}}
    return;
}

void surfaceElement::printElement(std::string string, FILE* pfile){
    fprintf(pfile, "%s in the order k, j, i, g, s\n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    for (int s = 0; s < _ns; s++) {
                        fprintf(pfile, "%f, ", getValue(s, g, i, j, k));
                    }
                    fprintf(pfile, "\n");
                }}}}
    return;
}

/**
 * Constructor
 * @param indices parameters that contain #cells x,y,z and #energy groups
 */
Loo::Loo(int *indices, double *k, double* albedo,
         void *phxyz, void *pflx, void *ptso,
         void *ptxs, void *pfxs, void *psxs, void *pp1sxs,
         void *pcur, void *pqcur, void *pfs)
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
      _previous_fission_source(_ng, _nx, _ny, _nz, ptso),
      _fission_source(_ng, _nx, _ny, _nz, pfs),
      _energy_integrated_fission_source(1, _nx, _ny, _nz),
      _total_source(_ng, _nx, _ny, _nz),
      /* energyElement */
      _nfiss_xs(_ng, _nx, _ny, _nz, pfxs),
      _scatt_xs(_ng, _nx, _ny, _nz, psxs),
      /* meshElement */
      _p1_scatt_xs(_ng, _nx, _ny, _nz, pp1sxs),
      _d(_ng, _nx, _ny, _nz),
      /* surfaceElement */
      _length(3, 1, _nx, _ny, _nz, phxyz),
      _area(3, 1, _nx, _ny, _nz),
      _current(_ns_3d, _ng, _nx, _ny, _nz, pcur),
      _quad_current(_ns_2d, _ng, _nx, _ny, _nz, pqcur),
      _quad_flux(_ns_2d, _ng, _nx, _ny, _nz),
      _old_quad_flux(_ns_2d, _ng, _nx, _ny, _nz),
      _quad_src_form_factor(_nt, _ng, _nx, _ny, _nz),
      _quad_src_total(_nt, _ng, _nx, _ny, _nz),
      _pfile(NULL)
{
    // FIXME: need to make sure the albedoes computed from cmfd are physical
    //printf("albedos: %f %f %f %f %f %f\n", _albedo[0], _albedo[1], _albedo[2],
    //     _albedo[3], _albedo[4], _albedo[5]);
    openLogFile();
    computeAreaVolume();
    computeTrackLength();
    generate2dTrack();
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

// Current parameters are from 1600 batches (200 inactive, 1400
// actives), 1 million neutrons/batch.
void Loo::readInReferenceParameters() {
    //_k = 1.339184; //1.3370953;

    double scalar_flux[] = 
        {0.00394678, 0.00103536, 0.00885761, 0.00246136, 0.01344737, 0.00374248, 0.01769299, 0.00492759, 0.02150122, 
         0.00598887, 0.02481962, 0.00691126, 0.02755156, 0.00767448, 0.02969247, 0.00826951, 0.03115085, 0.00867422, 
         0.03187850, 0.00887885, 0.03188963, 0.00887595, 0.03118538, 0.00868684, 0.02978393, 0.00829210, 0.02765423, 
         0.00769955, 0.02491412, 0.00694041, 0.02160858, 0.00601639, 0.01778203, 0.00495250, 0.01349390, 0.00375730, 
         0.00887795, 0.00246688, 0.00394774, 0.00103582};
    double previous_fission_source[] =
        {0.00017678, 0.00000000, 0.00041817, 0.00000000, 0.00063634, 
         0.00000000, 0.00083854, 0.00000000, 0.00102117, 0.00000000, 
         0.00118017, 0.00000000, 0.00131203, 0.00000000, 0.00141390, 
         0.00000000, 0.00148322, 0.00000000, 0.00151828, 0.00000000, 
         0.00151848, 0.00000000, 0.00148357, 0.00000000, 0.00141456, 
         0.00000000, 0.00131285, 0.00000000, 0.00118084, 0.00000000, 
         0.00102169, 0.00000000, 0.00083897, 0.00000000, 0.00063657, 
         0.00000000, 0.00041834, 0.00000000, 0.00017681, 0.00000000};

    double quad_current[] = 
        {0.000334, 0.000334, 0.000000, 0.000000, 0.000982, 0.000981, 0.000637, 0.000637, 0.000605, 0.000378, 0.000605, 0.000378, 0.000604, 0.000379, 0.000604, 0.000379, 
         0.000033, 0.000033, 0.000000, 0.000000, 0.000234, 0.000234, 0.000214, 0.000214, 0.000137, 0.000122, 0.000137, 0.000122, 0.000137, 0.000122, 0.000137, 0.000122, 
         0.000982, 0.000981, 0.000637, 0.000637, 0.001564, 0.001565, 0.001234, 0.001234, 0.001216, 0.001000, 0.001216, 0.001000, 0.001216, 0.000999, 0.001216, 0.000999, 
         0.000234, 0.000234, 0.000214, 0.000214, 0.000398, 0.000398, 0.000381, 0.000380, 0.000315, 0.000303, 0.000315, 0.000303, 0.000313, 0.000302, 0.000313, 0.000302, 
         0.001564, 0.001565, 0.001234, 0.001234, 0.002108, 0.002107, 0.001801, 0.001797, 0.001785, 0.001578, 0.001785, 0.001578, 0.001785, 0.001579, 0.001785, 0.001579, 
         0.000398, 0.000398, 0.000381, 0.000380, 0.000553, 0.000552, 0.000536, 0.000535, 0.000474, 0.000462, 0.000474, 0.000462, 0.000473, 0.000464, 0.000473, 0.000464, 
         0.002108, 0.002107, 0.001801, 0.001797, 0.002593, 0.002595, 0.002318, 0.002317, 0.002305, 0.002121, 0.002305, 0.002121, 0.002309, 0.002117, 0.002309, 0.002117, 
         0.000553, 0.000552, 0.000536, 0.000535, 0.000693, 0.000693, 0.000678, 0.000677, 0.000621, 0.000610, 0.000621, 0.000610, 0.000619, 0.000609, 0.000619, 0.000609, 
         0.002593, 0.002595, 0.002318, 0.002317, 0.003024, 0.003026, 0.002787, 0.002785, 0.002771, 0.002604, 0.002771, 0.002604, 0.002770, 0.002606, 0.002770, 0.002606, 
         0.000693, 0.000693, 0.000678, 0.000677, 0.000816, 0.000816, 0.000805, 0.000802, 0.000754, 0.000744, 0.000754, 0.000744, 0.000753, 0.000743, 0.000753, 0.000743, 
         0.003024, 0.003026, 0.002787, 0.002785, 0.003381, 0.003384, 0.003182, 0.003184, 0.003176, 0.003034, 0.003176, 0.003034, 0.003175, 0.003033, 0.003175, 0.003033, 
         0.000816, 0.000816, 0.000805, 0.000802, 0.000918, 0.000921, 0.000909, 0.000910, 0.000869, 0.000861, 0.000869, 0.000861, 0.000867, 0.000859, 0.000867, 0.000859, 
         0.003381, 0.003384, 0.003182, 0.003184, 0.003663, 0.003663, 0.003513, 0.003512, 0.003503, 0.003389, 0.003503, 0.003389, 0.003502, 0.003393, 0.003502, 0.003393, 
         0.000918, 0.000921, 0.000909, 0.000910, 0.001003, 0.001005, 0.000996, 0.000995, 0.000965, 0.000958, 0.000965, 0.000958, 0.000962, 0.000957, 0.000962, 0.000957, 
         0.003663, 0.003663, 0.003513, 0.003512, 0.003865, 0.003866, 0.003762, 0.003763, 0.003754, 0.003676, 0.003754, 0.003676, 0.003751, 0.003675, 0.003751, 0.003675, 
         0.001003, 0.001005, 0.000996, 0.000995, 0.001068, 0.001066, 0.001062, 0.001060, 0.001036, 0.001031, 0.001036, 0.001031, 0.001038, 0.001031, 0.001038, 0.001031, 
         0.003865, 0.003866, 0.003762, 0.003763, 0.003980, 0.003979, 0.003927, 0.003926, 0.003920, 0.003864, 0.003920, 0.003864, 0.003922, 0.003872, 0.003922, 0.003872, 
         0.001068, 0.001066, 0.001062, 0.001060, 0.001102, 0.001103, 0.001100, 0.001099, 0.001082, 0.001080, 0.001082, 0.001080, 0.001089, 0.001082, 0.001089, 0.001082, 
         0.003980, 0.003979, 0.003927, 0.003926, 0.003996, 0.003995, 0.003993, 0.003997, 0.003995, 0.003976, 0.003995, 0.003976, 0.003997, 0.003977, 0.003997, 0.003977, 
         0.001102, 0.001103, 0.001100, 0.001099, 0.001114, 0.001114, 0.001113, 0.001115, 0.001109, 0.001108, 0.001109, 0.001108, 0.001110, 0.001108, 0.001110, 0.001108, 
         0.003996, 0.003995, 0.003993, 0.003997, 0.003929, 0.003932, 0.003981, 0.003980, 0.003982, 0.003997, 0.003982, 0.003997, 0.003978, 0.003994, 0.003978, 0.003994, 
         0.001114, 0.001114, 0.001113, 0.001115, 0.001100, 0.001101, 0.001102, 0.001104, 0.001109, 0.001110, 0.001109, 0.001110, 0.001109, 0.001109, 0.001109, 0.001109, 
         0.003929, 0.003932, 0.003981, 0.003980, 0.003770, 0.003767, 0.003874, 0.003871, 0.003879, 0.003926, 0.003879, 0.003926, 0.003876, 0.003923, 0.003876, 0.003923, 
         0.001100, 0.001101, 0.001102, 0.001104, 0.001060, 0.001061, 0.001066, 0.001066, 0.001084, 0.001086, 0.001084, 0.001086, 0.001084, 0.001086, 0.001084, 0.001086, 
         0.003770, 0.003767, 0.003874, 0.003871, 0.003525, 0.003525, 0.003676, 0.003681, 0.003681, 0.003767, 0.003681, 0.003767, 0.003683, 0.003767, 0.003683, 0.003767, 
         0.001060, 0.001061, 0.001066, 0.001066, 0.001000, 0.001001, 0.001008, 0.001009, 0.001033, 0.001039, 0.001033, 0.001039, 0.001036, 0.001042, 0.001036, 0.001042, 
         0.003525, 0.003525, 0.003676, 0.003681, 0.003196, 0.003199, 0.003393, 0.003396, 0.003404, 0.003516, 0.003404, 0.003516, 0.003403, 0.003513, 0.003403, 0.003513, 
         0.001000, 0.001001, 0.001008, 0.001009, 0.000913, 0.000912, 0.000925, 0.000923, 0.000959, 0.000964, 0.000959, 0.000964, 0.000957, 0.000965, 0.000957, 0.000965, 
         0.003196, 0.003199, 0.003393, 0.003396, 0.002796, 0.002796, 0.003038, 0.003034, 0.003049, 0.003184, 0.003049, 0.003184, 0.003048, 0.003186, 0.003048, 0.003186, 
         0.000913, 0.000912, 0.000925, 0.000923, 0.000807, 0.000808, 0.000818, 0.000822, 0.000863, 0.000871, 0.000863, 0.000871, 0.000863, 0.000871, 0.000863, 0.000871, 
         0.002796, 0.002796, 0.003038, 0.003034, 0.002333, 0.002332, 0.002611, 0.002607, 0.002623, 0.002784, 0.002623, 0.002784, 0.002621, 0.002783, 0.002621, 0.002783, 
         0.000807, 0.000808, 0.000818, 0.000822, 0.000681, 0.000681, 0.000696, 0.000696, 0.000747, 0.000756, 0.000747, 0.000756, 0.000748, 0.000757, 0.000748, 0.000757, 
         0.002333, 0.002332, 0.002611, 0.002607, 0.001806, 0.001807, 0.002115, 0.002115, 0.002131, 0.002317, 0.002131, 0.002317, 0.002129, 0.002318, 0.002129, 0.002318, 
         0.000681, 0.000681, 0.000696, 0.000696, 0.000540, 0.000538, 0.000555, 0.000556, 0.000614, 0.000625, 0.000614, 0.000625, 0.000615, 0.000624, 0.000615, 0.000624, 
         0.001806, 0.001807, 0.002115, 0.002115, 0.001238, 0.001237, 0.001569, 0.001568, 0.001583, 0.001787, 0.001583, 0.001787, 0.001583, 0.001788, 0.001583, 0.001788, 
         0.000540, 0.000538, 0.000555, 0.000556, 0.000381, 0.000381, 0.000400, 0.000399, 0.000463, 0.000475, 0.000463, 0.000475, 0.000463, 0.000475, 0.000463, 0.000475, 
         0.001238, 0.001237, 0.001569, 0.001568, 0.000638, 0.000639, 0.000983, 0.000983, 0.001002, 0.001218, 0.001002, 0.001218, 0.001001, 0.001218, 0.001001, 0.001218, 
         0.000381, 0.000381, 0.000400, 0.000399, 0.000214, 0.000214, 0.000234, 0.000234, 0.000303, 0.000316, 0.000303, 0.000316, 0.000302, 0.000313, 0.000302, 0.000313, 
         0.000638, 0.000639, 0.000983, 0.000983, 0.000000, 0.000000, 0.000335, 0.000336, 0.000378, 0.000604, 0.000378, 0.000604, 0.000377, 0.000605, 0.000377, 0.000605, 
         0.000214, 0.000214, 0.000234, 0.000234, 0.000000, 0.000000, 0.000033, 0.000032, 0.000122, 0.000137, 0.000122, 0.000137, 0.000123, 0.000138, 0.000123, 0.000138};


    int c1 = 0, c2 = 0;
    for (int k = 0; k < _nz; ++k) {
        for (int j = 0; j < _ny; ++j) {
            for (int i = 0; i < _nx; ++i) {
                for (int g = 0; g < _ng; ++g) {
                    //_previous_fission_source.setValue(g, i, j, k,
                    //                                  previous_fission_source[c1]);
                    _scalar_flux.setValue(g, i, j, k, scalar_flux[c1]);
                    _old_scalar_flux.setValue(g, i, j, k, scalar_flux[c1]);
                    c1++;
                    for (int s = 0; s < _ns_2d; ++s) {
                        _quad_current.setValue(s, g, i, j, k, quad_current[c2]);
                        c2++;
                    }
                }}}}
    return;
}



/* process _scalar_flux and _current: the openmc generated
 * _scalar_flux and _current are volume-integrated and
 * area-integrated respectively */
void Loo::processFluxCurrent() {
  double scalar_flux, current, quad_current, volume, area;

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
		    for (int s = 0; s < _current.getNs(); s++) {
		      current = _current.getValue(s, g, i, j, k);
		      area = _area.getValue(s / 4, 0, i, j, k);
		      _current.setValue(s, g, i, j, k, current / area);
		    }

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

/* compute absorption xs: $\Sigma_{a,g} = \Sigma_{t,g} - \Sum_{g'}
 * \Sigma_{s, g\to g'} */
void Loo::processXs() {
    double tot_xs, trans_xs, abs_xs, scatt_xs, delta;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g1 = 0; g1 < _ng; g1++) {
                    scatt_xs = 0;
                    for (int g2 = 0; g2 < _ng; g2++) {
                        scatt_xs += _scatt_xs.getValue(g1, g2, i, j, k);
                    }

                    tot_xs = _total_xs.getValue(g1, i, j, k);
                    delta = _p1_scatt_xs.getValue(g1, i, j, k);
                    trans_xs = tot_xs - delta;
                    /* DEBUG: when reference parameters requested,
                     * overwrite trans_xs and delta */
                    if (REFERENCE_TRANS_XS) {
                        if (_ng == 1) trans_xs = 0.28011;
                        else if (_ng == 2) {
                            if (g1 == 0) trans_xs = 0.23148;
                            else trans_xs = 1.17371;
                        }
                        delta = tot_xs - trans_xs;
                    }
                    _total_xs.setValue(g1, i, j, k, trans_xs);

                    // absorption depends on the original total xs
                    abs_xs = tot_xs - scatt_xs;
                    _abs_xs.setValue(g1, i, j, k, abs_xs);
                    
                    // adjust in-group scattering s/t abs and scatt
                    // add up to trans xs
                    _scatt_xs.setValue(g1, g1, i, j, k, 
                                       _scatt_xs.getValue(g1, g1, i, j, k) - delta);

                    _d.setValue(g1, i, j, k, 1.0/3.0/trans_xs);
                }}}}
    return;
}

/* normalize tallies to remove the effect that
 * _previous_fission_source is accumulated from one less batch
 * compared with other tallies */
// DEBUG: making sure _previous_fission_source and fission source are
// normalized in a consistent manner
void Loo::normalizeTallies() {
    double sum, factor;

    // Debug: k passed in should not be zero
    if (_k < 1e-3) _k = 1.0;

    /* normalize the (m)-th order fission source to have average of
     * 1.0 (i.e., the sum should add up to be number of mesh
     * cells) */
    /* debug: the very first time LOO is turned on, old FS could be
       inf depends on how MC handles it */
    if (isinf(_previous_fission_source.getValue(0, 0, 0, 0))) {
      fprintf(_pfile, "found inf m-th FS, reset to 1.0");
      _previous_fission_source.one();
    }
    else {
      sum = _previous_fission_source.sum();
      if ((sum < 1e-5) or (sum > 1e5)) {
          _previous_fission_source.one();
    }}
    factor = (double) ( _nx * _ny * _nz) / _previous_fission_source.sum();
    _previous_fission_source.normalize(factor);

    /* compute fission source generated by MC. */
    normalizationByEnergyIntegratedFissionSourceAvg(1.0, true);
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
                            SIN_THETA_45 / P0);

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
    double xs, l, ex, src, src_form_factor, sum_quad_flux, out, in, fs, quad_src_total;
    double scattering_source, fission_source;
    int in_index[] = {13, 5, 4, 11, 10, 2, 3, 12};
    int out_index[] = {6, 14, 8, 7, 1, 9, 15, 0};

    // debug
    computeEnergyIntegratedFissionSource();

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    sum_quad_flux = 0;
                    l = _track_length.getValue(0, i, j, k);
                    fs = _previous_fission_source.getValue(g, i, j, k)
                        / _volume.getValue(0, i, j, k);
                    xs = _total_xs.getValue(g, i, j, k);
                    /* Debug */
                    if (xs < 1e-5) {
                        printf("(%d %d %d) g = %d has tiny xs %f,"
                               "this would cause infinite src.\n",
                               i, j, k, g, xs);
                    }

                    ex = exp(-xs * l);

                    /* computes m+1/2 scattering source */
                    scattering_source = 0;
                    fission_source = 0;
                    for (int g2 = 0; g2 < _ng; g2++) {
                        scattering_source +=
                            _scatt_xs.getValue(g2, g, i, j, k)
                            * _scalar_flux.getValue(g2, i, j, k);
                        fission_source += _nfiss_xs.getValue(g2, g, i, j, k)
                            * _scalar_flux.getValue(g2, i, j, k);
                    }

                    src = scattering_source + fs / _k;

                    double total = 0;
                    for (int t = 0; t < _nt; t++) {
                        in = _quad_flux.getValue(in_index[t], g, i, j, k);
                        out = _quad_flux.getValue(out_index[t], g, i, j, k);
                        quad_src_total = xs * (out - ex * in) / (1.0 - ex);

                        src_form_factor = quad_src_total / src * WT_Q;
                        total +=  src_form_factor;

                        _quad_src_form_factor.setValue(t, g, i, j, k,
                                                       src_form_factor);

                        /* Print out if negative quad src is *
                         generated. -1e-5 is used as the cutoff
                         because it looks like at least early on there
                         were multiple slightly negative values. */
                        if (src < -1e-5){
                            printf("A negative scattering quad src %e is "
                                   "generated for (%d %d %d) group %d "
                                   "track %d, out = %f, total src = %f, fs = %f, k = %f\n",
                                   src, i, j, k, g, t, out, quad_src_total, fs, _k);
                        }

                        /* sum of quad flux uses lagged fission source */
                        sum_quad_flux += quad_src_total / xs + (in - out) / (xs * l);
                        _quad_src_total.setValue(t,g, i, j, k, quad_src_total);
                    }

                    _sum_quad_flux.setValue(g, i, j, k, sum_quad_flux);
                }
            }
        }
    }

    _quad_src_form_factor.printElement("m+1/2 qs form factor", _pfile);
    //_sum_quad_flux.printElement("m+1/2 sum quad flux", _pfile);
}

/* iteratively solve the low-order problem using MOC (LOO) */
void Loo::executeLoo(){
    fprintf(_pfile, "data passed into LOO from openmc:\n");
    fprintf(_pfile, " _k = %f\n", _k);
    _previous_fission_source.printElement(" m-th fission source", _pfile);
    _energy_integrated_fission_source.printElement("fs", _pfile);
    _scalar_flux.printElement(" scalar flux", _pfile);
    //_total_xs.printElement(" total xs", _pfile);
    //_scatt_xs.printElement(" scat xs", _pfile);
    //_abs_xs.printElement(" abs xs", _pfile);
    //_p1_scatt_xs.printElement(" p1 scat xs", _pfile);
    //_d.printElement(" D", _pfile);
    //_nfiss_xs.printElement(" fission xs", _pfile);
    //_quad_current.printElement(" quad current", _pfile);
    checkBalance();
    
    /* loo iteration control */
    int min_loo_iter, max_loo_iter;
    double eps, eps_converged;

    /* memory allocation for data structure internal to this routine */
    meshElement fission_source (1, _nx, _ny, _nz);
    meshElement net_current (_ng, _nx, _ny, _nz);
    meshElement sum_quad_flux (_ng, _nx, _ny, _nz);
    surfaceElement quad_src (_nt, _ng, _nx, _ny, _nz);

    /* loop control variables: min, max number of loo sweeps to be performed */
    eps_converged = 1e-9;
    min_loo_iter = 10;
    max_loo_iter = 10000;

    /* save an old copy of energy-integrated fission source into
     * fission_source, so that when we update the term we can compute
     * L2 norm of relative change */
    // debug: call computequadsource an additional time here for debug
    // purpose
    computeQuadSource(quad_src);
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
        if (false) 
            computeScalarFlux(sum_quad_flux);
        else
            computeScalarFlux2(net_current);
        
        /* normalize scalar fluxes, quad fluxes, and leakage, by
         * calling computeEnergyIntegratedFissionSource */
        normalizationByEnergyIntegratedFissionSourceAvg(1.0, false);

        /* compute L2 norm of relative difference between successive
         * iterations, and compute k */
        eps = computeL2Norm(fission_source);
        computeK();

        double rms = 0;
        for (int i = 0; i < _nx; ++i) {
            for (int j = 0; j < _ny; ++j) {
                for (int k = 0; k < _nz; ++k) {
                    rms += pow(_scalar_flux.getValue(0, i, j, k) /
                               _old_scalar_flux.getValue(0, i, j, k) - 1.0, 2.0);
                }
            }
        }
        rms /= (double)(_nx * _ny * _nz);
        rms = pow(rms, 0.5);
        //fprintf(_pfile, "iter %d: k = %.7f, eps = %e, rms (phi) = %e\n",
        //        _loo_iter, _k, eps, rms);

        /* save _energy_integrated_fission_source into fission_source */
        _energy_integrated_fission_source.copyTo(fission_source);
	//if (_loo_iter % 20 == 0) {
	//  _energy_integrated_fission_source.printElement("fs", _pfile);
	//}
        /* check on convergence criteria */
        if ((eps < eps_converged) && (_loo_iter > min_loo_iter))
            break;
    }
    
    /* re-normalize _fission_source such that the average is 1.0, so
     * we can examine the solution generated by LOO */
    normalizationByEnergyIntegratedFissionSourceAvg(1.0, false);
    _energy_integrated_fission_source.printElement("fs", _pfile);
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
                 * _energy_integrated_fission_source */
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
    double scattering_source, fission_source, total_source, quad_source;
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g1 = 0; g1 < _ng; g1++) {
                    /* initialize source for this mesh this energy to be zero */
                    /* notice the two sources are not volume-integrated */
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


                    /* purpsoe: _fission_source passes fission source
                     * times volume to external files */
                    _fission_source.setValue(g1, i, j, k, fission_source
                                             * _volume.getValue(0, i, j, k));
                    
                    /* purpose: _total_source stores the current total
                     * source without volume */
                    total_source = scattering_source + fission_source / _k;
                    _total_source.setValue(g1, i, j, k, total_source);
                        
                    for (int t = 0; t < _nt; t++) {
                        quad_source = _quad_src_form_factor.getValue(t, g1, i, j, k)
                            / WT_Q * total_source;

                        quad_src.setValue(t, g1, i, j, k, quad_source);
                    }
                }
            }
        }
    }
    return;
}

/* the main sweeping routine, updating _quad_flux, sum_quad_flux,
 * net_current */
void Loo::sweep(meshElement& sum_quad_flux, meshElement& net_current,
                surfaceElement quad_src) {
    /* current angular flux */
    double psi, initial_psi;

    for (int g = 0; g < _ng; g++) {

        /* forward */
        for (int nl = 0; nl < _num_loop; nl++) {
            /* get initial psi for this loop of tracks. We always
             * start our track from the bottom boundary. That is, the
             * first track is mesh cell (nl, _ny - 1, 0)'s quad flux
             * 13 */
            psi = _quad_flux.getValue(13, g, nl, _ny - 1, 0) * _albedo[3];
            // debug
            initial_psi = psi;
            
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

            // debug
            if (false) { //fabs(initial_psi - psi) > 1e-5) {
                printf(" forward loop %d boundary fluxes has not converged: "
                       "%f -> %f\n", nl, initial_psi, psi);
            }
        }
        
        /* backward */
        for (int nl = 0; nl < _num_loop; nl++) {
            /* similar to forward loop, except track 12 instead of 13 */
            psi = _quad_flux.getValue(12, g, nl, _ny - 1, 0) * _albedo[3];

            // debug
            initial_psi = psi;
            
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

            // debug
            if (false) {//initial_psi - psi > -1e-5) {
                printf(" backward loop %d, boundary fluxes has not converged: "
                       "%f -> %f\n", nl, initial_psi, psi);
            }
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
//FIXME: might not even need dir as part of the logic here.
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
        /* dir = 0 means forward:  t = 6, 1 or 2, 4 or 3, 0 or 7;
         * dir = 1 means backward: t = 5, 2 or 1, 3 or 4, 7 or 0*/

        /* left geometry boundary in 2D */
        if ((s == 0) && (t == 6 - dir) && (i == 0))
            return true;

        /* right geometry boundary in 2D */
        if ((s == 1) && ((t == 1) || (t == 2)) && (i == _nx - 1))
            return true;

        // FIXME: top and bottom boundaries do not quite work for
        // non-reflective BCs
        /* top geometry boundary in 2D */
        if ((s == 2) && ((t == 4) || (t == 3)) && (j == 0))
            return true;

        /* bottom geometry boundary in 2D */
        if ((s == 3) && ((t == 7) || (t == 0)) && (j == _ny - 1))
            return true;
        return false;
    }
    return false;
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
void Loo::computeScalarFlux(meshElement sum_quad_flux){
    double phi_ratio, phi;

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
	        for (int g = 0; g < _ng; g++) {
                    phi_ratio = sum_quad_flux.getValue(g, i, j, k) /
		      _sum_quad_flux.getValue(g, i, j, k);
                    phi = _old_scalar_flux.getValue(g, i, j, k) * phi_ratio;
                    
                    //fprintf(_pfile, "(%d %d %d) flux update ratio = %e\n",
                    //       i, j, k, phi / _scalar_flux.getValue(g,i,j,k) - 1.0);
                    _scalar_flux.setValue(g, i, j, k, phi);
                }}}}
    return;
}

/* compute new mesh-cell averaged scalar flux _scalar_flux using LOO2 */
void Loo::computeScalarFlux2(meshElement net_current){
    double phi, vol, netcurrent;

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                vol = _volume.getValue(0, i, j, k);
                
                for (int g = 0; g < _ng; g++) {
                    netcurrent = net_current.getValue(g, i, j, k)
                        * SIN_THETA_45 * P0;

                    phi = (_total_source.getValue(g, i, j, k)
                           - netcurrent / vol) /
                        _total_xs.getValue(g, i, j, k);
                    
                    //fprintf(_pfile, "(%d %d %d) flux update ratio = %e\n",
                    //        i, j, k, phi / _scalar_flux.getValue(g,i,j,k) - 1.0);
                    _scalar_flux.setValue(g, i, j, k, phi);
                }}}}
    return;
}

/* normalize fission source, scalar flux, quad currents (only during
 * initialization), and leakage such that the average of mesh-cell
 * energy-integrated fission source is old_avg */
void Loo::normalizationByEnergyIntegratedFissionSourceAvg(double avg,
                                                          bool initialization) {
    double ratio;

    /* compute normalization factor such that the energy-integrated
     * fission source average is avg */
    ratio = avg / computeEnergyIntegratedFissionSource();

    /* purpose: normalize fission source, scalar flux, quad flux, and
     * leakage */
    _energy_integrated_fission_source.normalize(ratio);
    _scalar_flux.normalize(ratio);

    /* when in the initialization phase, we only normalize quad
     * current because quad flux has not been computed yet; after
     * initialization, we normalize quad flux */
    if (initialization) {
        _old_scalar_flux.normalize(ratio);
        _quad_current.normalize(ratio);
        _current.normalize(ratio);
        if (DIVIDE) _previous_fission_source.normalize(ratio);
    }
    else {
        _quad_flux.normalize(ratio);
    }

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

    _leakage *= SIN_THETA_45 * P0;
    _k = fission_total / (absorption_total + _leakage);

    return;
}

void Loo::checkBalance(){
    double fission, absorption, leakage, residual;

    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                fission = _energy_integrated_fission_source.getValue(0, i, j, k);
                absorption = 0;
                leakage = 0;
                for (int g = 0; g < _ng; g++) {
                    absorption += _abs_xs.getValue(g, i, j, k)
                        * _scalar_flux.getValue(g, i, j, k)
                        * _volume.getValue(0, i, j, k);

                    // computation of leakage from
                    // _quad_current. Convention: _quad_current
                    // leaving the cell is count as positive: 0,1,
                    // 6,7, 8,9, 14,15
                    leakage += _area.getValue(0, 0, i, j, k) *
                        (_quad_current.getValue(0, g, i, j, k)
                         + _quad_current.getValue(1, g, i, j, k)
                         - _quad_current.getValue(2, g, i, j, k)
                         - _quad_current.getValue(3, g, i, j, k)
                         - _quad_current.getValue(4, g, i, j, k)
                         - _quad_current.getValue(5, g, i, j, k)
                         + _quad_current.getValue(6, g, i, j, k)
                         + _quad_current.getValue(7, g, i, j, k));
                    leakage += _area.getValue(1, 0, i, j, k) *
                        (_quad_current.getValue(8, g, i, j, k)
                         + _quad_current.getValue(9, g, i, j, k)
                         - _quad_current.getValue(10, g, i, j, k)
                         - _quad_current.getValue(11, g, i, j, k)
                         - _quad_current.getValue(12, g, i, j, k)
                         - _quad_current.getValue(13, g, i, j, k)
                         + _quad_current.getValue(14, g, i, j, k)
                         + _quad_current.getValue(15, g, i, j, k));

                    residual = leakage + absorption - fission / _k;

                    /* if residual is larger than a certain threshold,
                     * print to screen */
                    if (true){ //fabs(residual / absorption) > 1e-4) {
                        fprintf(_pfile, "residual in (%d %d %d):"
                               " %.2e = %.2e + %.2e - %.2e / %.2e\n", i, j, k,
                               residual, leakage, absorption, fission,  _k);
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

double Loo::getRms() {
    return _rms;
}

double Loo::getK() {
    return _k;
}
