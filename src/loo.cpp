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

Loo* new_loo(int *indices, void *pflx, void *ptxs, void *pfxs, void *psxs,
             void *pcur, void *pqcur)
{
    Loo* loo = new Loo(indices, pflx, ptxs, pfxs, psxs, pcur, pqcur);
    return loo;
}

meshElement::meshElement(int ng, int nx, int ny, int nz, void *p)
    : _ng(ng), _nx(nx), _ny(ny), _nz(nz), _value((double*)p) { }

meshElement::~meshElement(){ }

double meshElement::getValue(int i, int j, int k, int g) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g < _ng);
    int index_f = g + i * _ng + j * _ng * _nx + k * _ng * _nx * _ny;
    return _value[index_f];
}

void Loo::printElement(meshElement element, std::string string){
    printf("%s ", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    printf("%f \n", element.getValue(g, i, j, k));
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

void Loo::printElement(surfaceElement element, std::string string){
    printf("%s \n", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                // FIXME: temporary, should be _ng
                for (int g = 0; g < 1; g++) {
                    // FIXME: temporary, should be _ns
                    for (int s = 0; s < element.getNs(); s++) {
                        printf("(%d %d %d) g = %d, s = %d: %f \n",
                               i, j, k, g, s, element.getValue(s, g, i, j, k));
                    }}}}}
}

/* element1 = quad_current, element2 = current*/
void Loo::verifyPartialCurrent(surfaceElement element1, surfaceElement element2){
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
}

/**
 * Constructor
 * @param indices parameters that contain #cells x,y,z and #energy groups
 */
Loo::Loo(int *indices, void *pflx, void *ptxs, void *pfxs, void *psxs,
         void *pcur, void *pqcur)
    : _nx(indices[0]),
      _ny(indices[1]),
      _nz(indices[2]),
      _ng(indices[3]),
      // FIXME: CMFD: ns = 12; LOO: ns = 16 in 2D, ns = 48 in 3D
      _ns_2d(16),
      _ns_3d(12),
      _num_loop(_nx),
      _num_track(4*_nx),
      _i_array(new int[_num_loop * _num_track]),
      _t_array(new int[_num_loop * _num_track]),
      _t_arrayb(new int[_num_loop *_num_track]),
      _old_flux(_ng, _nx, _ny, _nz, pflx),
      _total_xs(_ng, _nx, _ny, _nz, ptxs),
      _nfiss_xs(_ng, _nx, _ny, _nz, pfxs),
      _scatt_xs(_ng, _nx, _ny, _nz, psxs),
      _current(_ns_3d, _ng, _nx, _ny, _nz, pcur),
      _quad_current(_ns_2d, _ng, _nx, _ny, _nz, pqcur)
{
    generate2dTrack(_i_array, _t_array, _t_arrayb);

    printElement(_current, "current");
    printElement(_quad_current, "quad_current");
    //verifyPartialCurrent(_quad_current, _current);
}

/**
 *  Destructor: clear memory
 */
Loo::~Loo(){
}

void Loo::generate2dTrack(int *i_array, int *t_array, int *t_arrayb)
{
    int nl, nt, i, i_start;
    int nw = _nx;
    int len1 = _num_track / 2 - 1;
    int len2 = 1;
    for (nl = 0; nl < _num_loop; nl++)
    {
        i_start = (_ny - 1) * _nx + nl;

        for (nt = 0; nt < _num_track; nt++)
        {
            i = nl * _num_track + nt;

            /* First side: start with 0 index, goes like 0,5,0,5,... */
            if (nt < len1)
            {
                if (nt % 2 == 1)
                    i_start += 1;
                else if (nt != 0)
                    i_start -= nw;
                i_array[i] = i_start;


                if (nt % 2 == 0)
                {
                    t_array[i] = 0;
                    t_arrayb[i] = 1;
                }
                else
                {
                    t_array[i] = 5;
                    t_arrayb[i] = 4;
                }
            }
            /* 2nd side: start with odd index, goes like 2,7,... */
            else if (nt < len1 + len2)
            {
                if (nt % 2 == 0)
                    i_start -= nw;
                else if (nt != len1)
                    i_start -= 1;
                i_array[i] = i_start;

                if (nt % 2 == 1)
                {
                    t_array[i] = 2;
                    t_arrayb[i] = 3;
                }
                else
                {
                    t_array[i] = 7;
                    t_arrayb[i] = 6;
                }
            }
            /* 3rd side: start with even index, goes like 4,1,...*/
            else if (nt < len1 + len2 + len1)
            {
                if (nt % 2 == 1)
                    i_start -= 1;
                else if (nt != len1 + len2)
                    i_start += nw;
                i_array[i] = i_start;

                if (nt % 2 == 0)
                {
                    t_array[i] = 4;
                    t_arrayb[i] = 5;
                }
                else
                {
                    t_array[i] = 1;
                    t_arrayb[i] = 0;
                }
            }	
            /* last side */
            else
            {
                if (nt % 2 == 0)
                    i_start += nw;
                else if (nt != len1 + len2 + len1)
                    i_start += 1;
                i_array[i] = i_start;


                if (nt % 2 == 1)
                {
                    t_array[i] = 6;
                    t_arrayb[i] = 7;
                }
                else
                {
                    t_array[i] = 3;
                    t_arrayb[i] = 2;
                }
            }
        }
        len1 -= 2;
        len2 += 2;
    }

    return;
}

