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
             void *pcur)
{
    Loo* loo = new Loo(indices, pflx, ptxs, pfxs, psxs, pcur);
    return loo;
}

meshElement::meshElement(int nx, int ny, int nz, int ng, int ns, void *p)
    : _nx(nx),
      _ny(ny),
      _nz(nz),
      _ng(ng),
      _ns(ns),
      _value((double*)p)
{
}

meshElement::~meshElement(){
}

/* applies: _old_flux, _total_xs 
double meshElement::getValue(int i, int j, int k, int g) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g < _ng);
    int index_f = g + i * _ng + j * _ng * _nx + k * _ng * _nx * _ny;
    return _value[index_f];
}

/* applies: scattering & fission cross-sections */
double meshElement::getValue(int i, int j, int k, int g1, int g2) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g1 < _ng);
    assert(g2 < _ng);
    int index_f = g2 + g1 * _ng + i * _ng * _ng + j * _ng *  _ng * _nx
        + k * _ng * _ng * _nx * _ny;
    return _value[index_f];
}

surfaceElement::surfaceElement(int nx, int ny, int nz, int ng, int ns, void *p)
    : meshElement(nx, ny, nz, ng, ns, p) { }

double surfaceElement::getValue(int i, int j, int k, int g, int s) {
    assert(i < _nx);
    assert(j < _ny);
    assert(k < _nz);
    assert(g < _ng);
    assert(s < _ns);
    int index_f = s + g * _ns + i * _ns * _ng + j * _ns *  _ng * _nx
        + k * _ns * _ng * _nx * _ny;
    return _value[index_f];
}

/**
 * Constructor
 * @param indices parameters that contain #cells x,y,z and #energy groups
 */
Loo::Loo(int *indices, void *pflx, void *ptxs, void *pfxs, void *psxs, void *pcur)
    : _nx(indices[0]),
      _ny(indices[1]),
      _nz(indices[2]),
      _ng(indices[3]),
      _ns(12),
      _num_loop(_nx),
      _num_track(4*_nx),
      _i_array(new int[_num_loop * _num_track]),
      _t_array(new int[_num_loop * _num_track]),
      _t_arrayb(new int[_num_loop *_num_track]),
      _old_flux(_nx, _ny, _nz, _ng, _ns, pflx),
      _total_xs(_nx, _ny, _nz, _ng, _ns, ptxs),
      _nfiss_xs(_nx, _ny, _nz, _ng, _ns, pfxs),
      _scatt_xs(_nx, _ny, _nz, _ng, _ns, psxs),
      _current(_nx, _ny, _nz, _ng, _ns, pcur)
{
    generate2dTrack(_i_array, _t_array, _t_arrayb);
// variables untested:
    //printMeshElement(_total_xs, "total xs");
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

void Loo::printMeshElement(meshElement element, std::string string){
    printf("%s ", string.c_str());
    for (int k = 0; k < _nz; k++) {
        for (int j = 0; j < _ny; j++) {
            for (int i = 0; i < _nx; i++) {
                for (int g = 0; g < _ng; g++) {
                    printf("%f ", element.getValue(i,j,k,g,g));
                        }}}}
    printf("\n");
}


/** 
 *  Destructor: clear memory
 */
Loo::~Loo(){
}
