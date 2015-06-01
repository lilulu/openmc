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

class meshElement {
protected:
    int _nx, _ny, _nz, _ng;
    double *_value;

public:
    meshElement(int ng, int nx, int ny, int nz, void *p);
    virtual ~meshElement();
    double getValue(int g, int i, int j, int k);
};

class energyElement {
protected:
    int _nx, _ny, _nz, _ng;
    double *_value;

public:
    energyElement(int ng, int nx, int ny, int nz, void *p);
    virtual ~energyElement();
    double getValue(int g1, int g2, int i, int j, int k);
};


class surfaceElement {
protected:
    int _nx, _ny, _nz, _ng, _ns;
    double *_value;

public:
    surfaceElement(int ns, int ng, int nx, int ny, int nz, void *p);
    virtual ~surfaceElement();
    double getNs();
    double getValue(int s, int g, int i, int j, int k);
    void verifyPartialCurrent(surfaceElement element1, surfaceElement element2);
};


class Loo {
private:
    int _nx;
    int _ny;
    int _nz;
    int _ng;
    int _ns_2d;
    int _ns_3d;
    int _num_loop;
    int _num_track;
    int *_i_array;
    int *_t_array;
    int *_t_arrayb;
    meshElement _old_flux;
    meshElement _total_xs;
    energyElement _nfiss_xs;
    energyElement _scatt_xs;
    surfaceElement _current;
    surfaceElement _quad_current;

public:
    Loo(int *indices, void *pflx, void *ptxs, void *pfxs, void *psxs,
        void *pcur, void *pqcur);
    virtual ~Loo();

    void generate2dTrack(int *i_array, int *t_array, int *t_arrayb);

    // helper methods
    void printElement(meshElement element, std::string string);
    void printElement(surfaceElement element, std::string string);
    void verifyPartialCurrent(surfaceElement element1, surfaceElement element2);

};

extern "C" {
    Loo* new_loo(int *indices, void *pflx, void *ptxs, void *pfxs, void *psxs,
                 void *pcur, void *pqcur);
}
#endif /* LOO_H_ */
