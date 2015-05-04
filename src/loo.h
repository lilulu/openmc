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

class meshElement {
protected:
    int _nx, _ny, _nz, _ng, _ns;
    double *_value;

public:
    meshElement(int nx, int ny, int nz, int ng, int ns, void *p);
    virtual ~meshElement();
    double getValue(int i, int j, int k, int g);
    double getValue(int i, int j, int k, int g1, int g2);
};

class surfaceElement : public meshElement {
public:
    surfaceElement(int nx, int ny, int nz, int ng, int ns, void *p);
    double getValue(int i, int j, int k, int g, int s);
};


class Loo {
private:
    int _nx;
    int _ny;
    int _nz;
    int _ng;
    int _ns;
    int _num_loop;
    int _num_track;
    int *_i_array;
    int *_t_array;
    int *_t_arrayb;
    meshElement _old_flux;
    meshElement _total_xs;
    meshElement _nfiss_xs;
    meshElement _scatt_xs;
    surfaceElement _current;
    
public:
    Loo(int *indices, void *pflx, void *ptxs, void *pfxs, void *psxs, void *pcur);
    virtual ~Loo();

    void generate2dTrack(int *i_array, int *t_array, int *t_arrayb);

    // helper methods
    void printMeshElement(meshElement element, std::string string);

};

extern "C" {
    Loo* new_loo(int *indices, void *pflx, void *ptxs, void *pfxs, void *psxs,
        void *pcur);
}
#endif /* LOO_H_ */
