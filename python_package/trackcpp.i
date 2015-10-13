// TRACKCPP - Particle tracking code
// Copyright (C) 2015  LNLS Accelerator Physics Group
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

%module trackcpp

%{
#include <trackcpp/elements.h>
#include <trackcpp/kicktable.h>
#include <trackcpp/auxiliary.h>
#include <trackcpp/accelerator.h>
#include <trackcpp/pos.h>
#include <trackcpp/tracking.h>
#include <trackcpp/flat_file.h>
#include <trackcpp/optics.h>
#include "interface.h"
%}

%include "carrays.i"
%include "std_string.i"
%include "std_vector.i"
%include "stl.i"

namespace std {
    %template(CppStringVector) vector<string>;
    %template(CppDoubleVector) vector<double>;
    %template(CppElementVector) vector<Element>;
    %template(CppDoublePosVector) vector< Pos<double> >;
    %template(CppDoubleMatrix) vector< vector<double> >;
    %template(CppDoubleMatrixVector) vector< vector< vector<double> > >;
    %template(CppTwissVector) vector<Twiss>;
}

%inline %{
double c_array_get(double* v, int i) {
    return v[i];
}
void c_array_set(double* v, int i, double x) {
    v[i] = x;
}
double get_double_max() {
    return DBL_MAX;
}

%}


%include "elements.i"
%include "accelerator.i"
%include "optics.i"
%include "linalg.i"

%include "../include/trackcpp/kicktable.h"
%include "../include/trackcpp/auxiliary.h"
%include "../include/trackcpp/pos.h"
%include "../include/trackcpp/tracking.h"
//%include "../include/trackcpp/flat_file.h"
%include "interface.h"

%template(CppDoublePos) Pos<double>;
