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
#define SWIG_FILE_WITH_INIT
#include <trackcpp/elements.h>
#include <trackcpp/kicktable.h>
#include <trackcpp/auxiliary.h>
#include <trackcpp/accelerator.h>
#include <trackcpp/pos.h>
#include <trackcpp/tracking.h>
#include <trackcpp/flat_file.h>
#include <trackcpp/optics.h>
#include <trackcpp/diffusion_matrix.h>
#include <trackcpp/naff.h>
#include <trackcpp/linalg.h>
#include "interface.h"
%}

%include "carrays.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_list.i"
%include "stl.i"
%include "typemaps.i"
%include "std_except.i"


// Define a custom exception in Python
%pythoncode {
    class TrackingError(Exception):
        pass
}

%exception {
    try
    {
        $action
    }
    catch (const std::out_of_range &e)
    {
        PyObject *ex = PyObject_GetAttrString(
            PyImport_AddModule("trackcpp"), "TrackingError"
        );
        if (ex != NULL) PyErr_SetString(ex, e.what());
        else PyErr_SetString(PyExc_IndexError, e.what());
        SWIG_fail;
    }
    catch (...)
    {
        PyErr_SetString(PyExc_RuntimeError, "Unknown exception occurred");
        SWIG_fail;
    }
}


namespace std {
    %template(CppStringVector) vector<string>;
    %template(CppDoubleVector) vector<double>;
    %template(CppUnsigIntVector) vector<unsigned int>;
    %template(CppIntVector) vector<int>;
    %template(CppBoolVector) vector<bool>;
    %template(CppElementVector) vector<Element>;
    %template(CppDoublePosVector) vector< Pos<double> >;
    %template(CppDoubleMatrix) vector< vector<double> >;
    %template(CppDoubleMatrixVector) vector< vector< vector<double> > >;
    %template(CppTwissVector) vector<Twiss>;
    %template(CppMatrixVector) vector< Matrix >;
    %template(CppKicktableVector) list< Kicktable >;
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

%extend std::vector<double>
{
    double* data_()
    {
        return $self->data();
    }
}

%include "numpy.i"

%init %{
    import_array();
%}

// For ringpass and linepass
%apply (double* IN_ARRAY2, int DIM1, int DIM2 ){
    (double* orig_pos, int ni1, int ni2)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2 ) {
    (double* pos, int n1, int n2)}

// For calc_twiss
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2 ) {
    (double* twiss, int n1, int n2)}

//For find_m66 and diffusion_matrix
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 )
{
    (double *cumul_tm, int n1_tm, int n2_tm, int n3_tm)
}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 )
{
    (double *elem_tm, int n1_tm, int n2_tm, int n3_tm)
}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 )
{
    (double *bdiffmats, int n1_bd, int n2_bd, int n3_bd)
}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2 )
{
    (double *m66, int n1_m66, int n2_m66)
}

// For naff_general
%apply (double* IN_ARRAY2, int DIM1, int DIM2 )
{
    (double* re_in, int n1_re_in, int n2_re_in)
}
%apply (double* IN_ARRAY2, int DIM1, int DIM2 )
{
    (double* im_in, int n1_im_in, int n2_im_in)
}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2 )
{
    (double *ff_out, int n1_ff_out, int n2_ff_out)
}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2 )
{
    (double *re_out, int n1_re_out, int n2_re_out)
}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2 )
{
    (double *im_out, int n1_im_out, int n2_im_out)
}

// For Kicktable::getkicks
%apply double *OUTPUT { double &hkick__, double &vkick__ };

%include "../include/trackcpp/elements.h"
%include "../include/trackcpp/kicktable.h"
%include "../include/trackcpp/auxiliary.h"
%include "../include/trackcpp/accelerator.h"
%include "../include/trackcpp/pos.h"
%include "../include/trackcpp/tracking.h"
%include "../include/trackcpp/optics.h"
%include "../include/trackcpp/diffusion_matrix.h"
%include "../include/trackcpp/naff.h"
%include "../include/trackcpp/linalg.h"
%include "interface.h"

%template(CppDoublePos) Pos<double>;
