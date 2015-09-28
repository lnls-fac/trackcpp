%module trackcpp

%{
#include <trackcpp/elements.h>
#include <trackcpp/kicktable.h>
#include <trackcpp/auxiliary.h>
#include <trackcpp/accelerator.h>
#include <trackcpp/pos.h>
#include <trackcpp/tracking.h>
#include <trackcpp/flat_file.h>
#include "interface.h"
#include "elementswrapper.h"
%}
%include "carrays.i"
%include "std_string.i"
%include "std_vector.i"
%include "stl.i"

namespace std {
    %template(CppStringVector) vector<string>;
    %template(CppDoubleVector) vector<double>;
    //%template(CppDoublePVector) vector<double*>;
    %template(CppElementVector) vector<Element>;
    %template(CppDoublePosVector) vector< Pos<double> >;
    %template(CppDoubleMatrix) vector< vector<double> >;
    %template(CppDoubleMatrixVector) vector< vector< vector<double> > >;
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

%include "../include/trackcpp/kicktable.h"
%include "../include/trackcpp/auxiliary.h"
%include "../include/trackcpp/pos.h"
%include "../include/trackcpp/tracking.h"
%include "../include/trackcpp/flat_file.h"
%include "interface.h"
%include "elementswrapper.h"

%template(CppDoublePos) Pos<double>;
//%template(double_track_elementpass) track_elementpass<double>;
