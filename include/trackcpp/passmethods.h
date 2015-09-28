// The MIT License (MIT)
//
// Copyright (c) 2015 LNLS Accelerator Division
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef _PASS_METHODS_H
#define _PASS_METHODS_H

#include "accelerator.h"
#include "pos.h"
#include "auxiliary.h"

double get_magnetic_rigidity(const double energy);

template <typename T> Status::type pm_identity_pass              (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);
template <typename T> Status::type pm_drift_pass                 (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);
template <typename T> Status::type pm_str_mpole_symplectic4_pass (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);
template <typename T> Status::type pm_bnd_mpole_symplectic4_pass (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);
template <typename T> Status::type pm_corrector_pass             (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);
template <typename T> Status::type pm_cavity_pass                (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);
template <typename T> Status::type pm_thinquad_pass              (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);
template <typename T> Status::type pm_thinsext_pass              (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);
template <typename T> Status::type pm_kicktable_pass             (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);

#include "passmethods.hpp"

#endif
