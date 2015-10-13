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
