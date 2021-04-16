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

// Trackcpp passmethods are based on SLAC Andrei Terebilo AT version 1.3
// <http://www.slac.stanford.edu/grp/ssrl/spear/at/>.

#ifndef _PASS_METHODS_H
#define _PASS_METHODS_H

#include "accelerator.h"
#include "pos.h"
#include "auxiliary.h"

#define ATCOMPATIBLE 1

// constants for 4th-order symplectic integrator
const double DRIFT1 = 0.6756035959798286638e00;
const double DRIFT2 = -0.1756035959798286639e00;
const double KICK1 = 0.1351207191959657328e01;
const double KICK2 = -0.1702414383919314656e01;

// Physical constants used in the calculations
#ifdef ATCOMPATIBLE
  const double TWOPI = 6.28318530717959;   // 2*pi
  const double CGAMMA = 8.846056192e-05;   // cgamma, [m]/[GeV^3] Ref[1] (4.1)
  const double M0C2 = 5.10999060e5;        // Electron rest mass [eV]
  const double LAMBDABAR = 3.86159323e-13; // Compton wavelength/2pi [m]
  const double CER = 2.81794092e-15;       // Classical electron radius [m]
  const double CU = 1.323094366892892;     // 55/(24*sqrt(3)) factor
#else
  const double TWOPI_ = 2*M_PI;
  const double CGAMMA_ = 4*M_PI*electron_radius/pow(electron_rest_energy/electron_charge/1e9,3)/3;
  const double M0C2 = electron_rest_energy_MeV * 1e6;
  const double LAMBDABAR = reduced_planck_constant / light_speed / electron_mass;
  const double CER = electron_radius;
  const double CU = 55/(24*std::sqrt(3));
#endif


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
template <typename T> Status::type pm_matrix_pass                (Pos<T> &pos, const Element &elem, const Accelerator& accelerator);

#include "passmethods.hpp"

#endif
