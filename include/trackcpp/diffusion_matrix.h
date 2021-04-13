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

#ifndef _DIFFMAT
#define _DIFFMAT

#include "trackcpp.h"
#include "tracking.h"
#include "auxiliary.h"
#include "passmethods.hpp"

// constants for 4th-order symplectic integrator
#define DRIFT1 ( 0.6756035959798286638e00)
#define DRIFT2 (-0.1756035959798286639e00)
#define KICK1  ( 0.1351207191959657328e01)
#define KICK2  (-0.1702414383919314656e01)

// Physical constants used in the calculations
#define TWOPI (6.28318530717959)
#define CGAMMA (8.846056192e-05)  // [m]/[GeV^3] Ref[1] (4.1)
#define M0C2 (5.10999060e5)  // Electron rest mass [eV]
#define LAMBDABAR (3.86159323e-13)  // Compton wavelength/2pi [m]
#define CER (2.81794092e-15)  // Classical electron radius [m]
#define CU (1.323094366892892)  // 55/(24*sqrt(3)) factor


Status::type track_diffusionmatrix (const Accelerator& accelerator,
                                    const Pos<double>& fixed_point,
                                    const std::vector<Matrix>& tm,
                                    std::vector<Matrix>& bmat);

#endif
