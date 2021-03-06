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

Status::type track_diffusionmatrix (const Accelerator& accelerator,
                                    const Pos<double>& fixed_point,
                                    const std::vector<Matrix>& tm,
                                    std::vector<Matrix>& bmat);

#endif
