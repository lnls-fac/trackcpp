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

#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "dynap.h"
#include "accelerator.h"
#include "elements.h"
#include <vector>
#include <string>
#include <cstdlib>

Status::type print_closed_orbit      (const Accelerator& accelerator, const std::vector<Pos<double>>&    cod,  const std::string& filename = "cod_out.txt");
Status::type print_dynapgrid         (const Accelerator& accelerator, const std::vector<DynApGridPoint>& grid, const std::string& label, const std::string& filename = "dynap_out.txt", bool print_tunes=false);
Status::type print_tracking_ringpass (const Accelerator& accelerator, const std::vector<Pos<double>>& points, const std::string& filename = "track_linepass_out.txt");
Status::type print_tracking_linepass (const Accelerator& accelerator, const std::vector<Pos<double>>& points, const unsigned int start_element, const std::string& filename);
void         print_header            (FILE* fp);

#endif
