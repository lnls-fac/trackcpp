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

#ifndef _FLAT_FILE_H
#define _FLAT_FILE_H

#include "accelerator.h"
#include "elements.h"
#include "auxiliary.h"
#include <string>
#include <vector>
#include <iomanip>

struct FlatFileType {
    enum type_ {
      marker    = -1,
      drift     =  0,
      mpole     =  1,
      cavity    =  2,
      corrector =  3,
      thin_kick =  3,
      kicktable =  6
    };
  };

//Status::type read_flat_file(const std::string& filename, Accelerator& accelerator);

Status::type read_flat_file(std::string& filename, Accelerator& accelerator, bool file_flag = true);
Status::type write_flat_file(std::string& filename, const Accelerator& accelerator, bool file_flag = true);

#endif
