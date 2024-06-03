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

#ifndef _ACCELERATOR_H
#define _ACCELERATOR_H

#include "kicktable.h"
#include "elements.h"
#include <vector>
#include <string>

//struct Accelerator {
class Accelerator {
public:
  // energy < electron_rest_energy_eV -> energy = electron_rest_energy_eV:
  Accelerator(const double& energy=-1);
  double                  energy;              // [eV]
  bool                    cavity_on = false;
  int                     radiation_on = RadiationState::off;
  bool                    vchamber_on = false;
  int                     harmonic_number = 0;
  std::vector<Element>    lattice;
  std::string             lattice_version = "";

  bool operator==(const Accelerator& o) const;
  bool operator!=(const Accelerator& o) const { return !(*this == o); };
  bool isequal(const Accelerator& a) const { return *this == a; } // necessary for python_package
  double get_length() const;
  double get_time_aware_frac() const;
  friend std::ostream& operator<< (std::ostream &out, const Accelerator& a);
};

#endif
