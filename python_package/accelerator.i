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

class Accelerator {
public:
  Accelerator(const double& energy=-1);        // energy < electron_rest_energy -> energy = electron_rest_energy
  double                  energy;              // [eV]
  bool                    cavity_on;
  bool                    radiation_on;
  bool                    vchamber_on;
  int                     harmonic_number;
  std::vector<Element>    lattice;
  std::vector<Kicktable*> kicktables;

  bool operator==(const Accelerator& o) const;
  bool operator!=(const Accelerator& o) const { return !(*this == o); };
  bool isequal(const Accelerator& a) const { return *this == a; }

};
