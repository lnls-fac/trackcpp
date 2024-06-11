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

#include <trackcpp/accelerator.h>
#include <trackcpp/auxiliary.h>

Accelerator::Accelerator(const double& energy) {
  this->energy = (energy < electron_rest_energy_eV) ? electron_rest_energy_eV : energy;
}

double Accelerator::get_length() const {
  double length = 0.0;
  for(auto i=0; i<lattice.size(); ++i) length += lattice[i].length;
  return length;
}

bool Accelerator::operator==(const Accelerator& o) const {

  if (this->energy != o.energy) return false;
  if (this->cavity_on != o.cavity_on) return false;
  if (this->radiation_on != o.radiation_on) return false;
  if (this->vchamber_on != o.vchamber_on) return false;
  if (this->harmonic_number != o.harmonic_number) return false;
  if (this->lattice != o.lattice) return false;
  if (this->lattice_version != o.lattice_version) return false;

  return true;

}


std::ostream& operator<< (std::ostream &out, const Accelerator& a) {
  out <<              "energy         : " << a.energy;
  out << std::endl << "cavity_on      : " << a.cavity_on;
  out << std::endl << "radiation_on   : " << a.radiation_on;
  out << std::endl << "vchamber_on    : " << a.vchamber_on;
  out << std::endl << "harmonic_number: " << a.harmonic_number;
  out << std::endl << "lattice        : " << a.lattice.size() << " elements";
  out << std::endl << "lattice_version: " << a.lattice_version;
  return out;
}

double Accelerator::get_time_aware_elements_info(std::vector<unsigned int>& TAW_indices, std::vector<double>& TAW_positions) const {
  // for longitudinal kick before RF cavities
  TAW_positions.clear();
  TAW_indices.clear();
  size_t nr_elements = this->lattice.size();
  PassMethodsClass PMClass = PassMethodsClass();
  double s_pos = 0.0;
  TAW_positions.push_back(s_pos);
  for (size_t i = 0; i < nr_elements; i++) {
      auto elem = this->lattice[i];
      if (PMClass.is_time_aware_pm(elem.pass_method)) {
          s_pos += (elem.length*0.5);
          TAW_positions.push_back(s_pos);
          TAW_indices.push_back(i);
          s_pos += (elem.length*0.5);
      }
      else {
          s_pos += elem.length;
      }
  }
  TAW_positions.push_back(s_pos);
  return s_pos;
}
