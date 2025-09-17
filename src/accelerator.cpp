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

void Accelerator::get_time_aware_elements_info(
  std::vector<unsigned int>& time_aware_indices,
  std::vector<double>& time_aware_dl_kicks,
  unsigned int element_offset
) const {

  // for adjusting dl to keep the arrival-time in sync with the wall clock
  time_aware_indices.clear();
  time_aware_dl_kicks.clear();
  std::vector<double> time_aware_displacements = {};

  size_t nr_elements = this->lattice.size();
  PassMethodsClass PMClass = PassMethodsClass();

  double s_pos = 0.0;

  for (size_t i = 0; i < nr_elements; i++) {
      const Element& element = this->lattice[element_offset];

      if (PMClass.is_time_aware_pm(element.pass_method)) {
        time_aware_indices.push_back(i);
        s_pos += 0.5 * element.length;
        time_aware_displacements.push_back(s_pos);
        s_pos = 0.5 * element.length;
      }
      else {
        s_pos += element.length;
      }
      element_offset = (element_offset + 1) % nr_elements;
  }

  if (time_aware_indices.size() > 0) {
    time_aware_displacements[0] += s_pos;
  }

  //? NOTE : The diference between "line_length" and the sum of "time_aware_displacements" (cum_length) is on the order of
  //? 1e-15 ~ 1e-16. The propagation of this tiny error affects the tracking. The following line avoids losing precision.
  double cum_length = std::accumulate(time_aware_displacements.begin(), time_aware_displacements.end(), 0.0);
  double line_length = this->get_length();
  time_aware_displacements.back() += line_length - cum_length;

  double ddl = 0.0;
  for (size_t i=0; i<time_aware_indices.size(); i++) {
    ddl = light_speed*this->harmonic_number/this->lattice[time_aware_indices[i]].frequency - line_length;
    time_aware_dl_kicks.push_back(ddl * time_aware_displacements[i] / line_length);
  }

}
