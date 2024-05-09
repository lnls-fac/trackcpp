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
#include <cmath>
#include <algorithm>
 

Accelerator::Accelerator(const double energy) :
    energy(this->_energy),
    gamma_factor(this->_gamma_factor),
    beta_factor(this->_beta_factor),
    velocity(this->_velocity),
    brho(this->_brho) {
    setEnergy(energy);
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
  out << std::endl << "velocity       : " << a.velocity;
  out << std::endl << "beta_factor    : " << a.beta_factor;
  out << std::endl << "gamma_factor   : " << a.gamma_factor;
  return out;
}

void Accelerator::setEnergy(const double _energy_) {

  if (_energy_ <= electron_rest_energy_eV) {

    this->_energy = electron_rest_energy_eV;
    this->_velocity = 0.0;
    this->_gamma_factor = 1.0;
    this->_beta_factor = 0.0;
    this->_brho = 0.0;

  } else {

    this->_energy = _energy_;

    double gamma_ = _energy_ / electron_rest_energy_eV;
    double beta_ = sqrt(1.0 - (1.0 / (gamma_ * gamma_)));
    double velocity_ = beta_ * light_speed;
    double beam_rigidity_ = beta_ * _energy_ / light_speed;
    
    this->_velocity = velocity_;
    this->_gamma_factor = gamma_;
    this->_beta_factor = beta_;
    this->_brho = beam_rigidity_;
  }
  
}

void Accelerator::setGammaFactor(const double _gamma_) {

  if (_gamma_ <= 1.0) {

    this->_energy = electron_rest_energy_eV;
    this->_velocity = 0.0;
    this->_gamma_factor = 1.0;
    this->_beta_factor = 0.0;
    this->_brho = 0.0;

  } else {

    this->_gamma_factor = _gamma_;
    
    double energy_ = _gamma_ * electron_rest_energy_eV;
    double beta_ = sqrt(1.0 - (1.0 / (_gamma_ * _gamma_)));
    double velocity_ = beta_ * light_speed;
    double beam_rigidity_ = beta_ * energy_ / light_speed;
    
    this->_energy = energy_;
    this->_velocity = velocity_;
    this->_beta_factor = beta_;
    this->_brho = beam_rigidity_;
  }
  
}

void Accelerator::setBetaFactor(const double _beta_) {

  if (_beta_ <= 0.0) {

    this->_energy = electron_rest_energy_eV;
    this->_velocity = 0.0;
    this->_gamma_factor = 1.0;
    this->_beta_factor = 0.0;
    this->_brho = 0.0;

  } else {

    this->_beta_factor = _beta_;

    double gamma_ = 1.0 / sqrt(1.0 - (_beta_ * _beta_));
    double energy_ = gamma_ * electron_rest_energy_eV;
    double velocity_ = _beta_ * light_speed;
    double beam_rigidity_ = _beta_ * energy_ / light_speed;
    
    this->_energy = energy_;
    this->_velocity = velocity_;
    this->_gamma_factor = gamma_;
    this->_brho = beam_rigidity_;
  }
  
}

void Accelerator::setVelocity(const double _velocity_) {

  if (_velocity_ <= 0.0) {

    this->_energy = electron_rest_energy_eV;
    this->_velocity = 0.0;
    this->_gamma_factor = 1.0;
    this->_beta_factor = 0.0;
    this->_brho = 0.0;

  } else {

    this->_velocity = _velocity_;

    double beta_ = _velocity_ / light_speed;
    double gamma_ = 1.0 / sqrt(1.0 - (beta_ * beta_));
    double energy_ = gamma_ * electron_rest_energy_eV;
    double beam_rigidity_ = beta_ * energy_ / light_speed;

    this->_energy = energy_;
    this->_gamma_factor = gamma_;
    this->_beta_factor = beta_;
    this->_brho = beam_rigidity_;
  }
  
}

void Accelerator::setMagneticRigidity(const double _brho_) {

  if (_brho_ <= 0.0) {

    this->_energy = electron_rest_energy_eV;
    this->_velocity = 0.0;
    this->_gamma_factor = 1.0;
    this->_beta_factor = 0.0;
    this->_brho = 0.0;

  } else {

    this->_brho = _brho_;

    double k = _brho_ / (electron_rest_energy_eV / light_speed);
    double gamma_ = sqrt(1.0 + (k * k));
    double energy_ = gamma_ * electron_rest_energy_eV;
    double beta_ = sqrt(1.0 - (1.0 / (gamma_ * gamma_)));
    double velocity_ = beta_ * light_speed;

    this->_energy = energy_;
    this->_velocity = velocity_;
    this->_gamma_factor = gamma_;
    this->_beta_factor = beta_;
    
  }
  
}

double Accelerator::time_aware_fraction() const {
  int time_aware_count = 0;
  for (const Element &el : this->lattice) {
    if (el.pass_method == PassMethod::pm_cavity_pass){
      time_aware_count++;
    }
  }
  return 1.0/time_aware_count;
}