// The MIT License (MIT)
//
// Copyright (c) 2015 LNLS Accelerator Division
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <trackcpp/accelerator.h>
#include <trackcpp/auxiliary.h>

Accelerator::Accelerator(const double& energy) {
  this->energy = (energy < electron_rest_energy_MeV*1e6) ? electron_rest_energy_MeV*1e6 : energy;
}

bool Accelerator::operator==(const Accelerator& o) const {

  double                  energy;              // [eV]
  bool                    cavity_on;
  bool                    radiation_on;
  bool                    vchamber_on;
  int                     harmonic_number;
  std::vector<Element>    lattice;
  std::vector<Kicktable*> kicktables;

  if (this->energy != o.energy) return false;
  if (this->cavity_on != o.cavity_on) return false;
  if (this->radiation_on != o.radiation_on) return false;
  if (this->vchamber_on != o.vchamber_on) return false;
  if (this->harmonic_number != o.harmonic_number) return false;
  if (this->lattice != o.lattice) return false;

  // in the element-by-element comparison of the lattice kicktables
  // that are actually being referenced have already been compared.

  return true;

}


std::ostream& operator<< (std::ostream &out, const Accelerator& a) {
  out <<              "energy         : " << a.energy;
  out << std::endl << "cavity_on      : " << a.cavity_on;
  out << std::endl << "radiation_on   : " << a.radiation_on;
  out << std::endl << "vchamber_on    : " << a.vchamber_on;
  out << std::endl << "harmonic_number: " << a.harmonic_number;
  out << std::endl << "lattice        : " << a.lattice.size() << " elements";
  out << std::endl << "kicktables     : " << a.kicktables.size() << " elements";
  return out;
}
