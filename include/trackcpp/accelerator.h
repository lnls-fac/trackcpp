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

#ifndef _ACCELERATOR_H
#define _ACCELERATOR_H

#include "kicktable.h"
#include "elements.h"
#include <vector>

//struct Accelerator {
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
  bool isequal(const Accelerator& a) const { return *this == a; } // necessary for python_package

  friend std::ostream& operator<< (std::ostream &out, const Accelerator& a);

};

#endif
