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

#ifndef _OPTICS_H
#define _OPTICS_H

#include "elements.h"
#include "pos.h"
#include "accelerator.h"
#include "tracking.h"
#include "linalg.h"

class Twiss {

public:
  Pos<double> co;
  Vector etax;
  Vector etay;
  double mux, betax, alphax;
  double muy, betay, alphay;
  Twiss() : co(std::nan("")),
            etax(Vector({std::nan(""),std::nan("")})),
            etay(Vector({std::nan(""),std::nan("")})),
            mux(0), betax(std::nan("")), alphax(std::nan("")),
            muy(0), betay(std::nan("")), alphay(std::nan("")) {}

  bool isundef() const { return std::isnan(this->betax); }

friend std::ostream& operator<< (std::ostream &out, const Twiss& tw);

};

Status::type calc_twiss(const Accelerator& accelerator,
                        const Pos<double>& fixed_point,
                        Matrix& m66,
                        std::vector<Twiss>& twiss,
                        Twiss twiss0 = Twiss(),
                        bool closed_flag = false);

#endif
