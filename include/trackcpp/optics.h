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

#ifndef _OPTICS_H
#define _OPTICS_H

#include "elements.h"
#include "pos.h"
#include "accelerator.h"
#include "tracking.h"

class Twiss {

public:
  Pos<double> cox, coy;
  Pos<double> etax, etay;
  double mux, betax, alphax;
  double muy, betay, alphay;
  Twiss() {
    mux = muy = 0.0;
    betax = betay = std::nan("");
    alphax = alphay = std::nan("");
    etax.de = etay.de = 1.0;  // so that (etax * de) gives orbit with delta energy
  }
};

Status::type calc_twiss(const Accelerator& accelerator, const Pos<double>& fixed_point, Matrix& m66, std::vector<Twiss>& twiss);

#endif
