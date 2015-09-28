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

#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "dynap.h"
#include "accelerator.h"
#include "elements.h"
#include <vector>
#include <string>
#include <cstdlib>

Status::type print_closed_orbit      (const Accelerator& accelerator, const std::vector<Pos<double>>&    cod,  const std::string& filename = "cod_out.txt");
Status::type print_dynapgrid         (const Accelerator& accelerator, const std::vector<DynApGridPoint>& grid, const std::string& label, const std::string& filename = "dynap_out.txt", bool print_tunes=false);
Status::type print_tracking_ringpass (const Accelerator& accelerator, const std::vector<Pos<double>>& points, const std::string& filename = "track_linepass_out.txt");
Status::type print_tracking_linepass (const Accelerator& accelerator, const std::vector<Pos<double>>& points, const unsigned int start_element, const std::string& filename);
void         print_header            (FILE* fp);

#endif
