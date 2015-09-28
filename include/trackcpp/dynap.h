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

#ifndef _DYNAP_H
#define _DYNAP_H

#include "accelerator.h"
#include "elements.h"
#include "pos.h"
#include "auxiliary.h"
#include <vector>


struct DynApGridPoint {
  Pos<double>  p;              // Dislocation around closed-orbit
  unsigned int start_element;
  unsigned int lost_turn;
  unsigned int lost_element;
  Plane::type  lost_plane;     // Plane::no_plane,Plane::x,Plane::y,Plane::z
  double       nux1, nuy1;     // tunes at first half number of turns
  double       nux2, nuy2;     // tunes at second half number of turns
};

Status::type dynap_xy(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_x, double x_min, double x_max,
    unsigned int nrpts_y, double y_min, double y_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  );

Status::type dynap_ex(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_e, double e_min, double e_max,
    unsigned int nrpts_x, double x_min, double x_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  );

Status::type dynap_ma(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& e0,
    const double& e_tol,
    const double& s_min, const double& s_max,
    const std::vector<std::string>& fam_names,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  );

Status::type dynap_ma2(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& e_init,
    const double& e_delta,
    unsigned int nr_steps_back,
    double rescale,
    unsigned int nr_iterations,
    const double& s_min, const double& s_max,
    const std::vector<std::string>& fam_names,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  );




Status::type dynap_pxa(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& px0,
    const double& px_tol,
    const double& s_min, const double& s_max,
    const std::vector<std::string>& fam_names,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  );

  Status::type dynap_pya(
      const Accelerator& accelerator,
      std::vector<Pos<double> >& cod,
      unsigned int nr_turns,
      const Pos<double>& p0,
      const double& py0,
      const double& py_tol,
      const double& s_min, const double& s_max,
      const std::vector<std::string>& fam_names,
      bool calculate_closed_orbit,
      std::vector<DynApGridPoint>& grid,
      unsigned int nr_threads
    );

Status::type dynap_xyfmap(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_x, double x_min, double x_max,
    unsigned int nrpts_y, double y_min, double y_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  );

Status::type dynap_exfmap(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_e, double e_min, double e_max,
    unsigned int nrpts_x, double x_min, double x_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  );


#endif
