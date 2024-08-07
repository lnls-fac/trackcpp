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
  int lost_turn;
  int lost_element;
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

// Status::type dynap_ma(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& e0,
//     const double& e_tol,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads
//   );

Status::type dynap_acceptance(
    const std::string calc_type,
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

Status::type dynap_ma(
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

// Status::type dynap_pxa(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& px0,
//     const double& px_tol,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads
//   );

Status::type dynap_pxa(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& px_init,
    const double& px_delta,
    unsigned int nr_steps_back,
    double rescale,
    unsigned int nr_iterations,
    const double& s_min, const double& s_max,
    const std::vector<std::string>& fam_names,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  );

// Status::type dynap_pya(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& py0,
//     const double& py_tol,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads
//   );

Status::type dynap_pya(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& py_init,
    const double& py_delta,
    unsigned int nr_steps_back,
    double rescale,
    unsigned int nr_iterations,
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
