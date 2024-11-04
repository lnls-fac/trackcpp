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

#include <trackcpp/trackcpp.h>
#include <trackcpp/tracking.h>
#include <trackcpp/auxiliary.h>


// track_findm66
// -------------
// returns a vector with 6-d transfer matrices, one for each element
//
// inputs:
//    accelerator:  Structure representing the accelerator
//    fixed_point:  Pos representing calculated fixed_point.
//
// outputs:
//    tm:     vector of Matrix6 elements. Each component represents the accumulated
//            transfer matrix from the start of the lattice to the entrance of that element. The last element is the m66 of the line.
//    m66:    one-turn transfer matrix
//
//    v0:     const term of final map
//
//    RETURN:      status do tracking (see 'auxiliary.h')

Status::type track_findm66 (Accelerator& accelerator,
                            const Pos<double>& fixed_point,
                            std::vector<Matrix>& tm,
                            Matrix& m66,
                            Pos<double>& v0,
                            std::vector<unsigned int >& indices,
                            const double line_length,
                            const std::vector<unsigned int>& time_aware_element_indices,
                            const std::vector<double>& time_aware_element_positions) {

  Status::type status  = Status::success;
  const std::vector<Element>& lattice = accelerator.lattice;
  Pos<double> fp = fixed_point;

  const int radsts = accelerator.radiation_on;
  if (radsts == RadiationState::full){
    accelerator.radiation_on = RadiationState::damping;
  }
  // case no closed_orbit has been defined
  if (std::isnan(fp.rx)) fp = Pos<double>(0,0,0,0,0,0);

  int nr_elements  = lattice.size();
	std::vector<bool> indcs;
	indcs.reserve(nr_elements+1);
	for (unsigned int i=0; i<=nr_elements; ++i) indcs[i] = false;
	for (auto&& i: indices) if (i<=nr_elements) indcs[i] = true;

  Pos<Tpsa<6,1> > map;
  map.rx = Tpsa<6,1>(fp.rx, 0); map.px = Tpsa<6,1>(fp.px, 1);
  map.ry = Tpsa<6,1>(fp.ry, 2); map.py = Tpsa<6,1>(fp.py, 3);
  map.de = Tpsa<6,1>(fp.de, 4); map.dl = Tpsa<6,1>(fp.dl, 5);

  tm.clear(); tm.reserve(indices.size());
  unsigned int TAW_pivot = 0;
  double ddl = 0;
  for(unsigned int i=0; i<lattice.size(); ++i) {
    if (indcs[i]){
      Matrix m (6);
      m[0][0] = map.rx.c[1]; m[0][1] = map.rx.c[2]; m[0][2] = map.rx.c[3];
      m[0][3] = map.rx.c[4]; m[0][4] = map.rx.c[5]; m[0][5] = map.rx.c[6];
      m[1][0] = map.px.c[1]; m[1][1] = map.px.c[2]; m[1][2] = map.px.c[3];
      m[1][3] = map.px.c[4]; m[1][4] = map.px.c[5]; m[1][5] = map.px.c[6];
      m[2][0] = map.ry.c[1]; m[2][1] = map.ry.c[2]; m[2][2] = map.ry.c[3];
      m[2][3] = map.ry.c[4]; m[2][4] = map.ry.c[5]; m[2][5] = map.ry.c[6];
      m[3][0] = map.py.c[1]; m[3][1] = map.py.c[2]; m[3][2] = map.py.c[3];
      m[3][3] = map.py.c[4]; m[3][4] = map.py.c[5]; m[3][5] = map.py.c[6];
      m[4][0] = map.de.c[1]; m[4][1] = map.de.c[2]; m[4][2] = map.de.c[3];
      m[4][3] = map.de.c[4]; m[4][4] = map.de.c[5]; m[4][5] = map.de.c[6];
      m[5][0] = map.dl.c[1]; m[5][1] = map.dl.c[2]; m[5][2] = map.dl.c[3];
      m[5][3] = map.dl.c[4]; m[5][4] = map.dl.c[5]; m[5][5] = map.dl.c[6];
    tm.push_back(std::move(m));
    }
    // track through element
    if (i == time_aware_element_indices[TAW_pivot]) {
            ddl = light_speed*accelerator.harmonic_number/lattice[i].frequency - line_length;
            map.dl -= ddl * (time_aware_element_positions[TAW_pivot+1]-time_aware_element_positions[TAW_pivot]) / line_length;
            TAW_pivot++;
    }
    if ((status = track_elementpass (accelerator, lattice[i], map)) != Status::success) return status;

    if (i == time_aware_element_indices.back()) {
            map.dl -= ddl * (time_aware_element_positions[TAW_pivot+1]-time_aware_element_positions[TAW_pivot]) / line_length;
    }
  }

  m66 = Matrix(6);
  Matrix& m = m66;
  m[0][0] = map.rx.c[1]; m[0][1] = map.rx.c[2]; m[0][2] = map.rx.c[3];
  m[0][3] = map.rx.c[4]; m[0][4] = map.rx.c[5]; m[0][5] = map.rx.c[6];
  m[1][0] = map.px.c[1]; m[1][1] = map.px.c[2]; m[1][2] = map.px.c[3];
  m[1][3] = map.px.c[4]; m[1][4] = map.px.c[5]; m[1][5] = map.px.c[6];
  m[2][0] = map.ry.c[1]; m[2][1] = map.ry.c[2]; m[2][2] = map.ry.c[3];
  m[2][3] = map.ry.c[4]; m[2][4] = map.ry.c[5]; m[2][5] = map.ry.c[6];
  m[3][0] = map.py.c[1]; m[3][1] = map.py.c[2]; m[3][2] = map.py.c[3];
  m[3][3] = map.py.c[4]; m[3][4] = map.py.c[5]; m[3][5] = map.py.c[6];
  m[4][0] = map.de.c[1]; m[4][1] = map.de.c[2]; m[4][2] = map.de.c[3];
  m[4][3] = map.de.c[4]; m[4][4] = map.de.c[5]; m[4][5] = map.de.c[6];
  m[5][0] = map.dl.c[1]; m[5][1] = map.dl.c[2]; m[5][2] = map.dl.c[3];
  m[5][3] = map.dl.c[4]; m[5][4] = map.dl.c[5]; m[5][5] = map.dl.c[6];

  if (indcs[lattice.size()]) tm.push_back(m66);

  // constant term of the final map
  v0.rx = map.rx.c[0]; v0.px = map.px.c[0]; v0.ry = map.ry.c[0]; v0.py = map.py.c[0]; v0.de = map.de.c[0]; v0.dl = map.dl.c[0];

  accelerator.radiation_on = radsts;

  return status;

}


Status::type track_findm66 (Accelerator& accelerator,
                            const Pos<double>& fixed_point,
                            std::vector<Matrix>& tm,
                            Matrix& m66,
                            Pos<double>& v0,
                            const double line_length,
        const std::vector<unsigned int>& time_aware_element_indices,
        const std::vector<double>& time_aware_element_positions) {

  std::vector<unsigned int> indices;
  unsigned int nr_elements = accelerator.lattice.size();

  indices.reserve(nr_elements + 1);
	for (unsigned int i=0; i<=nr_elements; ++i) indices.push_back(i);

	return track_findm66 (accelerator, fixed_point, tm, m66, v0, indices, line_length, time_aware_element_indices, time_aware_element_positions);
}


Status::type track_findorbit6(
    Accelerator& accelerator,
    std::vector<Pos<double> >& closed_orbit,
    const Pos<double>& fixed_point_guess) {

  const std::vector<Element>& the_ring = accelerator.lattice;

  double delta        = 1e-9;              // [m],[rad],[dE/E]
  double tolerance    = 2.22044604925e-14;
  int    max_nr_iters = 50;

  const int radsts = accelerator.radiation_on;
  if (radsts == RadiationState::full){
    accelerator.radiation_on = RadiationState::damping;
  }

  // for longitudinal kick before RF cavities
  std::vector<double> TAW_positions;
  std::vector<unsigned int> TAW_indices;
  double accelerator_length = accelerator.get_time_aware_elements_info(TAW_indices, TAW_positions);

  // temporary vectors and matrices
  std::vector<Pos<double> > co(7,0);
  for(auto i=0; i<7; ++i) co[i] = fixed_point_guess;
  std::vector<Pos<double> > co2(7,0);
  std::vector<Pos<double> > D(7,0);
  std::vector<Pos<double> > M(6,0);
  Pos<double> dco(1.0,1.0,1.0,1.0,1.0,1.0);
  matrix6_set_identity_posvec(D, delta);

  int nr_iter = 0;
  while ((get_max(dco) > tolerance) and (nr_iter <= max_nr_iters)) {
    co = co + D;
    Pos<double> Ri = co[6];
    std::vector<Pos<double> > co2;
    unsigned int element_offset = 0;
    Plane::type lost_plane;
    Status::type status = Status::success;

    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[0], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[1], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[2], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[3], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[4], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[5], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[6], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));

    if (status != Status::success) {
      return Status::findorbit_one_turn_matrix_problem;
    }

    Pos<double> Rf = co2[6];

    M[0] = (co2[0] - Rf) / delta;
    M[1] = (co2[1] - Rf) / delta;
    M[2] = (co2[2] - Rf) / delta;
    M[3] = (co2[3] - Rf) / delta;
    M[4] = (co2[4] - Rf) / delta;
    M[5] = (co2[5] - Rf) / delta;

    Pos<double> b = Rf - Ri;
    std::vector<Pos<double> > M_1(6,0);
    matrix6_set_identity_posvec(M_1);
    M_1 = M_1 - M;
    dco = linalg_solve6_posvec(M_1, b);
    co[6] = dco + Ri;
    co[0] = co[6]; co[1] = co[6];
    co[2] = co[6]; co[3] = co[6];
    co[4] = co[6]; co[5] = co[6];
    nr_iter++;
  }

  if (nr_iter > max_nr_iters) {
    return Status::findorbit_not_converged;
  }

  // propagates fixed point throught the_ring
  closed_orbit.clear();
  unsigned int element_offset = 0;
  Plane::type lost_plane;
  track_linepass(
    accelerator, co[6], true, element_offset, closed_orbit, lost_plane, accelerator_length, TAW_indices, TAW_positions
  );
  accelerator.radiation_on = radsts;
  return Status::success;

}

Status::type track_findorbit4(
    Accelerator& accelerator,
    std::vector<Pos<double> >& closed_orbit,
    const Pos<double>& fixed_point_guess) {

  const std::vector<Element>& the_ring = accelerator.lattice;

  double delta        = 1e-9;              // [m],[rad],[dE/E]
  double tolerance    = 2.22044604925e-14;
  int    max_nr_iters = 50;

  const int radsts = accelerator.radiation_on;
  if (radsts == RadiationState::full){
    accelerator.radiation_on = RadiationState::damping;
  }

  // for longitudinal kick before RF cavities
  std::vector<double> TAW_positions;
  std::vector<unsigned int> TAW_indices;
  double accelerator_length = accelerator.get_time_aware_elements_info(TAW_indices, TAW_positions);

  // temporary vectors and matrices
  // std::vector<Pos<double> > co(7,0); // no initial guess
  std::vector<Pos<double> > co(7,fixed_point_guess);
  std::vector<Pos<double> > co2(7,0);
  std::vector<Pos<double> > D(7,0);
  std::vector<Pos<double> > M(6,0);
  Pos<double> dco(1.0,1.0,1.0,1.0,0.0,0.0);
  Pos<double> theta(0.0,0.0,0.0,0.0,0.0,0.0);
  matrix6_set_identity_posvec(D, delta);

  int nr_iter = 0;
  while ((get_max(dco) > tolerance) and (nr_iter <= max_nr_iters)) {
    co = co + D;
    Pos<double> Ri = co[6];
    std::vector<Pos<double> > co2;
    unsigned int element_offset = 0;
    Plane::type lost_plane;
    Status::type status = Status::success;
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[0], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[1], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[2], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[3], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    status = (Status::type) ((int) status | (int) track_linepass(
      accelerator, co[6], false, element_offset, co2, lost_plane, accelerator_length, TAW_indices, TAW_positions
    ));
    if (status != Status::success) {
      return Status::findorbit_one_turn_matrix_problem;
    }
    Pos<double> Rf = co2[4];
    M[0] = (co2[0] - Rf) / delta;
    M[1] = (co2[1] - Rf) / delta;
    M[2] = (co2[2] - Rf) / delta;
    M[3] = (co2[3] - Rf) / delta;
    Pos<double> b = Rf - Ri;
    std::vector<Pos<double> > M_1(6,0);
    matrix6_set_identity_posvec(M_1);
    M_1 = M_1 - M;
    dco = linalg_solve4_posvec(M_1, b);
    co[6] = dco + Ri;
    co[0] = co[6]; co[1] = co[6];
    co[2] = co[6]; co[3] = co[6];
    co[4] = co[6]; co[5] = co[6];
    nr_iter++;
  }

  if (nr_iter > max_nr_iters) {
    return Status::findorbit_not_converged;
  }

  // propagates fixed point throught the_ring
  closed_orbit.clear();
  unsigned int element_offset = 0;
  Plane::type lost_plane;
  track_linepass(
    accelerator, co[6], true, element_offset, closed_orbit, lost_plane, accelerator_length, TAW_indices, TAW_positions
  );
  accelerator.radiation_on = radsts;
  return Status::success;

}
