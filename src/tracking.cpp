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

#include <trackcpp/trackcpp.h>
#include <trackcpp/auxiliary.h>


// track_findm66
// -------------
// returns a vector with 6-d transfer matrices, one for each element
//
// inputs:
//    accelerator: strucrure representing the accelerator
//    closed_orbit:  Pos vector representing calculated closed orbit.
//
// outputs:
//    tm:     vector of Matrix6 elements. Each component represents the accumulated
//            transfer matrix from the start of the lattice to the entrance of that element.
//    m66:    one-turn transfer matrix
//
//    RETURN:      status do tracking (see 'auxiliary.h')

//Status::type track_findm66 (const Accelerator& accelerator, std::vector<Pos<double> >& closed_orbit, std::vector<Matrix>& tm, Matrix& m66) {
Status::type track_findm66 (const Accelerator& accelerator, std::vector<Pos<double> >& closed_orbit, std::vector<Matrix>& tm) {

  Status::type status  = Status::success;
  const std::vector<Element>& lattice = accelerator.lattice;

  tm.clear();

  std::vector<double> row0 = {0,0,0,0,0,0};

  // case no closed_orbit has been defined
  if (closed_orbit.size() != lattice.size()) {
    closed_orbit.clear();
    for(unsigned int i=0; i<lattice.size(); ++i) {
      closed_orbit.push_back(Pos<double>(0,0,0,0,0,0));
    }
  }

  Pos<Tpsa<6,1> > map;
  map.rx = Tpsa<6,1>(closed_orbit[0].rx, 0); map.px = Tpsa<6,1>(closed_orbit[0].px, 1);
  map.ry = Tpsa<6,1>(closed_orbit[0].ry, 2); map.py = Tpsa<6,1>(closed_orbit[0].py, 3);
  map.de = Tpsa<6,1>(closed_orbit[0].de, 4); map.dl = Tpsa<6,1>(closed_orbit[0].dl, 5);

  for(unsigned int i=0; i<lattice.size(); ++i) {

    // track through element
    if ((status = track_elementpass (lattice[i], map, accelerator)) != Status::success) return status;

    Matrix m = {row0,row0,row0,row0,row0,row0};

    m[0][0] = map.rx.c[1]; m[0][1] = map.rx.c[2];
    m[0][2] = map.rx.c[3]; m[0][3] = map.rx.c[4];
    m[0][4] = map.rx.c[5]; m[0][5] = map.rx.c[6];

    m[1][0] = map.px.c[1]; m[1][1] = map.px.c[2];
    m[1][2] = map.px.c[3]; m[1][3] = map.px.c[4];
    m[1][4] = map.px.c[5]; m[1][5] = map.px.c[6];

    m[2][0] = map.ry.c[1]; m[2][1] = map.ry.c[2];
    m[2][2] = map.ry.c[3]; m[2][3] = map.ry.c[4];
    m[2][4] = map.ry.c[5]; m[2][5] = map.ry.c[6];

    m[3][0] = map.py.c[1]; m[3][1] = map.py.c[2];
    m[3][2] = map.py.c[3]; m[3][3] = map.py.c[4];
    m[3][4] = map.py.c[5]; m[3][5] = map.py.c[6];

    m[4][0] = map.de.c[1]; m[4][1] = map.de.c[2];
    m[4][2] = map.de.c[3]; m[4][3] = map.de.c[4];
    m[4][4] = map.de.c[5]; m[4][5] = map.de.c[6];

    m[5][0] = map.dl.c[1]; m[5][1] = map.dl.c[2];
    m[5][2] = map.dl.c[3]; m[5][3] = map.dl.c[4];
    m[5][4] = map.dl.c[5]; m[5][5] = map.dl.c[6];

    tm.push_back(m);



  }

  return status;

}

Status::type track_findorbit6(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& closed_orbit,
    const Pos<double>& fixed_point_guess) {

  const std::vector<Element>& the_ring = accelerator.lattice;

  double delta        = 1e-9;              // [m],[rad],[dE/E]
  double tolerance    = 2.22044604925e-14;
  int    max_nr_iters = 50;

  // calcs longitudinal fixed point
  double L0 = latt_findspos(the_ring, 1+the_ring.size());
  double T0 = L0 / light_speed;
  std::vector<int>    cav_idx = latt_findcells_frequency(the_ring, 0, true);
  double frf = the_ring[cav_idx[0]].frequency;
  double fixedpoint = light_speed*((1.0*accelerator.harmonic_number)/frf - T0);

  // temporary vectors and matrices
  // std::vector<Pos<double> > co(7,0); // no initial guess
  std::vector<Pos<double> > co(7,fixed_point_guess);
  std::vector<Pos<double> > co2(7,0);
  std::vector<Pos<double> > D(7,0);
  std::vector<Pos<double> > M(6,0);
  Pos<double> dco(1.0,1.0,1.0,1.0,1.0,1.0);
  Pos<double> theta(0.0,0.0,0.0,0.0,0.0,0.0);
  theta.dl = fixedpoint;
  matrix6_set_identity_posvec(D, delta);

  int nr_iter = 0;
  while ((get_max(dco) > tolerance) and (nr_iter <= max_nr_iters)) {
    co = co + D;
    Pos<double> Ri = co[6];
    std::vector<Pos<double> > co2;
    unsigned int element_offset = 0;
    Plane::type lost_plane;
    Status::type status = Status::success;

    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[0], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[1], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[2], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[3], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[4], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[5], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[6], co2, element_offset, lost_plane, false));

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

    Pos<double> b = Rf - Ri - theta;
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
  track_linepass(accelerator, co[6], closed_orbit, element_offset, lost_plane, true);
  closed_orbit.pop_back(); // eliminates last element which is the same as first
  return Status::success;

}

Status::type track_findorbit4(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& closed_orbit,
    const Pos<double>& fixed_point_guess) {

  const std::vector<Element>& the_ring = accelerator.lattice;

  double delta        = 1e-9;              // [m],[rad],[dE/E]
  double tolerance    = 2.22044604925e-14;
  int    max_nr_iters = 50;

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
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[0], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[1], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[2], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[3], co2, element_offset, lost_plane, false));
    //status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[4], co2, element_offset, lost_plane, false));
    //status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[5], co2, element_offset, lost_plane, false));
    status = (Status::type) ((int) status | (int) track_linepass(accelerator, co[6], co2, element_offset, lost_plane, false));
    if (status != Status::success) {
      //printf("nr_iter: %i\n", nr_iter);
      //printf("element: %i\n", element_offset);
      //printf("plane: %i\n", lost_plane);
      return Status::findorbit_one_turn_matrix_problem;
    }
    //print(co2);
    Pos<double> Rf = co2[4];
    M[0] = (co2[0] - Rf) / delta;
    M[1] = (co2[1] - Rf) / delta;
    M[2] = (co2[2] - Rf) / delta;
    M[3] = (co2[3] - Rf) / delta;
    //M[4] = (co2[4] - Rf) / delta;
    //M[5] = (co2[5] - Rf) / delta;
    //print(M);
    Pos<double> b = Rf - Ri;
    std::vector<Pos<double> > M_1(6,0);
    matrix6_set_identity_posvec(M_1);
    M_1 = M_1 - M;
    dco = linalg_solve4_posvec(M_1, b);
    //printf("%.4e %.4e %.4e %.4e %.4e %.4e\n", Rf.rx, Rf.px, Rf.ry, Rf.py, Rf.de, Rf.dl);
    //printf("%.4e %.4e %.4e %.4e %.4e %.4e\n", b.rx, b.px, b.ry, b.py, b.de, b.dl);
    //printf("%.4e %.4e %.4e %.4e %.4e %.4e\n", dco.rx, dco.px, dco.ry, dco.py, dco.de, dco.dl);
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
  track_linepass(accelerator, co[6], closed_orbit, element_offset, lost_plane, true);
  closed_orbit.pop_back(); // eliminates last element which is the same as first
  return Status::success;

}
