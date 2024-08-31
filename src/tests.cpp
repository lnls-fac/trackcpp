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

#include "commands.h"
#include <trackcpp/trackcpp.h>
#include <ctime>
#include <chrono>
#include <cstdlib>

int test_printlattice(const Accelerator& accelerator) {
  latt_print(accelerator.lattice);
  return 0;
}

int test_linepass(const Accelerator& accelerator) {

  const std::vector<Element>& the_ring = accelerator.lattice;

  Pos<> pos;
  pos.rx = 1*0.00100; pos.px = 0*0.00001;
  pos.ry = 1*0.00010; pos.py = 0*0.00001;

  std::vector<Pos<> > new_pos;
  unsigned int element_offset = 0;
  Plane::type lost_plane;
  Status::type status = track_linepass(accelerator, pos, new_pos, element_offset, lost_plane, true, 0, {0,}, {0.0, 0.0, });
  std::cout << "status: " << string_error_messages[status] << std::endl;


  FILE *fp;
  fp = fopen("orbit_trackcpp.txt", "w");
  for(unsigned int i=1-1; i<new_pos.size(); ++i) {
  //for(unsigned int i=2000; i<292; ++i) {
    const Pos<>& c = new_pos[i];
    fprintf(stdout, "%03i: %15s  %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E \n", i+1, the_ring[i % the_ring.size()].fam_name.c_str(), c.rx, c.px, c.ry, c.py, c.de, c.dl);
    fprintf(fp, "%+23.16E %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E \n", c.rx, c.px, c.ry, c.py, c.de, c.dl);
  }
  fclose(fp);

  return 0;

}

int test_linepass_tpsa(const Accelerator& accelerator, const std::vector<Element>& the_ring) {

  const int order = 3;
  Pos<Tpsa<6,order> > tpsa;

  tpsa.rx = Tpsa<6,order>(0, 0); tpsa.px = Tpsa<6,order>(0, 1);
  tpsa.ry = Tpsa<6,order>(0, 2); tpsa.py = Tpsa<6,order>(0, 3);
  tpsa.de = Tpsa<6,order>(0, 4); tpsa.dl = Tpsa<6,order>(0, 5);
  std::vector<Pos<Tpsa<6,order> > > new_tpsa;
  unsigned int element_offset = 0;
  Plane::type lost_plane;
  track_linepass(accelerator, tpsa, false, element_offset, new_tpsa, lost_plane, 0, {0,}, {0.0, 0.0, });
  for(unsigned int i=0; i<new_tpsa.size(); ++i) {
    //const Pos<Tpsa<6,1> >& c = new_particles[i];
    //std::cout << c.rx << std::endl;
    //fprintf(stdout, "%03i: %15s  %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E\n", i+1, the_ring[i].fam_name.c_str(), c.rx, c.px, c.ry, c.py, c.de);
  }

  return 0;

}

int test_ringpass(const Accelerator& accelerator) {


  Pos<> pos;
  pos.rx = 0.00100;
  pos.ry = 0.00010;

  std::vector<Pos<double> > new_pos;
  unsigned int element_offset = 0
  unsigned int lost_turn = 0;
  Plane::type lost_plane;

  clock_t begin, end;
  double time_spent;
  begin = clock();
  Status::type status = track_ringpass(
    accelerator,
    pos,
    5000,
    true,
    element_offset,
    new_pos,
    lost_plane,
    lost_turn
  );
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  if (status != Status::success) {
    std::cerr << "problem" << std::endl;
  }

  std::cout << "tracking_time: " << time_spent << std::endl;
  for(unsigned int i=0; i<new_pos.size(); ++i) {
    const Pos<>& c = new_pos[i];
    fprintf(stdout, "%03i: %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E\n", i+1, c.rx, c.px, c.ry, c.py, c.de, c.dl);
  }


  return 0;

}

int test_findorbit4() {

  Accelerator accelerator;
  std::string fname("/home/fac_files/code/trackcpp/tests/si_v07_c05.txt");
  read_flat_file(fname, accelerator);
  accelerator.cavity_on = false;
  accelerator.radiation_on = RadiationState::off;
  accelerator.vchamber_on = false;

  accelerator.lattice[10].polynom_b[0] = 1e-3;

  std::vector<Pos<double> > orbit;
  Status::type status = track_findorbit4(accelerator, orbit);
  if (status != Status::success) {
    std::cerr << string_error_messages[status] << std::endl;
  } else {
    const Pos<>& c = orbit[0];
    fprintf(stdout, "closed_orbit: %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E\n", c.rx, c.px, c.ry, c.py, c.de);
  }
  return 0;

}

int test_findorbit6() {

  Accelerator accelerator;
  std::string fname("/home/fac_files/code/trackcpp/tests/si_v07_c05.txt");
  read_flat_file(fname, accelerator);
  accelerator.cavity_on = true;
  accelerator.radiation_on = RadiationState::damping;
  accelerator.vchamber_on = false;

  std::vector<Pos<double> > orbit;
  Status::type status = track_findorbit6(accelerator, orbit);
  if (status != Status::success) {
    std::cerr << string_error_messages[status] << std::endl;
  } else {
    const Pos<>& c = orbit[0];
    fprintf(stdout, "closed_orbit: %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E\n", c.rx, c.px, c.ry, c.py, c.de);
  }
  return 0;

}

int test_cmd_track_linepass() {

  std::vector<std::string> args = {
      "trackcpp",
      "track_linepass",
      "/home/ximenes/pytrack/sirius_v500_ac10_5_bare_in_TRACY.txt",
      "3e9",
      "864",
      "off",
      "off",
      "off",
      "0",
      "0.001",
      "0.0",
      "0.0",
      "0.0",
      "0.0",
      "0.0"
  };

  return cmd_track_linepass(args);

}

int test_cmd_dynap_xy() {

  std::vector<std::string> args = {
      "trackcpp",
      "dynap_xy_threads",
      "sirius-v10.txt",
      "3e9",
      "864",
      "on",
      "on",
      "on",
      "0.0",
      "5000",
      "4",
      "-0.015",
      "+0.015",
      "4",
      "0.0",
      "+0.0035",
      "4"
  };
  return cmd_dynap_xy(args);

}

int test_cmd_dynap_ex() {

  std::vector<std::string> args = {
      "trackcpp",
      "dynap_ex_threads",
      "sirius-v10.txt",
      "3e9",
      "864",
      "on",
      "on",
      "on",
      "0.001",
      "5000",
      "4",
      "-0.05",
      "+0.05",
      "4",
      "-0.015",
      "0.0",
      "4"
  };
  return cmd_dynap_ex(args);

}

int test_cmd_dynap_ma() {

  std::vector<std::string> args = {
      "trackcpp",
      "dynap_ma",
      "sirius-v10.txt",
      "3e9",
      "864",
      "on",
      "on",
      "on",
      "50",      // nr_turns
      "30e-6",   // y0
      "0.01",    // e_init
      "0.005",   // e_delta
      "1",       // nr_steps_back
      "0.2",     // rescale
      "3",       // nr_iterations
      "0",       // s_min [m]
      "200",     // s_max [m]
      "8",       // nr_threads
      "qf2", "qfa", "sda"
       };

  return cmd_dynap_acceptance(args);

}

int test_cmd_dynap_xyfmap() {

  std::vector<std::string> args = {
      "trackcpp",
      "dynap_xyfmap",
      "sirius-v10.txt", // flatfile
      "3e9",      // ebeam energy [eV]
      "864",      // harmonic_number
      "off",      // radiation_state
      "off",      // cavity_state
      "on",       // chamber_state
      "0.0",      // de
      "5004",     // nr_turns
      "4",        // nrpts_x
      "-0.015",   // x_min
      "+0.015",   // x_max
      "4",        // nrpts_y
      "0.0",      // y_min
      "+0.0035",  // y_max
      "8"         // nr_threads  (0: let routine decide)
  };
  return cmd_dynap_xyfmap(args);

}

int test_cmd_dynap_exfmap() {

  std::vector<std::string> args = {
      "trackcpp",
      "dynap_exfmap",
      "sirius-v10.txt", // flatfile
      "3e9",      // ebeam energy [eV]
      "864",      // harmonic_number
      "off",      // radiation_state
      "off",      // cavity_state
      "on",       // chamber_state
      "0.001",    // y [m]
      "5004",     // nr_turns
      "4",        // nrpts_e
      "0.0",      // e_min
      "0.05",     // e_max
      "4",        // nrpts_x
      "-0.015",   // x_min [m]
      "+0.015",   // x_max [m]
      "8"         // nr_threads  (0: let routine decide)
  };
  return cmd_dynap_exfmap(args);

}

int test_kicktable(Accelerator& accelerator) {

  Kicktable t;
  const Kicktable *ptrKicktable = nullptr;
  add_kicktable("/home/fac_files/code/python/trackcpp/pytrack/id_kicktable.txt");
  add_kicktable("/home/fac_files/code/python/trackcpp/pytrack/id_kicktable2.txt");
  return 0;

}

int test_read_flat_file(Accelerator& accelerator) {

  std::string fname("/home/ximenes/flatfile.txt");
  read_flat_file(fname, accelerator);
  std::cout << "nr_elements: " << accelerator.lattice.size() << std::endl;
  return 0;

}

int test_simple_drift() {

  Accelerator accelerator;

  accelerator.energy = 3e9; // [ev]
  accelerator.harmonic_number = 864;
  accelerator.radiation_on = RadiationState::off;
  accelerator.cavity_on = false;
  accelerator.vchamber_on = false;

  Element ds = Element::drift("ds", 1.0);

  accelerator.lattice.push_back(ds);

  Pos<double> pos(0.001,0.002,0.003,0.004,0.005,0.006);
  track_elementpass (accelerator, ds, pos);

  fprintf(stdout, "test_simple_drift\n");
  fprintf(stdout, "rx: %+.16f\n", pos.rx);
  fprintf(stdout, "px: %+.16f\n", pos.px);
  fprintf(stdout, "ry: %+.16f\n", pos.ry);
  fprintf(stdout, "py: %+.16f\n", pos.py);
  fprintf(stdout, "de: %+.16f\n", pos.de);
  fprintf(stdout, "dl: %+.16f\n", pos.dl);

  return EXIT_SUCCESS;
}

int test_simple_quadrupole() {

  Accelerator accelerator;

  accelerator.energy = 3e9; // [ev]
  accelerator.harmonic_number = 864;
  accelerator.radiation_on = RadiationState::off;
  accelerator.cavity_on = false;
  accelerator.vchamber_on = false;

  Element ds = Element::quadrupole("qs", 1.0, 2.0);

  accelerator.lattice.push_back(ds);

  Pos<double> pos(0.001,0.002,0.003,0.004,0.005,0.006);
  track_elementpass (accelerator, ds, pos);

  fprintf(stdout, "test_simple_quadrupole\n");
  fprintf(stdout, "rx: %+.16f\n", pos.rx);
  fprintf(stdout, "px: %+.16f\n", pos.px);
  fprintf(stdout, "ry: %+.16f\n", pos.ry);
  fprintf(stdout, "py: %+.16f\n", pos.py);
  fprintf(stdout, "de: %+.16f\n", pos.de);
  fprintf(stdout, "dl: %+.16f\n", pos.dl);

  return EXIT_SUCCESS;
}

int test_linepass2() {
  Accelerator accelerator;
  //sirius_v500(accelerator.lattice);
  accelerator.energy = 3e9; // [ev]
  accelerator.harmonic_number = 864;
  accelerator.radiation_on = RadiationState::off;
  accelerator.cavity_on = false;
  accelerator.vchamber_on = false;

  Pos<double> orig_pos;
  std::vector<Pos<double>> pos;
  unsigned int element_offset = 0;
  Plane::type lost_plane;
  bool trajectory = true;

  orig_pos.rx = 0.0001;
  orig_pos.px = 0.0001;

  Status::type status = track_linepass (
      accelerator,
      orig_pos,
      trajectory,
      element_offset,
      pos,
      lost_plane,
      0,
      {0,},
      {0.0, 0.0, }
  );

      for(unsigned int i=0; i<10; ++i) {
        std::cout << std::endl;
        fprintf(stdout, "rx: %+.16f\n", pos[i].rx);
        fprintf(stdout, "px: %+.16f\n", pos[i].px);
        fprintf(stdout, "ry: %+.16f\n", pos[i].ry);
        fprintf(stdout, "py: %+.16f\n", pos[i].py);
        fprintf(stdout, "de: %+.16f\n", pos[i].de);
        fprintf(stdout, "dl: %+.16f\n", pos[i].dl);
      }

      return EXIT_SUCCESS;

}

int test_flatfile() {

  Accelerator accelerator;
  accelerator.energy = 0.0;
  accelerator.harmonic_number = 0;
  accelerator.cavity_on = false;
  accelerator.radiation_on = RadiationState::off;
  accelerator.vchamber_on = false;


  std::string fname = "/home/afonso/flatfile.txt";
  read_flat_file(fname, accelerator);
  std::cout << "Energy: " << accelerator.energy << " eV" << '\n';
  std::cout << "Harmonic number: " << accelerator.harmonic_number << '\n';
  std::cout << "Cavity on: " << accelerator.cavity_on << '\n';
  std::cout << "Radiation on: " << accelerator.radiation_on << '\n';
  std::cout << "Vacuum chamber on: " << accelerator.vchamber_on << '\n';
  accelerator.vchamber_on = false;
  fname = "/home/afonso/newflatfile.txt";
  write_flat_file(fname, accelerator);

  return EXIT_SUCCESS;

}

int test_calc_twiss() {

  Status::type status;

  Accelerator accelerator;
  std::string fname("sirius-v12.txt");
  status = read_flat_file(fname, accelerator);
  if (status != Status::success) {
    std::cerr << "could not open flat_file!" << std::endl;
    return status;
  }

  accelerator.cavity_on = true;
  accelerator.radiation_on = RadiationState::damping;
  accelerator.vchamber_on = false;

  Pos<double> fixed_point_guess;
  std::vector<Pos<double>> closed_orbit;
  status = track_findorbit6(accelerator, closed_orbit, fixed_point_guess);
  if (status != Status::success) {
    std::cerr << "could not find 6d closed orbit" << std::endl;
  }

  auto start = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();
  auto diff = end - start;

  std::vector<Twiss> twiss;
  Matrix m66;


  start = std::chrono::steady_clock::now();
  status = calc_twiss(accelerator, closed_orbit[0], m66, twiss);
  end = std::chrono::steady_clock::now(); diff = end - start;

  std::cout << "calc_twiss: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

  if (status != Status::success) {
    std::cerr << "could not calculate twiss" << std::endl;
  }


  for(unsigned int i=0; i<twiss.size(); ++i) {
    //std::cout << twiss[i] << std::endl;
  }


}

int test_matrix_inversion() {


  // without radiation
  // Matrix m3 = {
  //   {  8.4960e-01,   9.5487e+00,   4.8134e-02,  -1.0850e+00,   3.3269e-06,  -2.7194e-07 },
  //   { -3.0764e-02,   8.4960e-01,  -1.3468e-02,  -1.9919e-02,   6.5802e-07,  -1.2677e-08 },
  //   { -1.9919e-02,  -1.0850e+00,   4.7728e-01,   4.4408e+00,  -7.6608e-08,   1.3834e-08 },
  //   { -1.3468e-02,   4.8134e-02,  -1.7740e-01,   4.7728e-01,   1.2033e-07,  -1.3208e-09 },
  //   {  1.2031e-08,  -1.4563e-07,   1.1888e-09,   4.9343e-09,   9.9992e-01,  -1.0472e-02 },
  //   {  6.5896e-07,   3.3188e-06,   1.2042e-07,  -7.6371e-08,   7.7621e-02,   9.9927e-01 },
  // };

  // with radiation
  Matrix m1({
    {  8.5097e-01,   9.5435e+00,   4.8249e-02,  -1.0855e+00,   1.1619e-05,   3.4703e-07 },
    { -3.0791e-02,   8.4784e-01,  -1.3454e-02,  -1.9875e-02,   6.6804e-06,  -6.2775e-08 },
    { -2.0008e-02,  -1.0838e+00,   4.7620e-01,   4.4448e+00,  -1.2785e-06,  -8.2503e-09 },
    { -1.3476e-02,   4.7875e-02,  -1.7726e-01,   4.7778e-01,   5.8921e-07,  -6.6216e-09 },
    {  2.8625e-07,   1.7946e-06,   5.3438e-08,  -3.6462e-08,   9.9961e-01,  -1.0303e-02 },
    {  5.9759e-06,   5.4213e-05,   4.6994e-07,  -3.8983e-06,   7.7609e-02,   9.9928e-01 },
  });

  Matrix m2(m1); m2.inverse();
  Matrix m3; m3.multiplication(m1,m2);
  Matrix m4(6); m4.eye();
  Matrix m5; m5.linear_combination(1,m4,-1,m3);

  std::cout << "M:" << std::endl; m1._print();
  std::cout << "inv(M):" << std::endl; m2._print();
  std::cout << "M*inv(M)-1:" << std::endl; m5._print();


  // //matrix_multiplication(m3,m1,m2);
  // Status::type status = m3.inverse();
  // if (status != Status::success) {
  //   std::cerr << "newton did not converge!" << std::endl;
  // }
  //
  // std::cout << "inv_m:" << std::endl;
  // m3._print();

}

int test_new_write_flat_file() {

  Status::type status;

  Accelerator accelerator, accelerator2;
  std::string fname("sirius-v12-small.txt");
  status = read_flat_file(fname, accelerator, true);
  if (status != Status::success) {
    std::cerr << "could not open flat_file!" << std::endl;
    return status;
  }


  accelerator.cavity_on = true;
  accelerator.radiation_on = RadiationState::damping;
  accelerator.vchamber_on = false;

  std::string a;

  status = write_flat_file(a, accelerator, false);
  if (status != Status::success) {
    std::cout << "could not serialize accelerator!" << std::endl;
    return status;
  } else {
    std::cout << "size of accelerator         (#elements) : " << accelerator.lattice.size() << std::endl;
    std::cout << "size of serialized accelerator (#chars) : " << a.size() << std::endl;
  }

  status = read_flat_file(a, accelerator2, false);
  if (status != Status::success) {
    std::cout << "could not unserialize accelerator!" << std::endl;
    return status;
  } else {
    if (accelerator == accelerator2) {
      std::cout << "unserialized accelerator is identical to original" << std::endl;
      return Status::success;
    } else {
      std::cout << "unserialized accelerator is NOT identical to original !" << std::endl;
      return Status::not_implemented;
    }
  }

}

int cmd_tests(const std::vector<std::string>& args) {

  //test_printlattice(accelerator);
  //test_findm66(accelerator);
  //test_linepass(accelerator);
  //test_ringpass(accelerator);
  //test_linepass_tpsa(the_ring);
  //test_findorbit4();
  //test_findorbit6();
  //test_read_flat_file(accelerator);
  //test_cmd_dynap_xy();
  //test_cmd_dynap_ex();
  //test_cmd_dynap_ma();
  //test_cmd_dynap_ma();
  //test_cmd_dynap_pxa();
  //test_cmd_dynap_xyfmap();
  //test_cmd_dynap_exfmap();
  //test_cmd_track_linepass();
  //test_kicktable(accelerator);
  //test_simple_drift();
  //test_simple_quadrupole();
  //test_linepass2();
  test_calc_twiss();
  //test_matrix_inversion();
  //test_new_write_flat_file();

  return 0;

}
