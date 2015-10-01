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

#include "commands.h"
#include <trackcpp/trackcpp.h>
#include <ctime>

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
  Status::type status = track_linepass(accelerator, pos, new_pos, element_offset, lost_plane, true);
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
  track_linepass(accelerator, tpsa, new_tpsa, element_offset, lost_plane, false);
  for(unsigned int i=0; i<new_tpsa.size(); ++i) {
    //const Pos<Tpsa<6,1> >& c = new_particles[i];
    //std::cout << c.rx << std::endl;
    //fprintf(stdout, "%03i: %15s  %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E\n", i+1, the_ring[i].fam_name.c_str(), c.rx, c.px, c.ry, c.py, c.de);
  }

  return 0;

}


#include <cstdlib>
int test_ringpass(const Accelerator& accelerator) {


  Pos<> pos;
  pos.rx = 0.00100;
  pos.ry = 0.00010;

  std::vector<Pos<double> > new_pos;
  unsigned int element_offset = 0, lost_turn = 0;
  Plane::type lost_plane;

  clock_t begin, end;
  double time_spent;
  begin = clock();
  Status::type status = track_ringpass(accelerator, pos, new_pos, 5000, lost_turn, element_offset, lost_plane, true);
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
  read_flat_file("/home/fac_files/code/trackcpp/tests/si_v07_c05.txt", accelerator);
  accelerator.cavity_on = false;
  accelerator.radiation_on = false;
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
  read_flat_file("/home/fac_files/code/trackcpp/tests/si_v07_c05.txt", accelerator);
  accelerator.cavity_on = true;
  accelerator.radiation_on = true;
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

// int test_dynap_xy(const Accelerator& accelerator) {
//
//   std::vector<Pos<double> > closed_orbit;
//   unsigned int nr_turns = 5000;
//   Pos<double> p0(0,0,0,0,0,0);
//   unsigned int nrpts_x = 10;
//   double x_min = -0.015, x_max = +0.015;
//   unsigned int nrpts_y = 10;
//   double y_min = 0, y_max = +0.0035;
//   std::vector<DynApGridPoint> points;
//   dynap_xy(accelerator, closed_orbit, nr_turns, p0, nrpts_x, x_min, x_max, nrpts_y, y_min, y_max, true, points);
//
//   return 0;
// }

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

// int test_cmd_dynap_ma_threads() {
//
//   std::vector<std::string> args = {
//       "trackcpp",
//       "dynap_ma_threads",
//       "sirius-v10.txt",
//       "3e9",
//       "864",
//       "on",
//       "on",
//       "off",
//       "5000",
//       "30e-6",   // y0
//       "0.01",    // e0
//       "0.0005",  // tol_e
//       "0",       // s_min [m]
//       "30",      // s_max [m]
//       "4",       // nr_threads
//       "qaf",
//       "qad"};
//
//   return cmd_dynap_ma_threads(args);
//
// }

// int test_cmd_dynap_xy() {
//
//   std::vector<std::string> args = {
//       "trackcpp",
//       "dynap_xy",
//       "/home/ximenes/pytrack/sirius_v500_ac10_5_bare_with_ids_in.txt",
//       "3e9",
//       "864",
//       "on",
//       "on",
//       "on",
//       "0.0",
//       "5000",
//       "4",
//       "-0.015",
//       "+0.015",
//       "4",
//       "0.0",
//       "+0.0035"
//   };
//   return cmd_dynap_xy(args);
//
// }
//
// int test_cmd_dynap_ex() {
//
//   std::vector<std::string> args = {
//       "trackcpp",
//       "dynap_ex",
//       "/home/ximenes/pytrack/sirius_v500_ac10_5_bare_in.txt",
//       "3e9",
//       "864",
//       "on",
//       "on",
//       "on",
//       "1e-6",
//       "5000",
//       "2",
//       "-0.05",
//       "+0.05",
//       "2",
//       "-0.015",
//       "+0.015"
//   };
//   return cmd_dynap_xy(args);
//
// }

// int test_cmd_dynap_ma() {
//
//   std::vector<std::string> args = {
//       "trackcpp",
//       "dynap_ma",
//       "sirius-v10.txt",
//       "3e9",
//       "864",
//       "on",
//       "on",
//       "off",
//       "50",      // nr_turns
//       "30e-6",   // y0
//       "0.01",    // e0
//       "0.005",   // tol_e
//       "0",       // s_min [m]
//       "200",     // s_max [m]
//       "4",       // nr_threads
//       "qf2", "qfa", "sda"
//        };
//
//   return cmd_dynap_ma(args);
//
// }

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


// int test_cmd_dynap_pxa() {
//
//   std::vector<std::string> args = {
//       "trackcpp",
//       "dynap_pxa",
//       "sirius-v10.txt",
//       "3e9",
//       "864",
//       "on",
//       "on",
//       "off",
//       "50",      // nr_turns
//       "30e-6",   // y0
//       "1e-6",    // px
//       "1e-6",    // tol_e
//       "0",       // s_min [m]
//       "200",     // s_max [m]
//       "4",       // nr_threads
//       "qf2", "qfa", "sda"
//        };
//
//   return cmd_dynap_pxa(args);
//
// }

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
  add_kicktable("/home/fac_files/code/python/trackcpp/pytrack/id_kicktable.txt", accelerator.kicktables, ptrKicktable);
  add_kicktable("/home/fac_files/code/python/trackcpp/pytrack/id_kicktable2.txt", accelerator.kicktables, ptrKicktable);
  return 0;

}

int test_read_flat_file(Accelerator& accelerator) {

  read_flat_file("/home/ximenes/flatfile.txt", accelerator);
  std::cout << "nr_elements: " << accelerator.lattice.size() << std::endl;
  return 0;

}

int test_simple_drift() {

  Accelerator accelerator;

  accelerator.energy = 3e9; // [ev]
  accelerator.harmonic_number = 864;
  accelerator.radiation_on = false;
  accelerator.cavity_on = false;
  accelerator.vchamber_on = false;

  Element ds = Element::drift("ds", 1.0);

  accelerator.lattice.push_back(ds);

  Pos<double> pos(0.001,0.002,0.003,0.004,0.005,0.006);
  track_elementpass (ds, pos, accelerator);

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
  accelerator.radiation_on = false;
  accelerator.cavity_on = false;
  accelerator.vchamber_on = false;

  Element ds = Element::quadrupole("qs", 1.0, 2.0);

  accelerator.lattice.push_back(ds);

  Pos<double> pos(0.001,0.002,0.003,0.004,0.005,0.006);
  track_elementpass (ds, pos, accelerator);

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
  accelerator.radiation_on = false;
  accelerator.cavity_on = false;
  accelerator.vchamber_on = false;

  Pos<double> orig_pos;
  std::vector<Pos<double>> pos;
  unsigned int element_offset = 0;
  Plane::type lost_plane;
  bool trajectory = true;

  orig_pos.rx = 0.0001;
  orig_pos.px = 0.0001;

  Status::type status = track_linepass (accelerator,
      orig_pos,              // initial electron coordinates
      pos,     // vector with electron coordinates from tracking at every element.
      element_offset,  // index of starting element for tracking
      lost_plane,       // return plane in which particle was lost, if the case.
      trajectory);

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
  accelerator.radiation_on = false;
  accelerator.vchamber_on = false;

  read_flat_file("/home/afonso/flatfile.txt", accelerator);
  std::cout << "Energy: " << accelerator.energy << " eV" << '\n';
  std::cout << "Harmonic number: " << accelerator.harmonic_number << '\n';
  std::cout << "Cavity on: " << accelerator.cavity_on << '\n';
  std::cout << "Radiation on: " << accelerator.radiation_on << '\n';
  std::cout << "Vacuum chamber on: " << accelerator.vchamber_on << '\n';
  accelerator.vchamber_on = false;
  write_flat_file("/home/afonso/newflatfile.txt", accelerator);

  return EXIT_SUCCESS;

}

int test_calc_twiss() {

  Accelerator accelerator;
  read_flat_file("sirius-v10.txt", accelerator);
  accelerator.cavity_on = true;
  accelerator.radiation_on = true;
  accelerator.vchamber_on = false;

  std::cout << "Energy: " << accelerator.energy << " eV" << '\n';
  std::cout << "Harmonic number: " << accelerator.harmonic_number << '\n';
  std::cout << "Cavity on: " << accelerator.cavity_on << '\n';
  std::cout << "Radiation on: " << accelerator.radiation_on << '\n';
  std::cout << "Vacuum chamber on: " << accelerator.vchamber_on << '\n';

  Pos<double> fixed_point_guess;
  std::vector<Pos<double>> closed_orbit;
  Status::type status = track_findorbit6(accelerator, closed_orbit, fixed_point_guess);
  if (status != Status::success) {
    std::cerr << "could not find 6d closed orbit" << std::endl;
  }

  std::vector<Twiss> twiss;
  Matrix m66;
  status = calc_twiss(accelerator, closed_orbit[0], m66, twiss);
  if (status != Status::success) {
    std::cerr << "could not calculate twiss" << std::endl;
  }

  for(unsigned int i=0; i<twiss.size(); ++i) {
    std::cout << twiss[i].mux << std::endl;
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

  return 0;

}