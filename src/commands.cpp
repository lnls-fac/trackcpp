/// TRACKCPP - Particle tracking code
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
#include <trackcpp/output.h>
#include <trackcpp/dynap.h>
#include <trackcpp/flat_file.h>
#include <trackcpp/tracking.h>
#include <trackcpp/lattice.h>
#include <trackcpp/accelerator.h>
#include <trackcpp/elements.h>
#include <trackcpp/auxiliary.h>
#include <string>
#include <cstdlib>
#include <iostream>

int cmd_dynap_xy(const std::vector<std::string>& args) {

  if ((args.size() < 16) or (args.size() > 17)) {
    std::cerr << "dynap_xy: invalid number of arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "[cmd_dynap_xy]" << std::endl << std::endl;

  std::string  flat_filename(args[2]);
  double       ring_energy     = std::atof(args[3].c_str());
  unsigned int harmonic_number = std::atoi(args[4].c_str());
  std::string  cavity_state(args[5]);
  std::string  radiation_state(args[6]);
  std::string  vchamber_state(args[7]);
  double       de = std::atof(args[8].c_str());
  unsigned int nr_turns = std::atoi(args[9].c_str());
  unsigned int x_nrpts = std::atoi(args[10].c_str());
  double       x_min = std::atof(args[11].c_str());
  double       x_max = std::atof(args[12].c_str());
  unsigned int y_nrpts = std::atoi(args[13].c_str());
  double       y_min = std::atof(args[14].c_str());
  double       y_max = std::atof(args[15].c_str());
  unsigned int nr_threads = 1;
  if (args.size() == 17) {
    nr_threads = std::atoi(args[16].c_str());
  }

  print_header(stdout);
  std::cout << std::endl;
  std::cout << "flat_filename   : " << flat_filename << std::endl;
  std::cout << "energy[eV]      : " << ring_energy << std::endl;
  std::cout << "harmonic_number : " << harmonic_number << std::endl;
  std::cout << "cavity_state    : " << cavity_state << std::endl;
  std::cout << "radiation_state : " << radiation_state << std::endl;
  std::cout << "vchamber_state  : " << vchamber_state << std::endl;
  std::cout << "de              : " << de << std::endl;
  std::cout << "nr_turns        : " << nr_turns << std::endl;
  std::cout << "x_nrpts         : " << x_nrpts << std::endl;
  std::cout << "x_min[m]        : " << x_min << std::endl;
  std::cout << "x_max[m]        : " << x_max << std::endl;
  std::cout << "y_nrpts         : " << y_nrpts << std::endl;
  std::cout << "y_min[m]        : " << y_min << std::endl;
  std::cout << "y_max[m]        : " << y_max << std::endl;
  std::cout << "nr_threads      : " << nr_threads << std::endl;

  std::cout << std::endl;
  std::cout << get_timestamp() << " begin timestamp" << std::endl;

  Accelerator accelerator;

  // reads flat file
  Status::type status = read_flat_file(flat_filename, accelerator);
  if (status == Status::file_not_found) {
    std::cerr << "dynap_xy: flat file not found!" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << get_timestamp() << " input file with flat lattice read." << std::endl;

  // builds accelerator
  accelerator.energy = ring_energy;
  accelerator.harmonic_number = harmonic_number;
  accelerator.cavity_on = (cavity_state == "on");
  accelerator.radiation_on = (radiation_state == "on") ? RadiationState::damping : RadiationState::off;
  accelerator.vchamber_on = (vchamber_state == "on");

  // calcs dynamical aperture
  std::vector<Pos<double> > cod;
  Pos<double> p0(0,0,0,0,de,0);
  std::vector<DynApGridPoint> grid;
  dynap_xy(accelerator, cod, nr_turns, p0, x_nrpts, x_min, x_max, y_nrpts, y_min, y_max, true, grid, nr_threads);

  // generates output files
  std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
  status = print_closed_orbit(accelerator, cod);
  if (status == Status::file_not_opened) return status;
  std::cout << get_timestamp() << " saving dynap_xy grid to file" << std::endl;
  status = print_dynapgrid (accelerator, grid, "[dynap_xy]", "dynap_xy_out.txt");
  if (status == Status::file_not_opened) return status;

  std::cout << get_timestamp() << " end timestamp" << std::endl;
  return EXIT_SUCCESS;

}

int cmd_dynap_ex(const std::vector<std::string>& args) {

  if ((args.size() < 16) or (args.size() > 17)) {
    std::cerr << "dynap_ex: invalid number of arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "[cmd_dynap_ex]" << std::endl << std::endl;

  std::string  flat_filename(args[2]);
  double       ring_energy     = std::atof(args[3].c_str());
  unsigned int harmonic_number = std::atoi(args[4].c_str());
  std::string  cavity_state(args[5]);
  std::string  radiation_state(args[6]);
  std::string  vchamber_state(args[7]);
  double       y = std::atof(args[8].c_str());
  unsigned int nr_turns = std::atoi(args[9].c_str());
  unsigned int e_nrpts = std::atoi(args[10].c_str());
  double       e_min = std::atof(args[11].c_str());
  double       e_max = std::atof(args[12].c_str());
  unsigned int x_nrpts = std::atoi(args[13].c_str());
  double       x_min = std::atof(args[14].c_str());
  double       x_max = std::atof(args[15].c_str());
  unsigned int nr_threads = 1;
  if (args.size() == 17) {
    nr_threads = std::atoi(args[16].c_str());
  }

  print_header(stdout);
  std::cout << std::endl;
  std::cout << "flat_filename   : " << flat_filename << std::endl;
  std::cout << "energy[eV]      : " << ring_energy << std::endl;
  std::cout << "harmonic_number : " << harmonic_number << std::endl;
  std::cout << "cavity_state    : " << cavity_state << std::endl;
  std::cout << "radiation_state : " << radiation_state << std::endl;
  std::cout << "vchamber_state  : " << vchamber_state << std::endl;
  std::cout << "y[m]            : " << y << std::endl;
  std::cout << "nr_turns        : " << nr_turns << std::endl;
  std::cout << "e_nrpts         : " << e_nrpts << std::endl;
  std::cout << "e_min           : " << e_min << std::endl;
  std::cout << "e_max           : " << e_max << std::endl;
  std::cout << "x_nrpts         : " << x_nrpts << std::endl;
  std::cout << "x_min[m]        : " << x_min << std::endl;
  std::cout << "x_max[m]        : " << x_max << std::endl;
  std::cout << "nr_threads      : " << nr_threads << std::endl;

  std::cout << std::endl;
  std::cout << get_timestamp() << " begin timestamp" << std::endl;

  Accelerator accelerator;

  // reads flat file
  Status::type status = read_flat_file(flat_filename, accelerator);
  if (status == Status::file_not_found) {
    std::cerr << "dynap_ex: flat file not found!" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << get_timestamp() << " input file with flat lattice read." << std::endl;

  // builds accelerator
  accelerator.energy = ring_energy;
  accelerator.harmonic_number = harmonic_number;
  accelerator.cavity_on = (cavity_state == "on");
  accelerator.radiation_on = (radiation_state == "on") ? RadiationState::damping : RadiationState::off;
  accelerator.vchamber_on = (vchamber_state == "on");

  // calcs dynamical aperture
  std::vector<Pos<double> > cod;
  Pos<double> p0(0,0,y,0,0,0);
  std::vector<DynApGridPoint> grid;
  dynap_ex(accelerator, cod, nr_turns, p0, e_nrpts, e_min, e_max, x_nrpts, x_min, x_max, true, grid, nr_threads);

  // generates output files
  std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
  status = print_closed_orbit(accelerator, cod);
  if (status == Status::file_not_opened) return status;
  std::cout << get_timestamp() << " saving dynap_ex grid to file" << std::endl;
  status = print_dynapgrid (accelerator, grid, "[dynap_ex]", "dynap_ex_out.txt");
  if (status == Status::file_not_opened) return status;

  std::cout << get_timestamp() << " end timestamp" << std::endl;
  return EXIT_SUCCESS;

}

// int cmd_dynap_ma(const std::vector<std::string>& args) {
//
//   if (args.size() < 15) {
//     std::cerr << "dynap_ma: invalid number of arguments!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   std::cout << std::endl;
//   std::cout << "[cmd_dynap_ma]" << std::endl << std::endl;
//
//   std::string  flat_filename(args[2]);
//   double       ring_energy      = std::atof(args[3].c_str());
//   unsigned int harmonic_number = std::atoi(args[4].c_str());
//   std::string  cavity_state(args[5].c_str());
//   std::string  radiation_state(args[6].c_str());
//   std::string  vchamber_state(args[7].c_str());
//   unsigned int nr_turns = std::atoi(args[8].c_str());
//   double       y0    = std::atof(args[9].c_str());
//   double       e0    = std::atof(args[10].c_str());
//   double       e_tol = std::atof(args[11].c_str());
//   double       s_min = std::atof(args[12].c_str());
//   double       s_max = std::atof(args[13].c_str());
//
//   unsigned int nr_threads = 1;
//   std::vector<std::string> fam_names;
//   for(unsigned int i=14; i<args.size(); ++i) {
//     unsigned int nr = std::atoi(args[i].c_str());
//     if (nr > 0) {
//       nr_threads = nr;
//     } else fam_names.push_back(args[i]);
//   }
//
//
//   print_header(stdout);
//   std::cout << std::endl;
//   std::cout << "flat_filename   : " << flat_filename << std::endl;
//   std::cout << "energy[eV]      : " << ring_energy << std::endl;
//   std::cout << "harmonic_number : " << harmonic_number << std::endl;
//   std::cout << "cavity_state    : " << cavity_state << std::endl;
//   std::cout << "radiation_state : " << radiation_state << std::endl;
//   std::cout << "vchamber_state  : " << vchamber_state << std::endl;
//   std::cout << "nr_turns        : " << nr_turns << std::endl;
//   std::cout << "y0[m]           : " << y0 << std::endl;
//   std::cout << "e0              : " << e0 << std::endl;
//   std::cout << "e_tol           : " << e_tol << std::endl;
//   std::cout << "s_min[m]        : " << s_min << std::endl;
//   std::cout << "s_max[m]        : " << s_max << std::endl;
//   std::cout << "nr_threads      : " << nr_threads << std::endl;
//   std::cout << "fam_names       : ";
//   for(unsigned int i=0; i<fam_names.size(); ++i) std::cout << fam_names[i] << " "; std::cout << std::endl;
//
//   std::cout << std::endl;
//   std::cout << get_timestamp() << " begin timestamp" << std::endl;
//
//   Accelerator accelerator;
//
//   // reads flat file
//   Status::type status = read_flat_file(flat_filename, accelerator);
//   if (status == Status::file_not_found) {
//     std::cerr << "dynap_ma: flat file not found!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   // builds accelerator
//   accelerator.energy = ring_energy;
//   accelerator.harmonic_number = harmonic_number;
//   accelerator.cavity_on = (cavity_state == "on");
//   accelerator.radiation_on = (radiation_state == "on");
//   accelerator.vchamber_on = (vchamber_state == "on");
//
//   // calcs dynamical aperture
//   std::vector<Pos<double> > cod;
//   Pos<double> p0(0,0,y0,0,0,0);
//   std::vector<DynApGridPoint> grid;
//   dynap_ma(accelerator, cod, nr_turns, p0, e0, e_tol, s_min, s_max, fam_names, true, grid, nr_threads);
//
//   // generates output files
//   std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
//   status = print_closed_orbit(accelerator, cod);
//   if (status == Status::file_not_opened) return status;
//   std::cout << get_timestamp() << " saving dynap_ma grid to file" << std::endl;
//   status = print_dynapgrid (accelerator, grid, "[dynap_ma]", "dynap_ma_out.txt");
//   if (status == Status::file_not_opened) return status;
//
//   std::cout << get_timestamp() << " end timestamp" << std::endl;
//   return EXIT_SUCCESS;
//
// }

int cmd_dynap_acceptance(const std::vector<std::string>& args) {

  std::string calc_type_str, p_init_str, p_delta_str;
  if (args[1] == "dynap_ma") {
    calc_type_str = "dynap_ma "; p_init_str = "e_init[rad] "; p_delta_str = "e_delta[rad] ";
  } else if (args[1] == "dynap_pxa") {
    calc_type_str = "dynap_pxa"; p_init_str = "px_init[rad]"; p_delta_str = "px_delta[rad]";
  } else if (args[1] == "dynap_pya") {
    calc_type_str = "dynap_pya"; p_init_str = "py_init[rad]"; p_delta_str = "py_delta[rad]";
  }

  if (args.size() < 17) {
    std::cerr << args[1] << ": invalid number of arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "[" << args[1] << "]" << std::endl << std::endl;

  std::string  flat_filename(args[2]);
  double       ring_energy      = std::atof(args[3].c_str());
  unsigned int harmonic_number = std::atoi(args[4].c_str());
  std::string  cavity_state(args[5].c_str());
  std::string  radiation_state(args[6].c_str());
  std::string  vchamber_state(args[7].c_str());
  unsigned int nr_turns = std::atoi(args[8].c_str());
  double       y0       = std::atof(args[9].c_str());
  double       p_init   = std::atof(args[10].c_str());
  double       p_delta  = std::atof(args[11].c_str());
  unsigned int nr_steps_back = std::atoi(args[12].c_str());
  double       rescale = std::atof(args[13].c_str());
  unsigned int nr_iterations = std::atoi(args[14].c_str());
  double       s_min    = std::atof(args[15].c_str());
  double       s_max    = std::atof(args[16].c_str());
  unsigned int nr_threads = 1;
  std::vector<std::string> fam_names;
  for(unsigned int i=17; i<args.size(); ++i) {
    unsigned int nr = std::atoi(args[i].c_str());
    if (nr > 0) {
      nr_threads = nr;
    } else fam_names.push_back(args[i]);
  }

  print_header(stdout);
  std::cout << std::endl;
  std::cout <<   "flat_filename   : " << flat_filename << std::endl;
  std::cout <<   "energy[eV]      : " << ring_energy << std::endl;
  std::cout <<   "harmonic_number : " << harmonic_number << std::endl;
  std::cout <<   "cavity_state    : " << cavity_state << std::endl;
  std::cout <<   "radiation_state : " << radiation_state << std::endl;
  std::cout <<   "vchamber_state  : " << vchamber_state << std::endl;
  std::cout <<   "nr_turns        : " << nr_turns << std::endl;
  std::cout <<   "y0[m]           : " << y0 << std::endl;
  std::cout << p_init_str << "    : " << p_init << std::endl;
  std::cout << p_delta_str << "   : " << p_delta << std::endl;
  std::cout <<   "nr_steps_back   : " << nr_steps_back << std::endl;
  std::cout <<   "rescale         : " << rescale << std::endl;
  std::cout <<   "nr_iterations   : " << nr_iterations << std::endl;
  std::cout <<   "s_min[m]        : " << s_min << std::endl;
  std::cout <<   "s_max[m]        : " << s_max << std::endl;
  std::cout <<   "nr_threads      : " << nr_threads << std::endl;
  std::cout <<   "fam_names       : ";
  for(unsigned int i=0; i<fam_names.size(); ++i) std::cout << fam_names[i] << " "; std::cout << std::endl;

  std::cout << std::endl;
  std::cout << get_timestamp() << " begin timestamp" << std::endl;

  Accelerator accelerator;

  // reads flat file
  Status::type status = read_flat_file(flat_filename, accelerator);
  if (status == Status::file_not_found) {
    std::cerr << args[1] << ": flat file not found!" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << get_timestamp() << " input file with flat lattice read." << std::endl;

  // builds accelerator
  accelerator.energy = ring_energy;
  accelerator.harmonic_number = harmonic_number;
  accelerator.cavity_on = (cavity_state == "on");
  accelerator.radiation_on = (radiation_state == "on") ? RadiationState::damping : RadiationState::off;
  accelerator.vchamber_on = (vchamber_state == "on");

  // calcs dynamical aperture
  std::vector<Pos<double> > cod;
  Pos<double> p0(0,0,y0,0,0,0);
  std::vector<DynApGridPoint> grid;
  //dynap_pxa(accelerator, cod, nr_turns, p0, p_init, p_delta, nr_steps_back, rescale, nr_iterations, s_min, s_max, fam_names, true, grid, nr_threads);
  dynap_acceptance(args[1], accelerator, cod, nr_turns, p0, p_init, p_delta, nr_steps_back, rescale, nr_iterations, s_min, s_max, fam_names, true, grid, nr_threads);

  // generates output files
  std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
  status = print_closed_orbit(accelerator, cod);
  if (status == Status::file_not_opened) return status;
  std::cout << get_timestamp() << " saving " << args[1] << " grid to file" << std::endl;
  status = print_dynapgrid (accelerator, grid, "["+args[1]+"]", args[1]+"_out.txt");
  if (status == Status::file_not_opened) return status;

  std::cout << get_timestamp() << " end timestamp" << std::endl;
  return EXIT_SUCCESS;

}



// int cmd_dynap_ma(const std::vector<std::string>& args) {
//
//   if (args.size() < 17) {
//     std::cerr << "dynap_ma: invalid number of arguments!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   std::cout << std::endl;
//   std::cout << "[cmd_dynap_ma]" << std::endl << std::endl;
//
//   std::string  flat_filename(args[2]);
//   double       ring_energy      = std::atof(args[3].c_str());
//   unsigned int harmonic_number = std::atoi(args[4].c_str());
//   std::string  cavity_state(args[5].c_str());
//   std::string  radiation_state(args[6].c_str());
//   std::string  vchamber_state(args[7].c_str());
//   unsigned int nr_turns = std::atoi(args[8].c_str());
//   double       y0      = std::atof(args[9].c_str());
//   double       e_init  = std::atof(args[10].c_str());
//   double       e_delta = std::atof(args[11].c_str());
//   unsigned int nr_steps_back = std::atoi(args[12].c_str());
//   double       rescale = std::atof(args[13].c_str());
//   unsigned int nr_iterations = std::atoi(args[14].c_str());
//   double       s_min = std::atof(args[15].c_str());
//   double       s_max = std::atof(args[16].c_str());
//   unsigned int nr_threads = 1;
//   std::vector<std::string> fam_names;
//   for(unsigned int i=17; i<args.size(); ++i) {
//     unsigned int nr = std::atoi(args[i].c_str());
//     if (nr > 0) {
//       nr_threads = nr;
//     } else fam_names.push_back(args[i]);
//   }
//
//   print_header(stdout);
//   std::cout << std::endl;
//   std::cout << "flat_filename   : " << flat_filename << std::endl;
//   std::cout << "energy[eV]      : " << ring_energy << std::endl;
//   std::cout << "harmonic_number : " << harmonic_number << std::endl;
//   std::cout << "cavity_state    : " << cavity_state << std::endl;
//   std::cout << "radiation_state : " << radiation_state << std::endl;
//   std::cout << "vchamber_state  : " << vchamber_state << std::endl;
//   std::cout << "nr_turns        : " << nr_turns << std::endl;
//   std::cout << "y0[m]           : " << y0 << std::endl;
//   std::cout << "e_init          : " << e_init << std::endl;
//   std::cout << "e_delta         : " << e_delta << std::endl;
//   std::cout << "nr_steps_back   : " << nr_steps_back << std::endl;
//   std::cout << "rescale         : " << rescale << std::endl;
//   std::cout << "nr_iterations   : " << nr_iterations << std::endl;
//   std::cout << "s_min[m]        : " << s_min << std::endl;
//   std::cout << "s_max[m]        : " << s_max << std::endl;
//   std::cout << "nr_threads      : " << nr_threads << std::endl;
//   std::cout << "fam_names       : ";
//   for(unsigned int i=0; i<fam_names.size(); ++i) std::cout << fam_names[i] << " "; std::cout << std::endl;
//
//   std::cout << std::endl;
//   std::cout << get_timestamp() << " begin timestamp" << std::endl;
//
//   Accelerator accelerator;
//
//   // reads flat file
//   Status::type status = read_flat_file(flat_filename, accelerator);
//   if (status == Status::file_not_found) {
//     std::cerr << "dynap_ma: flat file not found!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   // builds accelerator
//   accelerator.energy = ring_energy;
//   accelerator.harmonic_number = harmonic_number;
//   accelerator.cavity_on = (cavity_state == "on");
//   accelerator.radiation_on = (radiation_state == "on");
//   accelerator.vchamber_on = (vchamber_state == "on");
//
//   // calcs dynamical aperture
//   std::vector<Pos<double> > cod;
//   Pos<double> p0(0,0,y0,0,0,0);
//   std::vector<DynApGridPoint> grid;
//   dynap_ma(accelerator, cod, nr_turns, p0, e_init, e_delta, nr_steps_back, rescale, nr_iterations, s_min, s_max, fam_names, true, grid, nr_threads);
//
//   // generates output files
//   std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
//   status = print_closed_orbit(accelerator, cod);
//   if (status == Status::file_not_opened) return status;
//   std::cout << get_timestamp() << " saving dynap_ma grid to file" << std::endl;
//   status = print_dynapgrid (accelerator, grid, "[dynap_ma]", "dynap_ma2_out.txt");
//   if (status == Status::file_not_opened) return status;
//
//   std::cout << get_timestamp() << " end timestamp" << std::endl;
//   return EXIT_SUCCESS;
//
// }

// int cmd_dynap_pxa(const std::vector<std::string>& args) {
//
//   if (args.size() < 17) {
//     std::cerr << "dynap_pxa: invalid number of arguments!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   std::cout << std::endl;
//   std::cout << "[cmd_dynap_pxa]" << std::endl << std::endl;
//
//   std::string  flat_filename(args[2]);
//   double       ring_energy      = std::atof(args[3].c_str());
//   unsigned int harmonic_number = std::atoi(args[4].c_str());
//   std::string  cavity_state(args[5].c_str());
//   std::string  radiation_state(args[6].c_str());
//   std::string  vchamber_state(args[7].c_str());
//   unsigned int nr_turns = std::atoi(args[8].c_str());
//   double       y0       = std::atof(args[9].c_str());
//   double       px_init  = std::atof(args[10].c_str());
//   double       px_delta = std::atof(args[11].c_str());
//   unsigned int nr_steps_back = std::atoi(args[12].c_str());
//   double       rescale = std::atof(args[13].c_str());
//   unsigned int nr_iterations = std::atoi(args[14].c_str());
//   double       s_min    = std::atof(args[15].c_str());
//   double       s_max    = std::atof(args[16].c_str());
//   unsigned int nr_threads = 1;
//   std::vector<std::string> fam_names;
//   for(unsigned int i=17; i<args.size(); ++i) {
//     unsigned int nr = std::atoi(args[i].c_str());
//     if (nr > 0) {
//       nr_threads = nr;
//     } else fam_names.push_back(args[i]);
//   }
//
//   print_header(stdout);
//   std::cout << std::endl;
//   std::cout << "flat_filename   : " << flat_filename << std::endl;
//   std::cout << "energy[eV]      : " << ring_energy << std::endl;
//   std::cout << "harmonic_number : " << harmonic_number << std::endl;
//   std::cout << "cavity_state    : " << cavity_state << std::endl;
//   std::cout << "radiation_state : " << radiation_state << std::endl;
//   std::cout << "vchamber_state  : " << vchamber_state << std::endl;
//   std::cout << "nr_turns        : " << nr_turns << std::endl;
//   std::cout << "y0[m]           : " << y0 << std::endl;
//   std::cout << "px_init[rad]    : " << px_init << std::endl;
//   std::cout << "px_delta[rad]   : " << px_delta << std::endl;
//   std::cout << "nr_steps_back   : " << nr_steps_back << std::endl;
//   std::cout << "rescale         : " << rescale << std::endl;
//   std::cout << "nr_iterations   : " << nr_iterations << std::endl;
//   std::cout << "s_min[m]        : " << s_min << std::endl;
//   std::cout << "s_max[m]        : " << s_max << std::endl;
//   std::cout << "nr_threads      : " << nr_threads << std::endl;
//   std::cout << "fam_names       : ";
//   for(unsigned int i=0; i<fam_names.size(); ++i) std::cout << fam_names[i] << " "; std::cout << std::endl;
//
//   std::cout << std::endl;
//   std::cout << get_timestamp() << " begin timestamp" << std::endl;
//
//   Accelerator accelerator;
//
//   // reads flat file
//   Status::type status = read_flat_file(flat_filename, accelerator);
//   if (status == Status::file_not_found) {
//     std::cerr << "dynap_pxa: flat file not found!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   // builds accelerator
//   accelerator.energy = ring_energy;
//   accelerator.harmonic_number = harmonic_number;
//   accelerator.cavity_on = (cavity_state == "on");
//   accelerator.radiation_on = (radiation_state == "on");
//   accelerator.vchamber_on = (vchamber_state == "on");
//
//   // calcs dynamical aperture
//   std::vector<Pos<double> > cod;
//   Pos<double> p0(0,0,y0,0,0,0);
//   std::vector<DynApGridPoint> grid;
//   dynap_pxa(accelerator, cod, nr_turns, p0, px_init, px_delta, nr_steps_back, rescale, nr_iterations, s_min, s_max, fam_names, true, grid, nr_threads);
//
//   // generates output files
//   std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
//   status = print_closed_orbit(accelerator, cod);
//   if (status == Status::file_not_opened) return status;
//   std::cout << get_timestamp() << " saving dynap_pxa grid to file" << std::endl;
//   status = print_dynapgrid (accelerator, grid, "[dynap_pxa]", "dynap_pxa_out.txt");
//   if (status == Status::file_not_opened) return status;
//
//   std::cout << get_timestamp() << " end timestamp" << std::endl;
//   return EXIT_SUCCESS;
//
// }

// int cmd_dynap_pxa_old(const std::vector<std::string>& args) {
//
//   if (args.size() < 15) {
//     std::cerr << "dynap_pxa: invalid number of arguments!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   std::cout << std::endl;
//   std::cout << "[cmd_dynap_pxa]" << std::endl << std::endl;
//
//   std::string  flat_filename(args[2]);
//   double       ring_energy      = std::atof(args[3].c_str());
//   unsigned int harmonic_number = std::atoi(args[4].c_str());
//   std::string  cavity_state(args[5].c_str());
//   std::string  radiation_state(args[6].c_str());
//   std::string  vchamber_state(args[7].c_str());
//   unsigned int nr_turns = std::atoi(args[8].c_str());
//   double       y0     = std::atof(args[9].c_str());
//   double       px0    = std::atof(args[10].c_str());
//   double       px_tol = std::atof(args[11].c_str());
//   double       s_min  = std::atof(args[12].c_str());
//   double       s_max  = std::atof(args[13].c_str());
//   unsigned int nr_threads = 1;
//   std::vector<std::string> fam_names;
//   for(unsigned int i=14; i<args.size(); ++i) {
//     unsigned int nr = std::atoi(args[i].c_str());
//     if (nr > 0) {
//       nr_threads = nr;
//     } else fam_names.push_back(args[i]);
//   }
//
//
//   print_header(stdout);
//   std::cout << std::endl;
//   std::cout << "flat_filename   : " << flat_filename << std::endl;
//   std::cout << "energy[eV]      : " << ring_energy << std::endl;
//   std::cout << "harmonic_number : " << harmonic_number << std::endl;
//   std::cout << "cavity_state    : " << cavity_state << std::endl;
//   std::cout << "radiation_state : " << radiation_state << std::endl;
//   std::cout << "vchamber_state  : " << vchamber_state << std::endl;
//   std::cout << "nr_turns        : " << nr_turns << std::endl;
//   std::cout << "y0[m]           : " << y0 << std::endl;
//   std::cout << "px0[rad]        : " << px0 << std::endl;
//   std::cout << "px_tol          : " << px_tol << std::endl;
//   std::cout << "s_min[m]        : " << s_min << std::endl;
//   std::cout << "s_max[m]        : " << s_max << std::endl;
//   std::cout << "nr_threads      : " << nr_threads << std::endl;
//   std::cout << "fam_names       : ";
//   for(unsigned int i=0; i<fam_names.size(); ++i) std::cout << fam_names[i] << " "; std::cout << std::endl;
//
//   std::cout << std::endl;
//   std::cout << get_timestamp() << " begin timestamp" << std::endl;
//
//   Accelerator accelerator;
//
//   // reads flat file
//   Status::type status = read_flat_file(flat_filename, accelerator);
//   if (status == Status::file_not_found) {
//     std::cerr << "dynap_pxa: flat file not found!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   // builds accelerator
//   accelerator.energy = ring_energy;
//   accelerator.harmonic_number = harmonic_number;
//   accelerator.cavity_on = (cavity_state == "on");
//   accelerator.radiation_on = (radiation_state == "on");
//   accelerator.vchamber_on = (vchamber_state == "on");
//
//   // calcs dynamical aperture
//   std::vector<Pos<double> > cod;
//   Pos<double> p0(0,0,y0,0,0,0);
//   std::vector<DynApGridPoint> grid;
//   dynap_pxa(accelerator, cod, nr_turns, p0, px0, px_tol, s_min, s_max, fam_names, true, grid, nr_threads);
//
//   // generates output files
//   std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
//   status = print_closed_orbit(accelerator, cod);
//   if (status == Status::file_not_opened) return status;
//   std::cout << get_timestamp() << " saving dynap_pxa grid to file" << std::endl;
//   status = print_dynapgrid (accelerator, grid, "[dynap_pxa]", "dynap_pxa_out.txt");
//   if (status == Status::file_not_opened) return status;
//
//   std::cout << get_timestamp() << " end timestamp" << std::endl;
//   return EXIT_SUCCESS;
//
// }

// int cmd_dynap_pya(const std::vector<std::string>& args) {
//
//   if (args.size() < 15) {
//     std::cerr << "dynap_pya: invalid number of arguments!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   std::cout << std::endl;
//   std::cout << "[cmd_dynap_pya]" << std::endl << std::endl;
//
//   std::string  flat_filename(args[2]);
//   double       ring_energy      = std::atof(args[3].c_str());
//   unsigned int harmonic_number = std::atoi(args[4].c_str());
//   std::string  cavity_state(args[5].c_str());
//   std::string  radiation_state(args[6].c_str());
//   std::string  vchamber_state(args[7].c_str());
//   unsigned int nr_turns = std::atoi(args[8].c_str());
//   double       y0     = std::atof(args[9].c_str());
//   double       py0    = std::atof(args[10].c_str());
//   double       py_tol = std::atof(args[11].c_str());
//   double       s_min  = std::atof(args[12].c_str());
//   double       s_max  = std::atof(args[13].c_str());
//   unsigned int nr_threads = 1;
//   std::vector<std::string> fam_names;
//   for(unsigned int i=14; i<args.size(); ++i) {
//     unsigned int nr = std::atoi(args[i].c_str());
//     if (nr > 0) {
//       nr_threads = nr;
//     } else fam_names.push_back(args[i]);
//   }
//
//
//   print_header(stdout);
//   std::cout << std::endl;
//   std::cout << "flat_filename   : " << flat_filename << std::endl;
//   std::cout << "energy[eV]      : " << ring_energy << std::endl;
//   std::cout << "harmonic_number : " << harmonic_number << std::endl;
//   std::cout << "cavity_state    : " << cavity_state << std::endl;
//   std::cout << "radiation_state : " << radiation_state << std::endl;
//   std::cout << "vchamber_state  : " << vchamber_state << std::endl;
//   std::cout << "nr_turns        : " << nr_turns << std::endl;
//   std::cout << "y0[m]           : " << y0 << std::endl;
//   std::cout << "py0[rad]        : " << py0 << std::endl;
//   std::cout << "py_tol          : " << py_tol << std::endl;
//   std::cout << "s_min[m]        : " << s_min << std::endl;
//   std::cout << "s_max[m]        : " << s_max << std::endl;
//   std::cout << "nr_threads      : " << nr_threads << std::endl;
//   std::cout << "fam_names       : ";
//   for(unsigned int i=0; i<fam_names.size(); ++i) std::cout << fam_names[i] << " "; std::cout << std::endl;
//
//   std::cout << std::endl;
//   std::cout << get_timestamp() << " begin timestamp" << std::endl;
//
//   Accelerator accelerator;
//
//   // reads flat file
//   Status::type status = read_flat_file(flat_filename, accelerator);
//   if (status == Status::file_not_found) {
//     std::cerr << "dynap_pya: flat file not found!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   // builds accelerator
//   accelerator.energy = ring_energy;
//   accelerator.harmonic_number = harmonic_number;
//   accelerator.cavity_on = (cavity_state == "on");
//   accelerator.radiation_on = (radiation_state == "on");
//   accelerator.vchamber_on = (vchamber_state == "on");
//
//   // calcs dynamical aperture
//   std::vector<Pos<double> > cod;
//   Pos<double> p0(0,0,y0,0,0,0);
//   std::vector<DynApGridPoint> grid;
//   dynap_pya(accelerator, cod, nr_turns, p0, py0, py_tol, s_min, s_max, fam_names, true, grid, nr_threads);
//
//   // generates output files
//   std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
//   status = print_closed_orbit(accelerator, cod);
//   if (status == Status::file_not_opened) return status;
//   std::cout << get_timestamp() << " saving dynap_pya grid to file" << std::endl;
//   status = print_dynapgrid (accelerator, grid, "[dynap_pya]", "dynap_pya_out.txt");
//   if (status == Status::file_not_opened) return status;
//
//   std::cout << get_timestamp() << " end timestamp" << std::endl;
//   return EXIT_SUCCESS;
//
// }

// int cmd_dynap_pya(const std::vector<std::string>& args) {
//
//   if (args.size() < 17) {
//     std::cerr << "dynap_pya: invalid number of arguments!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   std::cout << std::endl;
//   std::cout << "[cmd_dynap_pya]" << std::endl << std::endl;
//
//   std::string  flat_filename(args[2]);
//   double       ring_energy      = std::atof(args[3].c_str());
//   unsigned int harmonic_number = std::atoi(args[4].c_str());
//   std::string  cavity_state(args[5].c_str());
//   std::string  radiation_state(args[6].c_str());
//   std::string  vchamber_state(args[7].c_str());
//   unsigned int nr_turns = std::atoi(args[8].c_str());
//   double       y0       = std::atof(args[9].c_str());
//   double       py_init  = std::atof(args[10].c_str());
//   double       py_delta = std::atof(args[11].c_str());
//   unsigned int nr_steps_back = std::atoi(args[12].c_str());
//   double       rescale = std::atof(args[13].c_str());
//   unsigned int nr_iterations = std::atoi(args[14].c_str());
//   double       s_min    = std::atof(args[15].c_str());
//   double       s_max    = std::atof(args[16].c_str());
//   unsigned int nr_threads = 1;
//   std::vector<std::string> fam_names;
//   for(unsigned int i=17; i<args.size(); ++i) {
//     unsigned int nr = std::atoi(args[i].c_str());
//     if (nr > 0) {
//       nr_threads = nr;
//     } else fam_names.push_back(args[i]);
//   }
//
//   print_header(stdout);
//   std::cout << std::endl;
//   std::cout << "flat_filename   : " << flat_filename << std::endl;
//   std::cout << "energy[eV]      : " << ring_energy << std::endl;
//   std::cout << "harmonic_number : " << harmonic_number << std::endl;
//   std::cout << "cavity_state    : " << cavity_state << std::endl;
//   std::cout << "radiation_state : " << radiation_state << std::endl;
//   std::cout << "vchamber_state  : " << vchamber_state << std::endl;
//   std::cout << "nr_turns        : " << nr_turns << std::endl;
//   std::cout << "y0[m]           : " << y0 << std::endl;
//   std::cout << "py_init[rad]    : " << py_init << std::endl;
//   std::cout << "py_delta[rad]   : " << py_delta << std::endl;
//   std::cout << "nr_steps_back   : " << nr_steps_back << std::endl;
//   std::cout << "rescale         : " << rescale << std::endl;
//   std::cout << "nr_iterations   : " << nr_iterations << std::endl;
//   std::cout << "s_min[m]        : " << s_min << std::endl;
//   std::cout << "s_max[m]        : " << s_max << std::endl;
//   std::cout << "nr_threads      : " << nr_threads << std::endl;
//   std::cout << "fam_names       : ";
//   for(unsigned int i=0; i<fam_names.size(); ++i) std::cout << fam_names[i] << " "; std::cout << std::endl;
//
//   std::cout << std::endl;
//   std::cout << get_timestamp() << " begin timestamp" << std::endl;
//
//   Accelerator accelerator;
//
//   // reads flat file
//   Status::type status = read_flat_file(flat_filename, accelerator);
//   if (status == Status::file_not_found) {
//     std::cerr << "dynap_pya: flat file not found!" << std::endl;
//     return EXIT_FAILURE;
//   }
//
//   // builds accelerator
//   accelerator.energy = ring_energy;
//   accelerator.harmonic_number = harmonic_number;
//   accelerator.cavity_on = (cavity_state == "on");
//   accelerator.radiation_on = (radiation_state == "on");
//   accelerator.vchamber_on = (vchamber_state == "on");
//
//   // calcs dynamical aperture
//   std::vector<Pos<double> > cod;
//   Pos<double> p0(0,0,y0,0,0,0);
//   std::vector<DynApGridPoint> grid;
//   dynap_pya(accelerator, cod, nr_turns, p0, py_init, py_delta, nr_steps_back, rescale, nr_iterations, s_min, s_max, fam_names, true, grid, nr_threads);
//
//   // generates output files
//   std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
//   status = print_closed_orbit(accelerator, cod);
//   if (status == Status::file_not_opened) return status;
//   std::cout << get_timestamp() << " saving dynap_pya grid to file" << std::endl;
//   status = print_dynapgrid (accelerator, grid, "[dynap_pya]", "dynap_pya_out.txt");
//   if (status == Status::file_not_opened) return status;
//
//   std::cout << get_timestamp() << " end timestamp" << std::endl;
//   return EXIT_SUCCESS;
//
// }

int cmd_dynap_xyfmap(const std::vector<std::string>& args) {

  if (args.size() != 17) {
    std::cerr << "dynap_xyfmap: invalid number of arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "[cmd_dynap_xyfmap]" << std::endl << std::endl;

  std::string  flat_filename(args[2]);
  double       ring_energy     = std::atof(args[3].c_str());
  unsigned int harmonic_number = std::atoi(args[4].c_str());
  std::string  cavity_state(args[5]);
  std::string  radiation_state(args[6]);
  std::string  vchamber_state(args[7]);
  double       de = std::atof(args[8].c_str());
  unsigned int nr_turns = std::atoi(args[9].c_str());
  unsigned int x_nrpts = std::atoi(args[10].c_str());
  double       x_min = std::atof(args[11].c_str());
  double       x_max = std::atof(args[12].c_str());
  unsigned int y_nrpts = std::atoi(args[13].c_str());
  double       y_min = std::atof(args[14].c_str());
  double       y_max = std::atof(args[15].c_str());
  unsigned int nr_threads = std::atoi(args[16].c_str());

  print_header(stdout);
  std::cout << std::endl;
  std::cout << "flat_filename   : " << flat_filename << std::endl;
  std::cout << "energy[eV]      : " << ring_energy << std::endl;
  std::cout << "harmonic_number : " << harmonic_number << std::endl;
  std::cout << "cavity_state    : " << cavity_state << std::endl;
  std::cout << "radiation_state : " << radiation_state << std::endl;
  std::cout << "vchamber_state  : " << vchamber_state << std::endl;
  std::cout << "de              : " << de << std::endl;
  std::cout << "nr_turns        : " << nr_turns << std::endl;
  std::cout << "x_nrpts         : " << x_nrpts << std::endl;
  std::cout << "x_min[m]        : " << x_min << std::endl;
  std::cout << "x_max[m]        : " << x_max << std::endl;
  std::cout << "y_nrpts         : " << y_nrpts << std::endl;
  std::cout << "y_min[m]        : " << y_min << std::endl;
  std::cout << "y_max[m]        : " << y_max << std::endl;
  std::cout << "nr_threads      : " << nr_threads << std::endl;

  std::cout << std::endl;
  std::cout << get_timestamp() << " begin timestamp" << std::endl;

  Accelerator accelerator;

  // reads flat file
  Status::type status = read_flat_file(flat_filename, accelerator);
  if (status == Status::file_not_found) {
    std::cerr << "dynap_xyfmap: flat file not found!" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << get_timestamp() << " input file with flat lattice read." << std::endl;

  // builds accelerator
  accelerator.energy = ring_energy;
  accelerator.harmonic_number = harmonic_number;
  accelerator.cavity_on = (cavity_state == "on");
  accelerator.radiation_on = (radiation_state == "on") ? RadiationState::damping : RadiationState::off;
  accelerator.vchamber_on = (vchamber_state == "on");

  // calcs dynamical aperture
  std::vector<Pos<double> > cod;
  Pos<double> p0(0,0,0,0,de,0);
  std::vector<DynApGridPoint> grid;
  dynap_xyfmap(accelerator, cod, nr_turns, p0, x_nrpts, x_min, x_max, y_nrpts, y_min, y_max, true, grid, nr_threads);

  // generates output files
  std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
  status = print_closed_orbit(accelerator, cod);
  if (status == Status::file_not_opened) return status;
  std::cout << get_timestamp() << " saving dynap_xyfmap grid to file" << std::endl;
  status = print_dynapgrid (accelerator, grid, "[dynap_fmap]", "dynap_xyfmap_out.txt", true);
  if (status == Status::file_not_opened) return status;

  std::cout << get_timestamp() << " end timestamp" << std::endl;
  return EXIT_SUCCESS;

}

int cmd_dynap_exfmap(const std::vector<std::string>& args) {

  if (args.size() != 17) {
    std::cerr << "dynap_exfmap: invalid number of arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "[cmd_dynap_exfmap]" << std::endl << std::endl;

  std::string  flat_filename(args[2]);
  double       ring_energy     = std::atof(args[3].c_str());
  unsigned int harmonic_number = std::atoi(args[4].c_str());
  std::string  cavity_state(args[5]);
  std::string  radiation_state(args[6]);
  std::string  vchamber_state(args[7]);
  double       y = std::atof(args[8].c_str());
  unsigned int nr_turns = std::atoi(args[9].c_str());
  unsigned int e_nrpts = std::atoi(args[10].c_str());
  double       e_min = std::atof(args[11].c_str());
  double       e_max = std::atof(args[12].c_str());
  unsigned int x_nrpts = std::atoi(args[13].c_str());
  double       x_min = std::atof(args[14].c_str());
  double       x_max = std::atof(args[15].c_str());
  unsigned int nr_threads = std::atoi(args[16].c_str());

  print_header(stdout);
  std::cout << std::endl;
  std::cout << "flat_filename   : " << flat_filename << std::endl;
  std::cout << "energy[eV]      : " << ring_energy << std::endl;
  std::cout << "harmonic_number : " << harmonic_number << std::endl;
  std::cout << "cavity_state    : " << cavity_state << std::endl;
  std::cout << "radiation_state : " << radiation_state << std::endl;
  std::cout << "vchamber_state  : " << vchamber_state << std::endl;
  std::cout << "y               : " << y << std::endl;
  std::cout << "nr_turns        : " << nr_turns << std::endl;
  std::cout << "e_nrpts         : " << e_nrpts << std::endl;
  std::cout << "e_min[m]        : " << e_min << std::endl;
  std::cout << "e_max[m]        : " << e_max << std::endl;
  std::cout << "x_nrpts         : " << x_nrpts << std::endl;
  std::cout << "x_min[m]        : " << x_min << std::endl;
  std::cout << "x_max[m]        : " << x_max << std::endl;
  std::cout << "nr_threads      : " << nr_threads << std::endl;

  std::cout << std::endl;
  std::cout << get_timestamp() << " begin timestamp" << std::endl;

  Accelerator accelerator;

  // reads flat file
  Status::type status = read_flat_file(flat_filename, accelerator);
  if (status == Status::file_not_found) {
    std::cerr << "dynap_exfmap: flat file not found!" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << get_timestamp() << " input file with flat lattice read." << std::endl;

  // builds accelerator
  accelerator.energy = ring_energy;
  accelerator.harmonic_number = harmonic_number;
  accelerator.cavity_on = (cavity_state == "on");
  accelerator.radiation_on = (radiation_state == "on") ? RadiationState::damping : RadiationState::off;
  accelerator.vchamber_on = (vchamber_state == "on");

  // calcs dynamical aperture
  std::vector<Pos<double> > cod;
  Pos<double> p0(0,0,y,0,0,0);
  std::vector<DynApGridPoint> grid;
  dynap_exfmap(accelerator, cod, nr_turns, p0, e_nrpts, e_min, e_max, x_nrpts, x_min, x_max, true, grid, nr_threads);

  // generates output files
  std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
  status = print_closed_orbit(accelerator, cod);
  if (status == Status::file_not_opened) return status;
  std::cout << get_timestamp() << " saving dynap_exfmap grid to file" << std::endl;
  status = print_dynapgrid (accelerator, grid, "[dynap_fmap]", "dynap_exfmap_out.txt", true);
  if (status == Status::file_not_opened) return status;

  std::cout << get_timestamp() << " end timestamp" << std::endl;
  return EXIT_SUCCESS;

}

int cmd_track_linepass(const std::vector<std::string>& args) {

  if (args.size() != 15) {
    std::cerr << "cmd_track_linepass: invalid number of arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "[cmd_track_linepass]" << std::endl << std::endl;

  std::string  flat_filename(args[2]);
  double       ring_energy     = std::atof(args[3].c_str());
  unsigned int harmonic_number = std::atoi(args[4].c_str());
  std::string  cavity_state(args[5]);
  std::string  radiation_state(args[6]);
  std::string  vchamber_state(args[7]);
  unsigned int start_element = std::atoi(args[8].c_str());
  double       rx0 = std::atof(args[9].c_str());
  double       px0 = std::atof(args[10].c_str());
  double       ry0 = std::atof(args[11].c_str());
  double       py0 = std::atof(args[12].c_str());
  double       de0 = std::atof(args[13].c_str());
  double       dl0 = std::atof(args[14].c_str());

  print_header(stdout);
  std::cout << std::endl;
  std::cout << "flat_filename   : " << flat_filename << std::endl;
  std::cout << "energy[eV]      : " << ring_energy << std::endl;
  std::cout << "harmonic_number : " << harmonic_number << std::endl;
  std::cout << "cavity_state    : " << cavity_state << std::endl;
  std::cout << "radiation_state : " << radiation_state << std::endl;
  std::cout << "vchamber_state  : " << vchamber_state << std::endl;
  std::cout << "start_element   : " << start_element << std::endl;
  std::cout << "rx0[m]          : " << rx0 << std::endl;
  std::cout << "px0[rad]        : " << px0 << std::endl;
  std::cout << "ry0[m]          : " << ry0 << std::endl;
  std::cout << "py0[rad]        : " << py0 << std::endl;
  std::cout << "de0             : " << de0 << std::endl;
  std::cout << "dl0[m]          : " << dl0 << std::endl;

  std::cout << std::endl;
  std::cout << get_timestamp() << " begin timestamp" << std::endl;

  Accelerator accelerator;

  // reads flat file
  Status::type status = read_flat_file(flat_filename, accelerator);
  if (status == Status::file_not_found) {
    std::cerr << "track_linepass: flat file not found!" << std::endl;
    return EXIT_FAILURE;
  }
  if (status == Status::flat_file_error) {
    std::cerr << "track_linepass: flat file with incorrect syntax!" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << get_timestamp() << " input file with flat lattice read." << std::endl;

  // builds accelerator
  accelerator.energy = ring_energy;
  accelerator.harmonic_number = harmonic_number;
  accelerator.cavity_on = (cavity_state == "on");
  accelerator.radiation_on = (radiation_state == "on") ? RadiationState::damping : RadiationState::off;
  accelerator.vchamber_on = (vchamber_state == "on");

  // does tracking
  Pos<double> pos(rx0,px0,ry0,py0,de0,dl0);
  std::vector<Pos<double>> pos_list;
  Plane::type lost_plane;
  unsigned int offset_element = start_element;
  // adjust dl to keep the arrival-time in sync with wall clock
  accelerator.update_time_aware_info();
  track_linepass(accelerator, pos, true, offset_element, pos_list, lost_plane
  );

  std::cout << get_timestamp() << " saving track_linepass data to file" << std::endl;
  status = print_tracking_linepass(accelerator, pos_list, start_element, "track_linepass_out.txt");
  if (status == Status::file_not_opened) return status;

  std::cout << get_timestamp() << " end timestamp" << std::endl;
  return EXIT_SUCCESS;

}

int cmd_help(const std::vector<std::string>& args) {

  std::vector<std::string> help_help = {
    "help",
    "====",
    "possible commands:",
    "",
    "help:           this help. accepts secong argument for detailed help",
    "dynap_xy:       calculates dynamical aperture in xy plane",
    "dynap_ex:       calculates dynamical aperture in ex plane",
    "dynap_ma:       calculates momentum acceptance along the ring",
    "dynap_pxa:      calculates maximum px along the ring",
    "dynap_pya:      calculates maximum px along the ring",
    "dynap_xyfmap:   calculates fmap in xy plane",
    "dynap_exfmap:   calculates fmap in ex plane",
    "track_linepass: does one turn tracking of a initial condition",
    "tests:          used for debugging and testing trackcpp"
  };

  std::vector<std::string> dynap_ma_help = {
    "dynap_ma",
    "========",
    "calculates momentum acceptance along the ring",
    "both positive and negative momentum acceptances are searched"
    "",
    "input parameters (in order):",
    "",
    "flat_filename   : name of the file with flat lattice",
    "energy[eV]      : energy of the electron beam",
    "harmonic_number : harmonic number of the lattice",
    "cavity_state    : on|off",
    "radiation_state : on|off",
    "vchamber_state  : on|off",
    "nr_turns        : number of turns to track for each initial condition",
    "y0[m]           : initial y position (typically a small number <> 0)",
    "e_init          : initial energy for search",
    "e_delta         : intial delta energy",
    "nr_steps_back   : number of current energy delta to step back when an unstable point if found",
    "rescale         : rescaling factor for e_delta in new iterations",
    "nr_iterations   : number of iterations",
    "s_min[m]        : calcs momentum acceptance limited at s=[s_min,s_max] and at elements belonging to 'fam_names'",
    "s_max[m]        : calcs momentum acceptance limited at s=[s_min,s_max]",
    "nr_threads      : number of threads to run in parallel",
    "fam_names       : calcs momentum acceptance limited at s=[s_min,s_max] and at elements belonging to 'fam_names'",
  };

  std::vector<std::string> default_text = {"trackcpp: no help available."};
  std::vector<std::string>* text = &default_text;
  if (args.size() == 2) {
    text = &help_help;
  } if (args.size() > 2) {
    if (args[2] == "dynap_ma") text = &dynap_ma_help;
  }
  print_header (stdout); std::cout << std::endl;
  for(unsigned int i=0; i<text->size(); ++i) std::cout << (*text)[i] << std::endl;

  return EXIT_SUCCESS;

}
