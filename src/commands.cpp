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
// #include <trackcpp/dynap.h>
#include <trackcpp/flat_file.h>
#include <trackcpp/tracking.h>
#include <trackcpp/lattice.h>
#include <trackcpp/accelerator.h>
#include <trackcpp/elements.h>
#include <trackcpp/auxiliary.h>
#include <string>
#include <cstdlib>
#include <iostream>

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
  print_header(); std::cout << std::endl;
  for(unsigned int i=0; i<text->size(); ++i) std::cout << (*text)[i] << std::endl;

  return EXIT_SUCCESS;

}
