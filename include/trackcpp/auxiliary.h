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

#ifndef _AUXILIARY_H
#define _AUXILIARY_H

#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <cmath>

struct PassMethod {
    enum type {
        pm_identity_pass                  = 0,
        pm_drift_pass                     = 1,
        pm_str_mpole_symplectic4_pass     = 2,
        pm_bnd_mpole_symplectic4_pass     = 3,
        pm_corrector_pass                 = 4,
        pm_cavity_pass                    = 5,
        pm_thinquad_pass                  = 6,
        pm_thinsext_pass                  = 7,
        pm_kicktable_pass                 = 8,
        pm_nr_pms                         = 9
    };
};

const std::vector<std::string> pm_dict = {
        "identity_pass",
        "drift_pass",
        "str_mpole_symplectic4_pass",
        "bnd_mpole_symplectic4_pass",
        "corrector_pass",
        "cavity_pass",
        "thinquad_pass",
        "thinsext_pass",
        "kicktable_pass"
};

struct Status {
    enum type {
        success = 0,
        passmethod_not_defined = 1,
        passmethod_not_implemented = 2,
        particle_lost = 3,
        inconsistent_dimensions = 4,
        uninitialized_memory = 5,
        findorbit_not_converged = 6,
        findorbit_one_turn_matrix_problem = 7,
        file_not_found = 8,
        file_not_opened = 9,
        kicktable_not_defined = 10,
        kicktable_out_of_range = 11,
        flat_file_error = 12,
        newton_not_converged = 13,
        not_implemented = 14,
    };
};

const std::vector<std::string> string_error_messages = {
        "success",
        "passmethod_not_defined",
        "passmethod_not_implemented",
        "particle_lost",
        "inconsistent_dimensions",
        "uninitialized_memory",
        "findorbit_not_converged",
        "findorbit_one_turn_matrix_problem",
        "file_not_found",
        "file_not_opened",
        "kicktable_not_defined",
        "kicktable_out_of_range",
        "flat_file_error",
        "newton_not_converged",
        "not_implemented",
};

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define VER_STR STR(VERSION)

const std::string string_version = "TRACKCPP version " + std::string(VER_STR) + " (" + std::string(__DATE__) + " " + std::string(__TIME__) + ")";

struct Plane {
	enum type {
		no_plane = 0,
		x = 1,
		y = 2,
		z = 3
	};
};

// template <typename T> class Pos;
// class Element;

const double light_speed              = 299792458;         // [m/s]   - definition
const double vacuum_permeability      = 4*M_PI*1e-7;       // [T.m/A] - definition
const double electron_charge          = 1.60217656535e-19; // [C]     - 2014-06-11
const double electron_mass            = 9.1093829140e-31;  // [Kg]    - 2014-06-11
const double electron_rest_energy     = electron_mass * pow(light_speed,2);             // [Kg.m^2/s^2] - derived
const double vaccum_permitticity      = 1/(vacuum_permeability * pow(light_speed,2));   // [V.s/(A.m)]  - derived
const double electron_rest_energy_MeV = (electron_rest_energy / electron_charge) / 1e6; // [MeV] - derived
const double electron_radius          = pow(electron_charge,2)/(4*M_PI*vaccum_permitticity*electron_rest_energy); // [m] - derived

template <typename T>
int sgn(T val) {
	if (val >= 0) return 1; else return -1;
}

bool isfinite(const double& v);
std::string get_timestamp();

//typedef double Vector4[4];
//typedef double Vector6[6];
//typedef double Matrix4[4][4];

//typedef std::vector<double> Vector;
//typedef std::vector<std::vector<double> > Matrix;

#undef sqr

#endif
