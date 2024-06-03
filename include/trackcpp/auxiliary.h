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
#include <random>

class PassMethodsClass {
public:
    static const int pm_identity_pass                  = 0;
    static const int pm_drift_pass                     = 1;
    static const int pm_str_mpole_symplectic4_pass     = 2;
    static const int pm_bnd_mpole_symplectic4_pass     = 3;
    static const int pm_corrector_pass                 = 4;
    static const int pm_cavity_0_pass                  = 5;
    static const int pm_thinquad_pass                  = 6;
    static const int pm_thinsext_pass                  = 7;
    static const int pm_kickmap_pass                   = 8;
    static const int pm_matrix_pass                    = 9;
    static const int pm_drift_g2l_pass                 = 10;
    static const int pm_cavity_1comp_pass              = 11;
    static const int pm_cavity_1frac_pass              = 12;
    static const int pm_cavity_2_pass                  = 13;
    static const int pm_nr_pms                         = 14;  // counter for number of passmethods
    PassMethodsClass() {
        passmethods.push_back("identity_pass");
        passmethods.push_back("drift_pass");
        passmethods.push_back("str_mpole_symplectic4_pass");
        passmethods.push_back("bnd_mpole_symplectic4_pass");
        passmethods.push_back("corrector_pass");
        passmethods.push_back("cavity_0_pass");
        passmethods.push_back("thinquad_pass");
        passmethods.push_back("thinsext_pass");
        passmethods.push_back("kicktable_pass");
        passmethods.push_back("matrix_pass");
        passmethods.push_back("drift_g2l_pass");
        passmethods.push_back("cavity_1comp_pass");
        passmethods.push_back("cavity_1frac_pass");
        passmethods.push_back("cavity_2_pass");
    }
    int size() const { return passmethods.size(); }
    std::string operator[](const int i) const { return passmethods[i]; }
    bool is_time_aware_pm(const int i) const {
        bool flag = false;
        switch (i)
        {
        case pm_cavity_0_pass:
            flag = true;
            break;
        case pm_cavity_1comp_pass:
            flag = true;
            break;
        case pm_cavity_1frac_pass:
            flag = true;
            break;
        case pm_cavity_2_pass:
            flag = true;
            break;
        default:
            break;
        }
        return flag;
    };
private:
    std::vector<std::string> passmethods;
};

// this is to be superseeded by PassMethodClass
struct PassMethod {
    enum type {
        pm_identity_pass                  = 0,
        pm_drift_pass                     = 1,
        pm_str_mpole_symplectic4_pass     = 2,
        pm_bnd_mpole_symplectic4_pass     = 3,
        pm_corrector_pass                 = 4,
        pm_cavity_0_pass                  = 5,
        pm_thinquad_pass                  = 6,
        pm_thinsext_pass                  = 7,
        pm_kickmap_pass                   = 8,
        pm_matrix_pass                    = 9,
        pm_drift_g2l_pass                 = 10,
        pm_cavity_1comp_pass              = 11,
        pm_cavity_1frac_pass              = 12,
        pm_cavity_2_pass                  = 13,
        pm_nr_pms                         = 14,
    };
};

// this is to be superseede by PassMethodClass
const std::vector<std::string> pm_dict = {
    // this vector has to have the same number of entries as enum above
    "identity_pass",
    "drift_pass",
    "str_mpole_symplectic4_pass",
    "bnd_mpole_symplectic4_pass",
    "corrector_pass",
    "cavity_0_pass",
    "thinquad_pass",
    "thinsext_pass",
    "kicktable_pass",
    "matrix_pass",
    "drift_g2l_pass",
    "cavity_1comp_pass",
    "cavity_1frac_pass",
    "cavity_2_pass",
};

struct RadiationState {
    enum type {
        off                    = 0,
        damping                = 1,
        full                   = 2,
    };
};

const std::vector<std::string> rad_dict = {
    "off",
    "damping",
    "full"
};

struct Distributions {
    enum type {
        normal =  0,
        uniform  =  1,
    };
};

const std::vector<std::string> distributions_dict = {
    "normal",
    "uniform",
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
        z = 3,
        xy = 4
    };
};

struct VChamberShape {
    // integer corresponds to p-norm number defining the normalized shape
    enum type {
        rhombus = 1,
        ellipse = 2,
        rectangle = 0,  // (inf-norm)
    };
};

// See https://en.wikipedia.org/wiki/International_System_of_Units

const double light_speed              = 299792458;         // [m/s]   - definition
const double electron_charge          = 1.602176634e-19;   // [C]     - definition
const double reduced_planck_constant  = 1.054571817e-34;   // [J.s]   - definition
const double electron_mass            = 9.1093837015e-31;  // [Kg]    - 2022-03-19 - https://physics.nist.gov/cgi-bin/cuu/Value?me|search_for=electron+mass
const double vacuum_permeability      = 1.25663706212e-6;  // [T.m/A] - 2022-03-19 - https://physics.nist.gov/cgi-bin/cuu/Value?mu0|search_for=vacuum+permeability
const double electron_rest_energy     = electron_mass * pow(light_speed,2);      // [Kg.m^2/s^2] - derived
const double electron_rest_energy_eV  = electron_rest_energy / electron_charge;  // [eV] - derived
const double vacuum_permitticity      = 1/(vacuum_permeability * pow(light_speed,2));   // [V.s/(A.m)]  - derived
const double electron_radius          = pow(electron_charge,2)/(4*M_PI*vacuum_permitticity*electron_rest_energy); // [m] - derived
const double id6[6] = {0, 0, 0, 0, 0, 0}; // translation identity
const double id66[36] = {
    1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1
}; // rotation identity

template <typename T>
int sgn(T val) {
    if (val >= 0) return 1; else return -1;
}


void set_random_seed(unsigned rnd_seed);
void set_random_seed_with_random_device();
void set_random_distribution(unsigned value);
double gen_random_number();


#if __GNUC__ < 6
    bool isfinite(const double& v);
#endif

std::string get_timestamp();

//typedef double Vector4[4];
//typedef double Vector6[6];
//typedef double Matrix4[4][4];

//typedef std::vector<double> Vector;
//typedef std::vector<std::vector<double> > Matrix;

#undef sqr

#endif
