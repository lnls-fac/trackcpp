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

#ifndef INTERFACE_H
#define INTERFACE_H

#include <vector>
#include <string>
#include <trackcpp/accelerator.h>
#include <trackcpp/elements.h>
#include <trackcpp/auxiliary.h>
#include <trackcpp/tracking.h>
#include <trackcpp/pos.h>

struct LinePassArgs {
    unsigned int element_offset;
    std::vector< unsigned int > indices;
    std::vector< unsigned int > lost_plane;
    std::vector< unsigned int > lost_element;
};

struct RingPassArgs : public LinePassArgs {
    unsigned int nr_turns;
    std::vector< unsigned int > lost_turn;
    unsigned int element_offset;
    std::vector< unsigned int > lost_plane;
    std::vector< unsigned int > lost_element;
    unsigned int trajectory;
};

struct String {
public:
  std::string data;
  String(const std::string& v = "") : data(v) {};
};

Status::type track_elementpass_wrapper (
        const Element& el,
        double *pos, int n1, int n2,
        const Accelerator& accelerator);

Status::type track_linepass_wrapper (
        const Accelerator& accelerator,
        double *orig_pos, int ni1, int ni2,
        double *pos, int n1, int n2,
        LinePassArgs& args);

Status::type track_ringpass_wrapper (
        const Accelerator& accelerator,
        double *orig_pos, int ni1, int ni2,
        double *pos, int n1, int n2,
        RingPassArgs& args);

Element marker_wrapper(const std::string& fam_name_);
Element bpm_wrapper(const std::string& fam_name_);
Element hcorrector_wrapper(const std::string& fam_name_, const double& length_, const double& hkick_);
Element vcorrector_wrapper(const std::string& fam_name_, const double& length_, const double& vkick_);
Element corrector_wrapper(const std::string& fam_name_, const double& length_, const double& hkick_, const double& vkick_);
Element drift_wrapper(const std::string& fam_name_, const double& length_);
Element matrix_wrapper(const std::string& fam_name_, const double& length_);
Element quadrupole_wrapper(const std::string& fam_name_, const double& length_, const double& K_, const int nr_steps_ = 1);
Element sextupole_wrapper(const std::string& fam_name_, const double& length_, const double& S_, const int nr_steps_ = 1);
Element rfcavity_wrapper(const std::string& fam_name_, const double& length_, const double& frequency_, const double& voltage_, const double& phase_lag);
Element rbend_wrapper(const std::string& fam_name_, const double& length_,
                      const double& angle_, const double& angle_in_, const double& angle_out_,
                      const double& gap_, const double& fint_in_, const double& fint_out_,
                      const std::vector<double>& polynom_a_, const std::vector<double>& polynom_b_,
                      const double& K_, const double& S_);
Status::type write_flat_file_wrapper(String& fname, const Accelerator& accelerator, bool file_flag = true);
Status::type read_flat_file_wrapper(String& fname, Accelerator& accelerator, bool file_flag = true);

#endif // INTERFACE_H
