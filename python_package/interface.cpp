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

#include "interface.h"
#include <trackcpp/flat_file.h>


Status::type track_elementpass_wrapper (
         const Element& el,
         Pos<double> &orig_pos,
         const Accelerator& accelerator) {
         return track_elementpass (el,
                                   orig_pos,
                                   accelerator);
}

Status::type track_linepass_wrapper(
        const Accelerator &accelerator,
        Pos<double> &orig_pos,
        std::vector< Pos<double> >& pos,
        LinePassArgs& args) {
    return track_linepass(accelerator,
                          orig_pos,
                          pos,
                          args.element_offset,
                          args.lost_plane,
                          args.trajectory);
}

Status::type track_ringpass_wrapper (
        const Accelerator& accelerator,
        Pos<double> &orig_pos,
        std::vector< Pos<double> >& pos,
        RingPassArgs& args) {
    return track_ringpass(accelerator,
                          orig_pos,
                          pos,
                          args.nr_turns,
                          args.lost_turn,
                          args.element_offset,
                          args.lost_plane,
                          args.trajectory);
}



Element marker_wrapper(const std::string &fam_name_) {
    return Element::marker(fam_name_);
}

Element bpm_wrapper(const std::string &fam_name_) {
    return Element::bpm(fam_name_);
}

Element hcorrector_wrapper(const std::string &fam_name_, const double &length_, const double &hkick_) {
    return Element::hcorrector(fam_name_, length_, hkick_);
}

Element vcorrector_wrapper(const std::string &fam_name_, const double &length_, const double &vkick_) {
    return Element::vcorrector(fam_name_, length_, vkick_);
}

Element corrector_wrapper(const std::string &fam_name_, const double &length_, const double &hkick_, const double &vkick_) {
    return Element::corrector(fam_name_, length_, hkick_, vkick_);
}

Element drift_wrapper(const std::string &fam_name_, const double &length_) {
    return Element::drift(fam_name_, length_);
}

Element matrix_wrapper(const std::string &fam_name_, const double &length_) {
    return Element::matrix(fam_name_, length_);
}

Element rbend_wrapper(const std::string& fam_name_, const double& length_,
					  const double& angle_, const double& angle_in_, const double& angle_out_,
					  const double& gap_, const double& fint_in_, const double& fint_out_,
					  const std::vector<double>& polynom_a_, const std::vector<double>& polynom_b_,
					  const double& K_, const double& S_) {
    return Element::rbend(fam_name_, length_, angle_, angle_in_, angle_out_, gap_, fint_in_, fint_out_, polynom_a_, polynom_b_, K_, S_);
}

Element quadrupole_wrapper(const std::string &fam_name_, const double &length_, const double &K_, const int nr_steps_) {
    return Element::quadrupole(fam_name_, length_, K_, nr_steps_);
}

Element sextupole_wrapper(const std::string &fam_name_, const double &length_, const double &S_, const int nr_steps_) {
    return Element::sextupole(fam_name_, length_, S_, nr_steps_);
}

Element rfcavity_wrapper(const std::string &fam_name_, const double &length_, const double &frequency_, const double &voltage_, const double &phase_lag_) {
    return Element::rfcavity(fam_name_, length_, frequency_, voltage_, phase_lag_);
}

Status::type read_flat_file_wrapper(String& fname, Accelerator& accelerator, bool file_flag) {
  return read_flat_file(fname.data, accelerator, file_flag);
}

Status::type write_flat_file_wrapper(String& fname, const Accelerator& accelerator, bool file_flag) {
  return write_flat_file(fname.data, accelerator, file_flag);
}
