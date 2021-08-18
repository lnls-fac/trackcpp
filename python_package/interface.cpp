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
         double *pos, int n1, int n2,
         const Accelerator& accelerator) {

    std::vector<Pos<double> > post;

    post.reserve(n2);
    for (unsigned int i=0; i<n2; ++i){
        post.emplace_back(
            pos[0*n2 + i], pos[1*n2 + i],
            pos[2*n2 + i], pos[3*n2 + i],
            pos[4*n2 + i], pos[5*n2 + i]);
    }
    Status::type status = track_elementpass(
        el, post, accelerator);

    for (unsigned int i=0; i<post.size(); ++i){
        pos[0*n2 + i] = post[i].rx; pos[1*n2 + i] = post[i].px;
        pos[2*n2 + i] = post[i].ry; pos[3*n2 + i] = post[i].py;
        pos[4*n2 + i] = post[i].de; pos[5*n2 + i] = post[i].dl;
    }
    return status;
}

Status::type track_linepass_wrapper(
        const Accelerator &accelerator,
        double *orig_pos, int ni1, int ni2,
        double *pos, int n1, int n2,
        LinePassArgs& args) {

    std::vector<Pos<double> > post;
    std::vector<Pos<double> > orig_post;

    orig_post.reserve(ni2);
    for (unsigned int i=0; i<ni2; ++i){
        orig_post.emplace_back(
            orig_pos[0*ni2 + i], orig_pos[1*ni2 + i],
            orig_pos[2*ni2 + i], orig_pos[3*ni2 + i],
            orig_pos[4*ni2 + i], orig_pos[5*ni2 + i]);
    }
    Status::type status = track_linepass(accelerator,
                          orig_post,
                          post,
                          args.element_offset,
                          args.lost_plane,
                          args.lost_element,
                          args.indices);
    for (unsigned int i=0; i<post.size(); ++i){
        pos[0*n2 + i] = post[i].rx; pos[1*n2 + i] = post[i].px;
        pos[2*n2 + i] = post[i].ry; pos[3*n2 + i] = post[i].py;
        pos[4*n2 + i] = post[i].de; pos[5*n2 + i] = post[i].dl;
    }
    return status;
}

Status::type track_ringpass_wrapper (
        const Accelerator& accelerator,
        double *orig_pos, int ni1, int ni2,
        double *pos, int n1, int n2,
        RingPassArgs& args) {

    std::vector<Pos<double> > post;
    std::vector<Pos<double> > orig_post;

    orig_post.reserve(ni2);
    for (unsigned int i=0; i<ni2; ++i){
        orig_post.emplace_back(
            orig_pos[0*ni2 + i], orig_pos[1*ni2 + i],
            orig_pos[2*ni2 + i], orig_pos[3*ni2 + i],
            orig_pos[4*ni2 + i], orig_pos[5*ni2 + i]);
    }
    Status::type status = track_ringpass(accelerator,
                          orig_post,
                          post,
                          args.nr_turns,
                          args.lost_turn,
                          args.element_offset,
                          args.lost_plane,
                          args.lost_element,
                          args.trajectory);
    for (unsigned int i=0; i<post.size(); ++i){
        pos[0*n2 + i] = post[i].rx; pos[1*n2 + i] = post[i].px;
        pos[2*n2 + i] = post[i].ry; pos[3*n2 + i] = post[i].py;
        pos[4*n2 + i] = post[i].de; pos[5*n2 + i] = post[i].dl;
    }
    return status;
}

Status::type calc_twiss_wrapper (
        const Accelerator& accelerator,
        const Pos<double>& fixed_point,
        Matrix& m66,
        double *twiss, int n1, int n2,
        Twiss twiss0) {

    std::vector<Pos<double> > post;
    std::vector<Pos<double> > orig_post;

    std::vector<Twiss> twiss_;

    Status::type status = calc_twiss(
        accelerator, fixed_point, m66, twiss_, twiss0);

    if (status != Status::success) return status;
    for (unsigned int i=0; i<n1; ++i){
        unsigned int j = i*n2;
        twiss[j] = twiss_[i].spos;
        j++; twiss[j] = twiss_[i].betax;
        j++; twiss[j] = twiss_[i].alphax;
        j++; twiss[j] = twiss_[i].mux;
        j++; twiss[j] = twiss_[i].betay;
        j++; twiss[j] = twiss_[i].alphay;
        j++; twiss[j] = twiss_[i].muy;
        j++; twiss[j] = twiss_[i].etax[0];
        j++; twiss[j] = twiss_[i].etax[1];
        j++; twiss[j] = twiss_[i].etay[0];
        j++; twiss[j] = twiss_[i].etay[1];
        j++; twiss[j] = twiss_[i].co.rx;
        j++; twiss[j] = twiss_[i].co.px;
        j++; twiss[j] = twiss_[i].co.ry;
        j++; twiss[j] = twiss_[i].co.py;
        j++; twiss[j] = twiss_[i].co.de;
        j++; twiss[j] = twiss_[i].co.dl;
    }
    return status;
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

Element kickmap_wrapper(const std::string& fam_name_,  const std::string& kicktable_fname_, const int nr_steps_, const double& rescale_length_, const double& rescale_kicks_) {
    return Element::kickmap(fam_name_, kicktable_fname_, nr_steps_, rescale_length_, rescale_kicks_);
}

Status::type read_flat_file_wrapper(String& fname, Accelerator& accelerator, bool file_flag) {
  return read_flat_file(fname.data, accelerator, file_flag);
}

Status::type write_flat_file_wrapper(String& fname, const Accelerator& accelerator, bool file_flag) {
  return write_flat_file(fname.data, accelerator, file_flag);
}

Status::type kicktable_getkicks_wrapper(const int& kicktable_idx, const double& rx, const double& ry, double& hkick__, double& vkick__) {
  return kicktable_getkicks(kicktable_idx, rx, ry, hkick__, vkick__);
}

Status::type track_findm66_wrapper(
    const Accelerator& accelerator,
    const Pos<double>& fixed_point,
    double *cumul_tm, int n1_tm, int n2_tm, int n3_tm,
    double *m66, int n1_m66, int n2_m66,
    Pos<double>& v0,
    std::vector<unsigned int >& indices) {

    std::vector<Matrix> vec_tm;
    Matrix vec_m66;

    Status::type status = track_findm66(
        accelerator, fixed_point, vec_tm, vec_m66, v0, indices);

    for (unsigned int i=0; i<vec_tm.size(); ++i){
        const Matrix& m = vec_tm[i];
        for (unsigned int j=0; j<n2_tm; ++j)
            for (unsigned int k=0; k<n3_tm; ++k)
                cumul_tm[(i*n2_tm + j)*n3_tm + k] = m[j][k];
        }
    for (unsigned int i=0; i<n1_m66; ++i)
        for (unsigned int j=0; j<n2_m66; ++j)
            m66[i*n2_m66 + j] = vec_m66[i][j];
}


Status::type track_diffusionmatrix_wrapper(
    const Accelerator& accelerator,
    const Pos<double>& fixed_point,
    double *cumul_tm, int n1_tm, int n2_tm, int n3_tm,
    double *bdiffmats, int n1_bd, int n2_bd, int n3_bd){

    std::vector<Matrix> vec_tm;
    std::vector<Matrix> vec_bd;

    for (unsigned int i=0; i<n1_tm; ++i){
        Matrix m (6);
        for (unsigned int j=0; j<n2_tm; ++j)
            for (unsigned int k=0; k<n3_tm; ++k)
                m[j][k] = cumul_tm[(i*n2_tm + j)*n3_tm + k];
        vec_tm.emplace_back(m);
    }

    Status::type status = track_diffusionmatrix(
                accelerator, fixed_point, vec_tm, vec_bd);

    for (unsigned int i=0; i<n1_bd; ++i){
        const Matrix& m = vec_bd[i];
        for (unsigned int j=0; j<n2_bd; ++j)
            for (unsigned int k=0; k<n3_bd; ++k)
                bdiffmats[(i*n2_bd + j)*n3_bd + k] = m[j][k];
    }
}


void naff_general_wrapper(
    double *re_in, int n1_re_in, int n2_re_in,
    double *im_in, int n1_im_in, int n2_im_in,
    bool is_real, int nr_ff, int win,
    double *ff_out, int n1_ff_out, int n2_ff_out,
    double *re_out, int n1_re_out, int n2_re_out,
    double *im_out, int n1_im_out, int n2_im_out){

    std::vector<double> vec_re_in;
    std::vector<double> vec_im_in;
    std::vector<double> vec_ff_out;
    std::vector<double> vec_re_out;
    std::vector<double> vec_im_out;

    vec_re_in.resize(n2_re_in);
    vec_im_in.resize(n2_re_in);
    vec_ff_out.resize(n2_ff_out, 0.0);
    vec_re_out.resize(n2_ff_out, 0.0);
    vec_im_out.resize(n2_ff_out, 0.0);
    for (unsigned int i=0; i<n1_re_in; ++i){
        for (unsigned int j=0; j<n2_re_in; ++j){
            vec_re_in[j] = re_in[i*n2_re_in + j];
            vec_im_in[j] = im_in[i*n2_re_in + j];
        }

        naff_general(
            vec_re_in, vec_im_in,
            is_real, nr_ff, win,
            vec_ff_out, vec_re_out, vec_im_out);

        for (unsigned int j=0; j<n2_ff_out; ++j){
            ff_out[i*n2_ff_out + j] = vec_ff_out[j];
            re_out[i*n2_ff_out + j] = vec_re_out[j];
            im_out[i*n2_ff_out + j] = vec_im_out[j];
        }
    }
}
