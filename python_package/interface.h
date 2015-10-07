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
    Plane::type  lost_plane;
    bool         trajectory;
};

struct RingPassArgs : public LinePassArgs {
    unsigned int nr_turns;
    unsigned int lost_turn;
    unsigned int element_offset;
    Plane::type  lost_plane;
    unsigned int trajectory;
};

struct String {
  public: std::string data;
};

Status::type track_elementpass_wrapper (
        const Element& el,
        Pos<double> &orig_pos,
        const Accelerator& accelerator);

Status::type track_linepass_wrapper (
        const Accelerator& accelerator,
        Pos<double>& orig_pos,
        std::vector< Pos<double> >& pos,
        LinePassArgs& args);

Status::type track_ringpass_wrapper (
        const Accelerator& accelerator,
        Pos<double> &orig_pos,
        std::vector< Pos<double> >& pos,
        RingPassArgs& args);

Element marker_wrapper(const std::string& fam_name_);
Element bpm_wrapper(const std::string& fam_name_);
Element hcorrector_wrapper(const std::string& fam_name_, const double& length_, const double& hkick_);
Element vcorrector_wrapper(const std::string& fam_name_, const double& length_, const double& vkick_);
Element corrector_wrapper(const std::string& fam_name_, const double& length_, const double& hkick_, const double& vkick_);
Element drift_wrapper(const std::string& fam_name_, const double& length_);
Element quadrupole_wrapper(const std::string& fam_name_, const double& length_, const double& K_, const int nr_steps_ = 1);
Element sextupole_wrapper(const std::string& fam_name_, const double& length_, const double& S_, const int nr_steps_ = 1);
Element rfcavity_wrapper(const std::string& fam_name_, const double& length_, const double& frequency_, const double& voltage_);
Element rbend_wrapper(const std::string& fam_name_, const double& length_,
                      const double& angle_, const double& angle_in_, const double& angle_out_,
                      const double& gap_, const double& fint_in_, const double& fint_out_,
                      const std::vector<double>& polynom_a_, const std::vector<double>& polynom_b_,
                      const double& K_, const double& S_);
Status::type write_flat_file_wrapper(String& fname, const Accelerator& accelerator, bool file_flag = true);
Status::type read_flat_file_wrapper(String& fname, Accelerator& accelerator, bool file_flag = true);

#endif // INTERFACE_H
