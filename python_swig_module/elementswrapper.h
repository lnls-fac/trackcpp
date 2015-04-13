#ifndef ELEMENTSWRAPPER_H
#define ELEMENTSWRAPPER_H

#include "elements.h"

Element marker_wrapper(const std::string& fam_name_);
Element bpm_wrapper(const std::string& fam_name_);
Element hcorrector_wrapper(const std::string& fam_name_, const double& length_, const double& hkick_);
Element vcorrector_wrapper(const std::string& fam_name_, const double& length_, const double& vkick_);
Element corrector_wrapper(const std::string& fam_name_, const double& length_, const double& hkick_, const double& vkick_);
Element drift_wrapper(const std::string& fam_name_, const double& length_);
Element rbend_wrapper(const std::string& fam_name_, const double& length_,
					  const double& angle_, const double& angle_in_, const double& angle_out_,
					  const double& gap_, const double& fint_in_, const double& fint_out_,
					  const std::vector<double>& polynom_a_, const std::vector<double>& polynom_b_,
					  const double& K_, const double& S_);
Element quadrupole_wrapper(const std::string& fam_name_, const double& length_, const double& K_, const int nr_steps_ = 1);
Element sextupole_wrapper(const std::string& fam_name_, const double& length_, const double& S_, const int nr_steps_ = 1);
Element rfcavity_wrapper(const std::string& fam_name_, const double& length_, const double& frequency_, const double& voltage_);

#endif // ELEMENTSWRAPPER_H
