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

#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "kicktable.h"
#include "auxiliary.h"
#include "linalg.h"
#include <vector>
#include <string>
#include <fstream>
#include <cfloat>


class Element {
public:

//  fam_name(fam_name_), pass_method(PassMethod::pm_drift_pass),
//    nr_steps(1), length(length_),
//    hkick(0), vkick(0),
//    angle(0), angle_in(0), angle_out(0),
//    gap(0), fint_in(0), fint_out(0),
//    thin_KL(0), thin_SL(0),
//    frequency(0), voltage(0), phase_lag(0),
//    polynom_a(default_polynom), polynom_b(default_polynom),
//    hmax(DBL_MAX), vmax(DBL_MAX)

  std::string   fam_name;
  int           pass_method = PassMethod::pm_drift_pass;
  double        length      = 0;  // [m]
  int           nr_steps    = 1;  //
  int           vchamber    = VChamberShape::rectangle;  // vchamber type
  double        hmin        = -DBL_MAX;  // [m]
  double        hmax        =  DBL_MAX;  // [m]
  double        vmin        = -DBL_MAX;  // [m]
  double        vmax        =  DBL_MAX;  // [m]
  double        hkick       = 0;  // [rad]
  double        vkick       = 0;  // [rad]
  double        angle       = 0;  // [rad]
  double        angle_in    = 0;  // [rad]
  double        angle_out   = 0;  // [rad]
  double        gap         = 0;  // [m]
  double        fint_in     = 0;
  double        fint_out    = 0;
  double        thin_KL     = 0;  // [1/m]
  double        thin_SL     = 0;  // [1/m²]
  double        frequency   = 0;  // [Hz]
  double        voltage     = 0;  // [V]
  double        phase_lag   = 0;  // [rad]
  int           kicktable_idx = -1;   // index of kickmap object in kicktable_list
  double        rescale_kicks = 1.0;  // for kickmaps
  double        ks = 0;       // [1/m]
  double        kx = 0;       // [1/m]
  double        s0 = 0;       // [m]
  std::vector<std::vector<double>> coefs;
  
  std::vector<double> polynom_a = default_polynom;
  std::vector<double> polynom_b = default_polynom;
  Matrix              matrix66 = Matrix(6);

  double t_in[6];
  double t_out[6];
  double r_in[36];
  double r_out[36];

  bool has_t_in = false;
  bool has_t_out = false;
  bool has_r_in = false;
  bool has_r_out = false;

  void reflag_t_in(void);
  void reflag_t_out(void);
  void reflag_r_in(void);
  void reflag_r_out(void);

  const std::string& get_pass_method();
  void set_pass_method(const std::string& pass_method_);


  // default constructor (builds a drift-type element)
  Element(const std::string& fam_name_ = "", const double& length_ = 0);
  Element(const Element &) = default;

  static const std::vector<double> default_polynom;

  // front-end routines for typed element creation
  static Element marker     (const std::string& fam_name_);
  static Element bpm        (const std::string& fam_name_);
  static Element hcorrector (const std::string& fam_name_, const double& length_, const double& hkick_);
  static Element vcorrector (const std::string& fam_name_, const double& length_, const double& vkick_);
  static Element corrector  (const std::string& fam_name_, const double& length_, const double& hkick_, const double& vkick_);
  static Element drift      (const std::string& fam_name_, const double& length_);
  static Element drift_g2l  (const std::string& fam_name_, const double& length_);
  static Element matrix     (const std::string& fam_name_, const double& length_);
  static Element rbend      (const std::string& fam_name_, const double& length_,
                             const double& angle_, const double& angle_in_ = 0, const double& angle_out_ = 0,
                             const double& gap_ = 0, const double& fint_in_ = 0, const double& fint_out_ = 0,
                             const std::vector<double>& polynom_a_ = default_polynom, const std::vector<double>& polynom_b_ = default_polynom,
                             const double& K_ = 0, const double& S_ = 0, const int nr_steps_ = 20);
  static Element quadrupole (const std::string& fam_name_, const double& length_, const double& K_, const int nr_steps_ = 10);
  static Element sextupole  (const std::string& fam_name_, const double& length_, const double& S_, const int nr_steps_ = 5);
  static Element rfcavity   (const std::string& fam_name_, const double& length_, const double& frequency_, const double& voltage_, const double& phase_lag_);
  static Element kickmap    (const std::string& fam_name_, const std::string& kicktable_fname_, const int nr_steps_ = 20, const double& rescale_length_ = 1.0, const double& rescale_kicks_ = 1.0);
  static Element field3d    (const std::string& fam_name_, const double& length_, const double& s0_, const double& kx_, const double& ks_, const std::vector<std::vector<double>>& coefs_, const int nr_steps_ = 40);

  bool operator==(const Element& o) const;
  bool operator!=(const Element& o) const { return !(*this == o); };

  friend std::ostream& operator<< (std::ostream &out, const Element& el);

};

void initialize_marker(Element& element);
void initialize_corrector(Element& element, const double& hkick, const double& vkick);
void initialize_drift(Element& element);
void initialize_drift_g2l(Element& element);
void initialize_matrix(Element& element);
void initialize_rbend(Element& element, const double& angle, const double& angle_in, const double& angle_out,
            const double& gap, const double& fint_in, const double& fint_out,
            const std::vector<double>& polynom_a, const std::vector<double>& polynom_b,
            const double& K, const double& S, const int nr_steps);
void initialize_quadrupole(Element& element, const double& K, const int& nr_steps);
void initialize_sextupole(Element& element, const double& S, const int& nr_steps);
void initialize_rfcavity(Element& element, const double& frequency, const double& voltage, const double& phase_lag);
void initialize_kickmap(Element& element, const int& kicktable_idx, const int& nr_steps, const double &rescale_kicks);
void initialize_field3d(Element& element, const double& s0_, const double& kx_, const double& ks_, const std::vector<std::vector<double>>& coefs_, const int nr_steps_);

#endif
