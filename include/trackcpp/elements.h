// The MIT License (MIT)
//
// Copyright (c) 2015 LNLS Accelerator Division
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "kicktable.h"
#include "auxiliary.h"
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
//    frequency(0), voltage(0),
//    polynom_a(default_polynom), polynom_b(default_polynom),
//    hmax(DBL_MAX), vmax(DBL_MAX)

  std::string   fam_name;
  int           pass_method = PassMethod::pm_drift_pass;
  double        length      = 0;  // [m]
  int           nr_steps    = 1;  //
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

  std::vector<double> polynom_a = default_polynom;
  std::vector<double> polynom_b = default_polynom;

  const Kicktable*    kicktable   = nullptr;

  double              t_in[6],  t_out[6];
  double              r_in[36], r_out[36];

  const std::string& get_pass_method();
  void set_pass_method(const std::string& pass_method_);


  // default constructor (builds a drift-type element)
  Element(const std::string& fam_name_ = "", const double& length_ = 0);

  static const std::vector<double> default_polynom;

  // front-end routines for typed element creation
  static Element marker     (const std::string& fam_name_);
  static Element bpm        (const std::string& fam_name_);
  static Element hcorrector (const std::string& fam_name_, const double& length_, const double& hkick_);
  static Element vcorrector (const std::string& fam_name_, const double& length_, const double& vkick_);
  static Element corrector  (const std::string& fam_name_, const double& length_, const double& hkick_, const double& vkick_);
  static Element drift      (const std::string& fam_name_, const double& length_);
  static Element rbend      (const std::string& fam_name_, const double& length_,
                             const double& angle_, const double& angle_in_ = 0, const double& angle_out_ = 0,
                             const double& gap_ = 0, const double& fint_in_ = 0, const double& fint_out_ = 0,
                             const std::vector<double>& polynom_a_ = default_polynom, const std::vector<double>& polynom_b_ = default_polynom,
                             const double& K_ = 0, const double& S_ = 0, const int nr_steps_ = 20);
  static Element quadrupole (const std::string& fam_name_, const double& length_, const double& K_, const int nr_steps_ = 10);
  static Element sextupole  (const std::string& fam_name_, const double& length_, const double& S_, const int nr_steps_ = 5);
  static Element rfcavity   (const std::string& fam_name_, const double& length_, const double& frequency_, const double& voltage_);

  bool operator==(const Element& o) const;
  bool operator!=(const Element& o) const { return !(*this == o); };

  friend std::ostream& operator<< (std::ostream &out, const Element& el);

};

void initialize_marker(Element& element);
void initialize_corrector(Element& element, const double& hkick, const double& vkick);
void initialize_drift(Element& element);
void initialize_rbend(Element& element, const double& angle, const double& angle_in, const double& angle_out,
            const double& gap, const double& fint_in, const double& fint_out,
            const std::vector<double>& polynom_a, const std::vector<double>& polynom_b,
            const double& K, const double& S, const int nr_steps);
void initialize_quadrupole(Element& element, const double& K, const int& nr_steps);
void initialize_sextupole(Element& element, const double& S, const int& nr_steps);
void initialize_rfcavity(Element& element, const double& frequency, const double& voltage);

#endif
