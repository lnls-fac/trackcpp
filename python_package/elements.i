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

class Element {
public:
    std::string fam_name;
    int pass_method;
    double length;
    int nr_steps;
    double hkick;
    double vkick;
    double angle;
    double angle_in;
    double angle_out;
    double gap;
    double fint_in;
    double fint_out;
    double thin_KL;
    double thin_SL;
    double frequency;
    double voltage;
    std::vector<double> polynom_a;
    std::vector<double> polynom_b;
    const Kicktable* kicktable;
    double hmin;
    double hmax;
    double vmin;
    double vmax;
    double t_in[6], t_out[6];
    double r_in[36], r_out[36];

    Element(const std::string& fam_name_, const double& length_);
    Element(const Element&);

    const std::string& get_pass_method();
    void set_pass_method(const std::string& pass_method_);
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
