// TRACKC++
// ========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		Tue Dec 10 17:57:20 BRST 2013


#include "auxiliary.h"
#include "elements.h"
#include <cfloat>

const std::vector<double> Element::default_polynom = std::vector<double>(3,0);


// default constructor (constructs a drift element)
Element::Element(const std::string& fam_name_, const double& length_) :
    fam_name(fam_name_), length(length_) {
	for(unsigned int i=0; i<6; i++) {
		t_in[i] = t_out[i] = 0.0;
		for(unsigned int j=0; j<6; ++j) {
			if (i == j) {
				r_in[i*6+j] = r_out[i*6+j] = 1.0;
			} else {
				r_in[i*6+j] = r_out[i*6+j] = 0.0;
			}
		}
	}

}

const std::string& Element::get_pass_method() {
    return pm_dict[pass_method];
}

void Element::set_pass_method(const std::string &pass_method_) {
    int i;
    for(i = 0; i<pm_dict.size(); i++)
        if (pm_dict[i] == pass_method_)
            break;
    if (i < pm_dict.size())
        pass_method = i;
}


Element Element::marker (const std::string& fam_name_) {
	Element e = Element(fam_name_, 0);
    initialize_marker(e);
	return e;
}

Element Element::bpm (const std::string& fam_name_) {
    Element e = Element(fam_name_, 0);
    initialize_marker(e);
    return e;
}

Element Element::drift (const std::string& fam_name_, const double& length_) {
	Element e = Element(fam_name_, length_);
    initialize_drift(e);
	return e;
}

Element Element::hcorrector(const std::string& fam_name_, const double& length_, const double& hkick_) {
	Element e = Element(fam_name_, length_);
    initialize_corrector(e, hkick_, 0.0);
	return e;
}

Element Element::vcorrector(const std::string& fam_name_, const double& length_, const double& vkick_) {
	Element e = Element(fam_name_, length_);
    initialize_corrector(e, 0.0, vkick_);
	return e;
}

Element Element::corrector(const std::string& fam_name_, const double& length_, const double& hkick_, const double& vkick_) {
	Element e = Element(fam_name_, length_);
    initialize_corrector(e, hkick_, vkick_);
	return e;
}

Element Element::quadrupole (const std::string& fam_name_, const double& length_, const double& K_, const int nr_steps_) {
	Element e = Element(fam_name_, length_);
    initialize_quadrupole(e, K_, nr_steps_);
	return e;
}

Element Element::sextupole (const std::string& fam_name_, const double& length_, const double& S_, const int nr_steps_) {
	Element e = Element(fam_name_, length_);
    initialize_sextupole(e, S_, nr_steps_);
	return e;
}

Element Element::rbend (const std::string& fam_name_, const double& length_,
		const double& angle_, const double& angle_in_, const double& angle_out_,
		const double& gap_, const double& fint_in_, const double& fint_out_,
		const std::vector<double>& polynom_a_, const std::vector<double>& polynom_b_,
		const double& K_, const double& S_, const int nr_steps_) {
    Element e = Element(fam_name_, length_);
    initialize_rbend(e, angle_, angle_in_, angle_out_, gap_, fint_in_, fint_out_, polynom_a_, polynom_b_, K_, S_, nr_steps_);
	return e;
}

Element Element::rfcavity (const std::string& fam_name_, const double& length_, const double& frequency_, const double& voltage_) {
	Element e = Element(fam_name_, length_);
    initialize_rfcavity(e, frequency_, voltage_);
	return e;
}



void print_polynom(std::ostream& out, const std::string& label, const std::vector<double>& polynom) {
	int order = 0;
	for(unsigned int i=0; i<polynom.size(); ++i) {
		if (polynom[i] != 0) order = i+1;
	}
	if (order > 0) out << label;
	for(int i=0; i<order; ++i) {
		out << polynom[i] << " ";
	}
	if (order > 0) out << std::endl;
}

std::ostream& operator<< (std::ostream &out, const Element& el) {

	                      out << "fam_name      : " << el.fam_name << std::endl;
	if (el.length != 0)   out << "length        : " << el.length << std::endl;
	                      out << "pass_method   : " << pm_dict[el.pass_method] << std::endl;
	if (el.nr_steps > 1)  out << "nr_steps      : " << el.nr_steps << std::endl;
	if (el.thin_KL != 0)  out << "thin_KL       : " << el.thin_KL << std::endl;
	if (el.thin_SL != 0)  out << "thin_SL       : " << el.thin_SL << std::endl;
	if (el.angle != 0)    out << "bending_angle : " << el.angle << std::endl;
	if (el.angle != 0)    out << "entrance_angle: " << el.angle_in << std::endl;
	if (el.angle != 0)    out << "exit_angle    : " << el.angle_out << std::endl;
	if ((el.gap != 0) and ((el.fint_in != 0) or (el.fint_out != 0))) {
		                  out << "gap           : " << el.gap << std::endl;
		                  out << "fint_in       : " << el.fint_in << std::endl;
		                  out << "fint_out      : " << el.fint_out << std::endl;
	}
	print_polynom(        out,   "polynom_a     : ", el.polynom_a);
	print_polynom(        out,   "polynom_b     : ", el.polynom_b);
	if (el.frequency != 0)out << "frequency     : " << el.frequency << std::endl;
	if (el.voltage != 0)  out << "voltage       : " << el.voltage << std::endl;
	return out;
}

void initialize_marker(Element &element) {
    element.pass_method = PassMethod::pm_identity_pass;
}

void initialize_corrector(Element &element, const double &hkick, const double &vkick) {
    element.pass_method = PassMethod::pm_corrector_pass;
    element.hkick = hkick;
    element.vkick = vkick;
}

void initialize_drift(Element &element) {
    element.pass_method = PassMethod::pm_drift_pass;
}

void initialize_rbend(Element& element, const double& angle, const double& angle_in, const double& angle_out, const double& gap, const double& fint_in, const double& fint_out, const std::vector<double>& polynom_a, const std::vector<double>& polynom_b, const double& K, const double& S, const int nr_steps) {
    element.pass_method = PassMethod::pm_bnd_mpole_symplectic4_pass;
    element.angle = angle;
    element.angle_in = angle_in;
    element.angle_out = angle_out;
    element.gap = gap;
    element.fint_in = fint_in;
    element.fint_out = fint_out;
    element.polynom_a = polynom_a;
    element.polynom_b = polynom_b;
    element.polynom_b[1] = K;
    element.polynom_b[2] = S;
    element.nr_steps = nr_steps;
}



void initialize_quadrupole(Element &element, const double &K, const int &nr_steps) {
    element.pass_method = PassMethod::pm_str_mpole_symplectic4_pass;
    element.polynom_b[1] = K;
    element.nr_steps = nr_steps;
}

void initialize_sextupole(Element &element, const double &S, const int &nr_steps) {
    element.pass_method = PassMethod::pm_str_mpole_symplectic4_pass;
    element.polynom_b[2] = S;
    element.nr_steps = nr_steps;
}

void initialize_rfcavity(Element &element, const double &frequency, const double &voltage) {
    element.pass_method = PassMethod::pm_cavity_pass;
    element.frequency = frequency;
    element.voltage = voltage;
}
