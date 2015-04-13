#include "flat_file.h"
#include "elements.h"
#include "auxiliary.h"
#include <fstream>
#include <string>
#include <sstream>

Status::type read_flat_file(const std::string& filename, Accelerator& accelerator) {
	return read_flat_file_trackcpp(filename, accelerator);
}

static void synchronize_polynomials(Element& e);


Status::type read_flat_file_trackcpp(const std::string& filename, Accelerator& accelerator) {

	std::ifstream fp(filename.c_str());
	if (fp.fail()) return Status::file_not_found;

	accelerator.lattice.clear();

	Element e;
	std::string line;
	unsigned int line_count = 0;
	while (std::getline(fp, line)) {
		line_count++;
		std::istringstream ss(line);
		std::string cmd;
		ss >> cmd;
		if (cmd[0] == '#') continue;
		if (cmd.compare("fam_name") == 0) {
			if (e.fam_name.compare("") != 0) {
				accelerator.lattice.push_back(e);
				e = Element();
			}
			ss >> e.fam_name;
			continue;
		}
		if (cmd.compare("length")      == 0) { ss >> e.length;    continue; }
		if (cmd.compare("hmax")        == 0) { ss >> e.hmax;      continue; }
		if (cmd.compare("vmax")        == 0) { ss >> e.vmax;      continue; }
		if (cmd.compare("hkick")       == 0) { ss >> e.hkick;     continue; }
		if (cmd.compare("vkick")       == 0) { ss >> e.vkick;     continue; }
		if (cmd.compare("nr_steps")    == 0) { ss >> e.nr_steps;  continue; }
		if (cmd.compare("angle")       == 0) { ss >> e.angle;     continue; }
		if (cmd.compare("gap")         == 0) { ss >> e.gap;       continue; }
		if (cmd.compare("fint_in")     == 0) { ss >> e.fint_in;   continue; }
		if (cmd.compare("fint_out")    == 0) { ss >> e.fint_out;  continue; }
		if (cmd.compare("voltage")     == 0) { ss >> e.voltage;   continue; }
		if (cmd.compare("frequency")   == 0) { ss >> e.frequency; continue; }
		if (cmd.compare("angle_in")    == 0) { ss >> e.angle_in;  continue; }
		if (cmd.compare("angle_out")   == 0) { ss >> e.angle_out; continue; }
		if (cmd.compare("t_in")      == 0) { for(auto i=0; i<6; ++i) ss >> e.t_in[i];  continue; }
		if (cmd.compare("t_out")     == 0) { for(auto i=0; i<6; ++i) ss >> e.t_out[i]; continue; }
		if (cmd.compare("rx|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[0*6+i]; continue; }
		if (cmd.compare("px|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[1*6+i]; continue; }
		if (cmd.compare("ry|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[2*6+i]; continue; }
		if (cmd.compare("py|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[3*6+i]; continue; }
		if (cmd.compare("de|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[4*6+i]; continue; }
		if (cmd.compare("dl|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[5*6+i]; continue; }
		if (cmd.compare("rx|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[0*6+i]; continue; }
		if (cmd.compare("px|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[1*6+i]; continue; }
		if (cmd.compare("ry|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[2*6+i]; continue; }
		if (cmd.compare("py|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[3*6+i]; continue; }
		if (cmd.compare("de|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[4*6+i]; continue; }
		if (cmd.compare("dl|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[5*6+i]; continue; }
		if (cmd.compare("pass_method") == 0) {
			std::string pass_method; ss >> pass_method;
			bool found_pm = false;
			for(unsigned int i = 0; i<((unsigned int)PassMethod::pm_nr_pms); ++i) {
				if (pass_method.compare(pm_dict[i]) == 0) {
					e.pass_method = i;
                    if (pass_method.compare("kicktable_pass") == 0) {
                        Status::type status = add_kicktable(e.fam_name + ".txt", accelerator.kicktables, e.kicktable);
				        if (status != Status::success) {
                            return status;
                        } else {
                        }
                    }
					found_pm = true;
					break;
				}
			}
			if (found_pm) continue;
			return Status::passmethod_not_defined;
		}
		if (cmd.compare("polynom_a") == 0) {
			std::vector<unsigned int> order;
			std::vector<double> multipole;
			unsigned int size = 0;
			while (not ss.eof()) {
				unsigned int o; double m; ss >> o >> m;
				if (ss.eof()) break;
				order.push_back(o); multipole.push_back(m);
				if (o+1 > size) size = o+1;
			}
			if (size > 0) {
				e.polynom_a.resize(size, 0);
				for(unsigned int i=0; i<order.size(); ++i) e.polynom_a[order[i]] = multipole[i];
			}
			synchronize_polynomials(e);
			continue;
		}
		if (cmd.compare("polynom_b") == 0) {
			std::vector<unsigned int> order;
			std::vector<double> multipole;
			unsigned int size = 0;
			while (not ss.eof()) {
				unsigned int o; double m; ss >> o >> m;
				if (ss.eof()) break;
				order.push_back(o); multipole.push_back(m);
				if (o+1 > size) size = o+1;
			}
			if (size > 0) {
				e.polynom_b.resize(size, 0);
				for(unsigned int i=0; i<order.size(); ++i) e.polynom_b[order[i]] = multipole[i];
			}
			synchronize_polynomials(e);
			continue;
		}
		if (line.size()<2) continue;
		return Status::flat_file_error;
		//std::cout << line_count << ": " << line << std::endl;
	}
	fp.close();

	if (e.fam_name.compare("") != 0) {
		accelerator.lattice.push_back(e);
	}
	return Status::success;

}

Status::type read_flat_file_trackcpp2(const std::string& filename, Accelerator& accelerator) {

	std::ifstream fp(filename);
	if (fp.fail()) return Status::file_not_found;

	accelerator.lattice.clear();

	Element e;
	while (not fp.eof()) {
		std::string cmd, line;
		unsigned int index = 0;
		fp >> cmd;
		//std::cout << cmd << std::endl;
		if (cmd[0] == '#') {
			std::getline(fp, line);
			continue;
		}
		if (cmd.compare("fam_name") == 0) {
			if (e.fam_name.compare("") != 0) {
				accelerator.lattice.push_back(e);
				e = Element();
			} 
			fp >> e.fam_name;
			continue;
		}
		if (cmd.compare("index")       == 0) { fp >> index;       continue; }
		if (cmd.compare("length")      == 0) { fp >> e.length;    continue; }
		if (cmd.compare("hmax")        == 0) { fp >> e.hmax;      continue; }
		if (cmd.compare("vmax")        == 0) { fp >> e.vmax;      continue; }
		if (cmd.compare("hkick")       == 0) { fp >> e.hkick;     continue; }
		if (cmd.compare("vkick")       == 0) { fp >> e.vkick;     continue; }
		if (cmd.compare("nr_steps")    == 0) { fp >> e.nr_steps;  continue; }
		if (cmd.compare("angle")       == 0) { fp >> e.angle;     continue; }
		if (cmd.compare("gap")         == 0) { fp >> e.gap;       continue; }
		if (cmd.compare("fint_in")     == 0) { fp >> e.fint_in;   continue; }
		if (cmd.compare("fint_out")    == 0) { fp >> e.fint_out;  continue; }
		if (cmd.compare("angle_in")    == 0) { fp >> e.angle_in;  continue; }
		if (cmd.compare("angle_out")   == 0) { fp >> e.angle_out; continue; }
		if (cmd.compare("t_in")      == 0) { for(auto i=0; i<6; ++i) fp >> e.t_in[i];  continue; }
		if (cmd.compare("t_out")     == 0) { for(auto i=0; i<6; ++i) fp >> e.t_out[i]; continue; }
		if (cmd.compare("rx|r_in")   == 0) { for(auto i=0; i<6; ++i) fp >> e.r_in[0*6+i]; continue; }
		if (cmd.compare("px|r_in")   == 0) { for(auto i=0; i<6; ++i) fp >> e.r_in[1*6+i]; continue; }
		if (cmd.compare("ry|r_in")   == 0) { for(auto i=0; i<6; ++i) fp >> e.r_in[2*6+i]; continue; }
		if (cmd.compare("py|r_in")   == 0) { for(auto i=0; i<6; ++i) fp >> e.r_in[3*6+i]; continue; }
		if (cmd.compare("de|r_in")   == 0) { for(auto i=0; i<6; ++i) fp >> e.r_in[4*6+i]; continue; }
		if (cmd.compare("dl|r_in")   == 0) { for(auto i=0; i<6; ++i) fp >> e.r_in[5*6+i]; continue; }
		if (cmd.compare("rx|r_out")  == 0) { for(auto i=0; i<6; ++i) fp >> e.r_out[0*6+i]; continue; }
		if (cmd.compare("px|r_out")  == 0) { for(auto i=0; i<6; ++i) fp >> e.r_out[1*6+i]; continue; }
		if (cmd.compare("ry|r_out")  == 0) { for(auto i=0; i<6; ++i) fp >> e.r_out[2*6+i]; continue; }
		if (cmd.compare("py|r_out")  == 0) { for(auto i=0; i<6; ++i) fp >> e.r_out[3*6+i]; continue; }
		if (cmd.compare("de|r_out")  == 0) { for(auto i=0; i<6; ++i) fp >> e.r_out[4*6+i]; continue; }
		if (cmd.compare("dl|r_out")  == 0) { for(auto i=0; i<6; ++i) fp >> e.r_out[5*6+i]; continue; }
		if (cmd.compare("pass_method") == 0) {
			std::string pass_method; fp >> pass_method;
			bool found_pm = false;
			for(unsigned int i = 0; i<((unsigned int)PassMethod::pm_nr_pms); ++i) {
				if (pass_method.compare(pm_dict[i]) == 0) {
					e.pass_method = i;
					found_pm = true;
					break;
				}
			}
			if (found_pm) continue;
			return Status::passmethod_not_defined;
		}
		if (cmd.compare("polynom_a") == 0) {
			std::getline(fp, line);
			std::istringstream ss(line);
			std::vector<unsigned int> order;
			std::vector<double> multipole;
			unsigned int size = 0;
			while (not ss.eof()) {
				unsigned int o; double m; ss >> o >> m;
				if (ss.eof()) break;
				order.push_back(o); multipole.push_back(m);
				if (o+1 > size) size = o+1;
			}
			if (size > 0) {
				e.polynom_a.resize(size, 0);
				for(unsigned int i=0; i<order.size(); ++i) e.polynom_a[order[i]] = multipole[i];
			}
			synchronize_polynomials(e);
			continue;
		}
		if (cmd.compare("polynom_b") == 0) {
			std::getline(fp, line);
			std::istringstream ss(line);
			std::vector<unsigned int> order;
			std::vector<double> multipole;
			unsigned int size = 0;
			while (not ss.eof()) {
				unsigned int o; double m; ss >> o >> m;
				if (ss.eof()) break;
				order.push_back(o); multipole.push_back(m);
				if (o+1 > size) size = o+1;
			}
			if (size > 0) {
				e.polynom_b.resize(size, 0);
				for(unsigned int i=0; i<order.size(); ++i) e.polynom_b[order[i]] = multipole[i];
			}
			synchronize_polynomials(e);
			continue;
		}
	}
	accelerator.lattice.push_back(e);
	fp.close();
	
	return Status::success;

}

void synchronize_polynomials(Element& e) {

	unsigned int size = (e.polynom_a.size() > e.polynom_b.size()) ? e.polynom_a.size() : e.polynom_b.size();
	e.polynom_a.resize(size, 0);
	e.polynom_b.resize(size, 0);

}

static void read_polynomials(std::ifstream& fp, Element& e);

Status::type read_flat_file_tracy(const std::string& filename, Accelerator& accelerator) {

	std::ifstream fp(filename);
	if (fp.fail()) return Status::file_not_found;

	int Fnum, Knum, idx, type, method;
	accelerator.lattice.clear();
	while (not fp.eof()) {

		Element e;

		// reads header
		fp >> e.fam_name >> Fnum >> Knum >> idx; if (fp.eof()) break;
		if (e.fam_name == "prtmfile:") return Status::flat_file_error;
		fp >> type >> method >> e.nr_steps;
		if (e.nr_steps < 1) e.nr_steps = 1;
		fp >> e.hmax >> e.hmax >> e.vmax >> e.vmax;

		// tracy starts with "begin" zero-length drift element...
		if (e.fam_name == "begin") {
			fp >> e.length;
			continue;
		}

		// element-type specifics
		switch (type) {
			case FlatFileType::marker:
			{
				e.pass_method = PassMethod::pm_identity_pass;
			}; break;
			case FlatFileType::drift:
			{
				e.pass_method = PassMethod::pm_drift_pass;
				fp >> e.length;
			}; break;
			case FlatFileType::corrector:
			{
				e.pass_method = PassMethod::pm_corrector_pass;
				double tmpdbl; fp >> tmpdbl >> tmpdbl >> tmpdbl;
				int    tmpint; fp >> tmpint >> tmpint;
				fp >> tmpint >> e.hkick >> e.vkick;
				e.hkick = - e.hkick; // AT idiosyncrasies...
			}; break;
			case FlatFileType::cavity:
			{
				e.pass_method = PassMethod::pm_cavity_pass;
				int hnumber; double energy;
				fp >> e.voltage >> e.frequency >> hnumber >> energy;
				e.voltage *= energy; e.frequency *= light_speed / (2*M_PI);
				accelerator.harmonic_number = hnumber;
				accelerator.energy = energy;
			}; break;
			case FlatFileType::mpole:
			{
				double PdTPar, PdTerr;
				fp >> e.t_out[0] >> e.t_out[2] >> PdTPar >> PdTerr;
				fp >> e.length >> e.angle >> e.angle_in >> e.angle_out >> e.gap;
				e.angle *= e.length; e.angle_in *= M_PI/180.0; e.angle_out *= M_PI/180.0;
				if (e.angle != 0)
					e.pass_method = PassMethod::pm_bnd_mpole_symplectic4_pass;
				else
					e.pass_method = PassMethod::pm_str_mpole_symplectic4_pass;
				read_polynomials(fp, e);
				e.t_in[0] = -e.t_out[0]; e.t_in[2] = -e.t_out[2];
				double ang = M_PI*(PdTPar+PdTerr)/180.0;
				double C = cos(ang);
				double S = sin(ang);
				e.r_in [0*6+0] =  C; e.r_in [0*6+2] =  S;
				e.r_in [2*6+0] = -S; e.r_in [2*6+2] =  C;
				e.r_in [1*6+1] =  C; e.r_in [1*6+3] =  S;
				e.r_in [3*6+1] = -S; e.r_in [3*6+3] =  C;
				e.r_out[0*6+0] =  C; e.r_out[0*6+2] = -S;
				e.r_out[2*6+0] =  S; e.r_out[2*6+2] =  C;
				e.r_out[1*6+1] =  C; e.r_out[1*6+3] = -S;
				e.r_out[3*6+1] =  S; e.r_out[3*6+3] =  C;
			}; break;
			case FlatFileType::kicktable:
			{
				e.pass_method = PassMethod::pm_kicktable_pass;
				double tmpdbl; std::string filename;
				fp >> tmpdbl >> tmpdbl >> filename;
				Status::type status = add_kicktable(filename, accelerator.kicktables, e.kicktable);
				if (status == Status::success) {
					e.length = e.kicktable->length;
					//std::cout << accelerator.lattice.size() << " " << e.fam_name << ": " << e.kicktable << " " << e.kicktable->x_nrpts << std::endl;
				} else return status;

			}; break;
			default:
				break;

		}


		accelerator.lattice.push_back(e); idx++;

	};

	return Status::success;

}

static void read_polynomials(std::ifstream& fp, Element& e) {
	unsigned int nr_monomials, n_design, order;
	e.polynom_a = std::vector<double>(Element::default_polynom);
	e.polynom_b = std::vector<double>(Element::default_polynom);
	fp >> nr_monomials >> n_design;
	for(unsigned int i=0; i<nr_monomials; ++i) {
		fp >> order;
		if (order > e.polynom_b.size()) {
			e.polynom_a.resize(order,0);
			e.polynom_b.resize(order,0);
		}
		fp >> e.polynom_b[order-1] >> e.polynom_a[order-1];
	}
}
