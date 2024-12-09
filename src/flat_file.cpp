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

#include <trackcpp/flat_file.h>
#include <trackcpp/elements.h>
#include <trackcpp/auxiliary.h>
#include <trackcpp/linalg.h>
#include <fstream>
#include <string>
#include <sstream>

static const int hw = 18; // header field width
static const int pw = 16; // parameter field width
static const int np = 17; // number precision
static bool read_boolean_string(std::istringstream& ss);
static int process_rad_property(std::istringstream& ss);
static std::string get_boolean_string(bool value);
static bool has_matrix66(const Matrix& r);
static bool has_polynom(const std::vector<double>& p);
static void write_6d_vector(std::ostream& fp, const std::string& label, const double* t);
static void write_6d_vector(std::ostream& fp, const std::string& label, const std::vector<double>& t);
static void write_polynom(std::ostream& fp, const std::string& label, const std::vector<double>& p);
static void synchronize_polynomials(Element& e);
static void read_polynomials(std::ifstream& fp, Element& e);
static void write_flat_file_trackcpp(std::ostream& fp, const Accelerator& accelerator);
static Status::type read_flat_file_trackcpp(std::istream&, Accelerator& accelerator);


// -- implementation of API --

Status::type read_flat_file(std::string& filename, Accelerator& accelerator, bool file_flag) {
  if (file_flag) {
    std::ifstream fp(filename.c_str());
    if (fp.fail()) return Status::file_not_found;
    read_flat_file_trackcpp(fp, accelerator);
    fp.close();
    return Status::success;
  } else {
    std::stringstream fp;
    fp.str(filename);
    read_flat_file_trackcpp(fp, accelerator);
    //filename = fp.str();
    return Status::success;
  }

}

Status::type write_flat_file(std::string& filename, const Accelerator& accelerator, bool file_flag) {

  if (file_flag) {
    std::ofstream fp(filename.c_str());
    if (fp.fail()) return Status::file_not_found;
    write_flat_file_trackcpp(fp, accelerator);
    fp.close();
    return Status::success;
  } else {
    std::stringstream fp;
    write_flat_file_trackcpp(fp, accelerator);
    filename = fp.str(); // very costly. no other way?!
    return Status::success;
  }

}

void write_flat_file_trackcpp(std::ostream& fp, const Accelerator& accelerator) {

  fp.setf(std::ios_base::left | std::ios_base::scientific | std::ios_base::uppercase);
  fp.precision(np);

  fp << std::setw(hw) << "% energy" << accelerator.energy << " eV\n";
  fp << std::setw(hw) << "% harmonic_number" << accelerator.harmonic_number << "\n";
  fp << std::setw(hw) << "% cavity_on" << get_boolean_string(accelerator.cavity_on) << "\n";
  fp << std::setw(hw) << "% radiation_on" << accelerator.radiation_on << "\n";
  fp << std::setw(hw) << "% vchamber_on" << get_boolean_string(accelerator.vchamber_on) << "\n";
  fp << std::setw(hw) << "% lattice_version" << accelerator.lattice_version << "\n";
  fp << '\n';

  fp.setf(std::ios_base::showpos);

  for (auto i=0; i<accelerator.lattice.size(); ++i) {
    const Element& e = accelerator.lattice[i];
    fp.unsetf(std::ios_base::showpos);
    fp << "### " << std::setw(4) << std::setfill('0') << std::right << i << std::setfill(' ') << std::left << " ###\n";
    fp.setf(std::ios_base::showpos);
    fp << std::setw(pw) << "fam_name" << e.fam_name << '\n';
    fp << std::setw(pw) << "length" << e.length << '\n';
    fp << std::setw(pw) << "pass_method" << pm_dict[e.pass_method] << '\n';
    if (pm_dict[e.pass_method].compare("kicktable_pass") == 0) {
      unsigned int idx = e.kicktable_idx;
      const Kicktable& ktable = kicktable_list[idx];
      fp << std::setw(pw) << "kicktable_fname" << ktable.filename << '\n';
    }
    if (e.nr_steps != 1) {
      fp.unsetf(std::ios_base::showpos);
      fp << std::setw(pw) << "nr_steps" << e.nr_steps << '\n';
      fp.setf(std::ios_base::showpos);
    }
    if (has_polynom(e.polynom_a)) write_polynom(fp, "polynom_a", e.polynom_a);
    if (has_polynom(e.polynom_b)) write_polynom(fp, "polynom_b", e.polynom_b);
    if (e.vchamber != VChamberShape::rectangle) { fp << std::setw(pw) << "vchamber" << e.vchamber << '\n'; }
    if (e.hmin != -DBL_MAX) { fp << std::setw(pw) << "hmin" << e.hmin << '\n'; }
    if (e.hmax != DBL_MAX) { fp << std::setw(pw) << "hmax" << e.hmax << '\n'; }
    if (e.vmin != -DBL_MAX) { fp << std::setw(pw) << "vmin" << e.vmin << '\n'; }
    if (e.vmax != DBL_MAX) { fp << std::setw(pw) << "vmax" << e.vmax << '\n'; }
    if (e.hkick != 0) { fp << std::setw(pw) << "hkick" << e.hkick << '\n'; }
    if (e.vkick != 0) { fp << std::setw(pw) << "vkick" << e.vkick << '\n'; }
    if (e.angle != 0) { fp << std::setw(pw) << "angle" << e.angle << '\n'; }
    if (e.gap != 0) { fp << std::setw(pw) << "gap" << e.gap << '\n'; }
    if (e.fint_in != 0) { fp << std::setw(pw) << "fint_in" << e.fint_in << '\n'; }
    if (e.fint_out != 0) { fp << std::setw(pw) << "fint_out" << e.fint_out << '\n'; }
    if (e.voltage != 0) { fp << std::setw(pw) << "voltage" << e.voltage << '\n'; }
    if (e.frequency != 0) { fp << std::setw(pw) << "frequency" << e.frequency << '\n'; }
    if (e.phase_lag != 0) { fp << std::setw(pw) << "phase_lag" << e.phase_lag << '\n'; }
    if (e.angle_in != 0) { fp << std::setw(pw) << "angle_in" << e.angle_in << '\n'; }
    if (e.angle_out != 0) { fp << std::setw(pw) << "angle_out" << e.angle_out << '\n'; }
    if (e.rescale_kicks != 1.0) { fp << std::setw(pw) << "rescale_kicks" << e.rescale_kicks << '\n'; }
    if (e.has_t_in) write_6d_vector(fp, "t_in", e.t_in);
    if (e.has_t_out) write_6d_vector(fp, "t_out", e.t_out);
    if (e.has_r_in) {
      write_6d_vector(fp, "rx|r_in", &e.r_in[6*0]);
      write_6d_vector(fp, "px|r_in", &e.r_in[6*1]);
      write_6d_vector(fp, "ry|r_in", &e.r_in[6*2]);
      write_6d_vector(fp, "py|r_in", &e.r_in[6*3]);
      write_6d_vector(fp, "de|r_in", &e.r_in[6*4]);
      write_6d_vector(fp, "dl|r_in", &e.r_in[6*5]);
    }
    if (e.has_r_out) {
      write_6d_vector(fp, "rx|r_out", &e.r_out[6*0]);
      write_6d_vector(fp, "px|r_out", &e.r_out[6*1]);
      write_6d_vector(fp, "ry|r_out", &e.r_out[6*2]);
      write_6d_vector(fp, "py|r_out", &e.r_out[6*3]);
      write_6d_vector(fp, "de|r_out", &e.r_out[6*4]);
      write_6d_vector(fp, "dl|r_out", &e.r_out[6*5]);
    }
    if (has_matrix66(e.matrix66)) {
      write_6d_vector(fp, "rx|matrix66", e.matrix66[0]);
      write_6d_vector(fp, "px|matrix66", e.matrix66[1]);
      write_6d_vector(fp, "ry|matrix66", e.matrix66[2]);
      write_6d_vector(fp, "py|matrix66", e.matrix66[3]);
      write_6d_vector(fp, "de|matrix66", e.matrix66[4]);
      write_6d_vector(fp, "dl|matrix66", e.matrix66[5]);
    }

    fp << '\n';
  }

}

Status::type read_flat_file_trackcpp(std::istream& fp, Accelerator& accelerator) {

  PassMethodsClass pass_methods;

  accelerator.lattice.clear();

  Element e;
  std::string line;
  bool found_hmin = false;
  bool found_vmin = false;
  unsigned int line_count = 0;
  while (std::getline(fp, line)) {
    //std::cout << line << std::endl;
    line_count++;
    // std::cout << line_count << ": " << line << std::endl;
    std::istringstream ss(line);
    std::string cmd;
    ss >> cmd;
    if (cmd[0] == '#') continue;
    if (cmd[0] == '%') {
      ss >> cmd;
      if (cmd.compare("energy") == 0) { ss >> accelerator.energy; continue; }
      if (cmd.compare("harmonic_number") == 0) { ss >> accelerator.harmonic_number; continue; }
      if (cmd.compare("cavity_on") == 0) { accelerator.cavity_on = read_boolean_string(ss); continue; }
      // radiation_on needs its own processing function due backward compatibility
      if (cmd.compare("radiation_on") == 0){ accelerator.radiation_on = process_rad_property(ss); continue; }
      if (cmd.compare("vchamber_on") == 0) { accelerator.vchamber_on = read_boolean_string(ss); continue; }
      if (cmd.compare("lattice_version") == 0) { ss >> accelerator.lattice_version; continue; }
      continue;
    }
    if (cmd.compare("fam_name") == 0) {
      if (e.fam_name.compare("") != 0) {
        accelerator.lattice.push_back(e);
        e = Element();
      }
      ss >> e.fam_name;
      continue;
    }
    if (cmd.compare("length")      == 0) { ss >> e.length;    continue; }
    if (cmd.compare("vchamber")    == 0) { ss >> e.vchamber; continue; }
    if (cmd.compare("hmin")        == 0) { ss >> e.hmin; found_hmin = true; continue; }
    if (cmd.compare("hmax")        == 0) {
      ss >> e.hmax;
      if (not found_hmin){ e.hmin = -e.hmax;}
      found_hmin = false;
      continue;
      }
    if (cmd.compare("vmin")        == 0) { ss >> e.vmin; found_vmin = true; continue; }
    if (cmd.compare("vmax")        == 0) {
      ss >> e.vmax;
      if (not found_vmin){ e.vmin = -e.vmax;}
      found_vmin = false;
      continue;
      }
    if (cmd.compare("hkick")       == 0) { ss >> e.hkick;     continue; }
    if (cmd.compare("vkick")       == 0) { ss >> e.vkick;     continue; }
    if (cmd.compare("nr_steps")    == 0) { ss >> e.nr_steps;  continue; }
    if (cmd.compare("angle")       == 0) { ss >> e.angle;     continue; }
    if (cmd.compare("gap")         == 0) { ss >> e.gap;       continue; }
    if (cmd.compare("fint_in")     == 0) { ss >> e.fint_in;   continue; }
    if (cmd.compare("fint_out")    == 0) { ss >> e.fint_out;  continue; }
    if (cmd.compare("voltage")     == 0) { ss >> e.voltage;   continue; }
    if (cmd.compare("frequency")   == 0) { ss >> e.frequency; continue; }
    if (cmd.compare("phase_lag")   == 0) { ss >> e.phase_lag; continue; }
    if (cmd.compare("angle_in")    == 0) { ss >> e.angle_in;  continue; }
    if (cmd.compare("angle_out")   == 0) { ss >> e.angle_out; continue; }
    if (cmd.compare("rescale_kicks")   == 0) { ss >> e.rescale_kicks; continue; }
    if (cmd.compare("t_in")      == 0) { for(auto i=0; i<6; ++i) ss >> e.t_in[i];  continue; }
    e.reflag_t_in();

    if (cmd.compare("t_out")     == 0) { for(auto i=0; i<6; ++i) ss >> e.t_out[i]; continue; }
    e.reflag_t_out();

    if (cmd.compare("rx|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[0*6+i]; continue; }
    if (cmd.compare("px|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[1*6+i]; continue; }
    if (cmd.compare("ry|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[2*6+i]; continue; }
    if (cmd.compare("py|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[3*6+i]; continue; }
    if (cmd.compare("de|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[4*6+i]; continue; }
    if (cmd.compare("dl|r_in")   == 0) { for(auto i=0; i<6; ++i) ss >> e.r_in[5*6+i]; continue; }
    e.reflag_r_in();

    if (cmd.compare("rx|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[0*6+i]; continue; }
    if (cmd.compare("px|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[1*6+i]; continue; }
    if (cmd.compare("ry|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[2*6+i]; continue; }
    if (cmd.compare("py|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[3*6+i]; continue; }
    if (cmd.compare("de|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[4*6+i]; continue; }
    if (cmd.compare("dl|r_out")  == 0) { for(auto i=0; i<6; ++i) ss >> e.r_out[5*6+i]; continue; }
    e.reflag_r_out();

    if (cmd.compare("rx|matrix66")  == 0) { for(auto i=0; i<6; ++i) ss >> e.matrix66[0][i]; continue; }
    if (cmd.compare("px|matrix66")  == 0) { for(auto i=0; i<6; ++i) ss >> e.matrix66[1][i]; continue; }
    if (cmd.compare("ry|matrix66")  == 0) { for(auto i=0; i<6; ++i) ss >> e.matrix66[2][i]; continue; }
    if (cmd.compare("py|matrix66")  == 0) { for(auto i=0; i<6; ++i) ss >> e.matrix66[3][i]; continue; }
    if (cmd.compare("de|matrix66")  == 0) { for(auto i=0; i<6; ++i) ss >> e.matrix66[4][i]; continue; }
    if (cmd.compare("dl|matrix66")  == 0) { for(auto i=0; i<6; ++i) ss >> e.matrix66[5][i]; continue; }
    if (cmd.compare("pass_method") == 0) {
      std::string pass_method; ss >> pass_method;
      bool found_pm = false;
      for(unsigned int i = 0; i<pass_methods.size(); ++i) {
        if (pass_method.compare(pass_methods[i]) == 0) {
          e.pass_method = i;
          found_pm = true;
          break;
        }
      }
      if (found_pm) continue;
      return Status::passmethod_not_defined;
    }
    if (cmd.compare("kicktable_fname") == 0) {
      std::string fname; ss >> fname;
      int idx = add_kicktable(fname);
      if (idx < 0) {
        return Status::file_not_found;
      } else {
        e.kicktable_idx = idx;
        continue;
      }
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
  }

  if (e.fam_name.compare("") != 0) {
    accelerator.lattice.push_back(e);
  }

  return Status::success;

}

static Status::type read_flat_file_tracy(const std::string& filename, Accelerator& accelerator) {

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
    fp >> e.hmin >> e.hmax >> e.vmin >> e.vmax;


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
        // NOTE: which one?
        // e.pass_method = PassMethod::pm_drift_g2l_pass;
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

        e.reflag_t_in();
        e.reflag_t_out();
        e.reflag_r_in();
        e.reflag_r_out();

      }; break;
      case FlatFileType::kicktable:
      {
        e.pass_method = PassMethod::pm_kickmap_pass;
        double tmpdbl; std::string filename;
        fp >> tmpdbl >> tmpdbl >> filename;
        int idx = add_kicktable(filename);
        if (idx >= 0) {
          e.kicktable_idx = idx;
          e.length = kicktable_list[e.kicktable_idx].length;
          //std::cout << accelerator.lattice.size() << " " << e.fam_name << ": " << e.kicktable << " " << e.kicktable->x_nrpts << std::endl;
        } else return Status::file_not_found;

      }; break;
      default:
        break;

    }

    accelerator.lattice.push_back(e); idx++;

  };

  return Status::success;

}


// -- static auxiliary functions --

static void synchronize_polynomials(Element& e) {

  unsigned int size = (e.polynom_a.size() > e.polynom_b.size()) ? e.polynom_a.size() : e.polynom_b.size();
  e.polynom_a.resize(size, 0);
  e.polynom_b.resize(size, 0);

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

static bool read_boolean_string(std::istringstream& ss) {
  std::string s;
  ss >> s;
  if (s.compare("true") == 0)
    return true;
  else
    return false;
}

static int process_rad_property(std::istringstream& ss){
  std::string s;
  ss >> s;
  if (std::isdigit(s[0])){
    return std::stoi(s);
  }else if(s.compare("true") == 0){
    return 1;
  }else{
    return 0;
  }
}

static std::string get_boolean_string(bool value) {
  if (value)
    return "true";
  else
    return "false";
}

static bool has_matrix66(const Matrix& m) {
  for (int i=0; i<6; ++i)
    for (int j=0; j<6; ++j)
      if ((i != j) & (m[i][j] != 0.0))
        return true;
      else if ((i == j) & (m[i][j] != 1.0))
        return true;
  return false;
}

static bool has_polynom(const std::vector<double>& p) {
  for (int i=0; i<p.size(); ++i)
    if (p[i] != 0)
      return true;

  return false;
}

static void write_6d_vector(std::ostream& fp, const std::string& label, const double* t) {
  fp << std::setw(pw) << label;
  for (int i=0; i<6; ++i)
    fp << t[i] << "  ";
  fp << '\n';
}

static void write_6d_vector(std::ostream& fp, const std::string& label, const std::vector<double>& t) {
  fp << std::setw(pw) << label;
  for (int i=0; i<6; ++i)
    fp << t[i] << "  ";
  fp << '\n';
}

static void write_polynom(std::ostream& fp, const std::string& label, const std::vector<double>& p) {
  fp << std::setw(pw) << label;
  for (int i=0; i<p.size(); ++i)
    if (p[i] != 0) {
      fp.unsetf(std::ios_base::showpos);
      fp << i << ' ';
      fp.setf(std::ios_base::showpos);
      fp << p[i] << ' ';
    }
  fp << '\n';
}
