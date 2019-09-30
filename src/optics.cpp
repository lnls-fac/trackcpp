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

#include <trackcpp/optics.h>
#include <trackcpp/linalg.h>

double get_magnetic_rigidity(const double energy) {
    double gamma = (energy/1e6) / (electron_rest_energy_MeV);
    double beta  = sqrt(1 - 1/(gamma*gamma));
    double b_rho = beta * energy / light_speed; // [T.m]
    return b_rho;
}


// std::ostream& operator<< (std::ostream &out, const Twiss& tw) {
//   char buffer[255];
//   std::cout << "aqui" << std::endl;
//
//   out << "rx    :"; sprintf(buffer,  "%+.4e, ", tw.co.rx); out << buffer;
//   out << "ry    :"; sprintf(buffer,  "%+.4e, ", tw.co.ry); out << buffer;
//   out << std::endl;
//   out << "px    :"; sprintf(buffer,  "%+.4e, ", tw.co.px); out << buffer;
//   out << "py    :"; sprintf(buffer,  "%+.4e, ", tw.co.py); out << buffer;
//   out << std::endl;
//   out << "de    :"; sprintf(buffer,  "%+.4e, ", tw.co.de); out << buffer;
//   out << "dl    :"; sprintf(buffer,  "%+.4e, ", tw.co.dl); out << buffer;
//   out << std::endl;
//   out << "etax  :"; sprintf(buffer,  "%+.4e, ", tw.etax[0]); out << buffer;
//   out << "etay  :"; sprintf(buffer,  "%+.4e, ", tw.etay[0]); out << buffer;
//   out << std::endl;
//   out << "etalx :"; sprintf(buffer,  "%+.4e, ", tw.etax[0]); out << buffer;
//   out << "etaly :"; sprintf(buffer,  "%+.4e, ", tw.etay[0]); out << buffer;
//   out << std::endl;
//   out << "mux   :"; sprintf(buffer,  "%+.4e, ", tw.mux); out << buffer;
//   out << "muy   :"; sprintf(buffer,  "%+.4e, ", tw.muy); out << buffer;
//   out << std::endl;
//   out << "betax :"; sprintf(buffer,  "%+.4e, ", tw.betax); out << buffer;
//   out << "betay :"; sprintf(buffer,  "%+.4e, ", tw.betay); out << buffer;
//   out << std::endl;
//   out << "alphax:"; sprintf(buffer,  "%+.4e, ", tw.alphax); out << buffer;
//   out << "alphay:"; sprintf(buffer,  "%+.4e  ", tw.alphay); out << buffer;
//   return out;
// }

//#define TIMEIT
#ifdef TIMEIT
#include <chrono>
#endif

Status::type calc_twiss(const Accelerator& accelerator,
                        const Pos<double>& fixed_point,
                        Matrix& m66,
                        std::vector<Twiss>& twiss,
                        Twiss twiss0,
                        bool closed_flag) {

#ifdef TIMEIT
  auto start = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();
  auto diff = end - start;
#endif

#ifdef TIMEIT
  start = std::chrono::steady_clock::now();
#endif

  // propagates fixed point to find closed orbit
  Pos<double> fp = fixed_point;
  std::vector<Pos<double>> closed_orbit;
  Plane::type lost_plane;
  unsigned int element_offset = 0;
  Status::type status = track_linepass(accelerator, fp, closed_orbit, element_offset, lost_plane, true);
  if (status != Status::success) return status;


#ifdef TIMEIT
  end = std::chrono::steady_clock::now(); diff = end - start;
  std::cout << "calc_twiss::close_orbit: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
#endif

#ifdef TIMEIT
  start = std::chrono::steady_clock::now();
#endif

  std::vector<Matrix> atm;
  Pos<double> v0;
  status = track_findm66(accelerator, closed_orbit[0], atm, m66, v0);
  if (status != Status::success) return status;

#ifdef TIMEIT
  end = std::chrono::steady_clock::now(); diff = end - start;
  std::cout << "calc_twiss::findm66: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
#endif

#ifdef TIMEIT
  start = std::chrono::steady_clock::now();
#endif

  const double dpp = 1e-8;
  Pos<double> fpp = fixed_point;
  fpp.de += dpp;
  std::vector<Pos<double>> codp;

  // initial twiss parameters
  twiss.clear(); twiss.reserve(atm.size());
  if (twiss0.isundef()) { // case of periodic solution
    twiss0.spos = 0.0;
    // --- beta functions
    double sin_mux = sgn(m66[0][1]) * std::sqrt(-m66[0][1]*m66[1][0]-pow(m66[0][0]-m66[1][1],2)/4.0);
    double sin_muy = sgn(m66[2][3]) * std::sqrt(-m66[2][3]*m66[3][2]-pow(m66[2][2]-m66[3][3],2)/4.0);
    twiss0.alphax = (m66[0][0]-m66[1][1])/2.0/sin_mux;
    twiss0.alphay = (m66[2][2]-m66[3][3])/2.0/sin_muy;
    twiss0.betax  =  m66[0][1]/sin_mux;
    twiss0.betay  =  m66[2][3]/sin_muy;
    // --- closed orbit
    twiss0.co = closed_orbit[0];

    // // --- dispersion function based on eta = (1 - M)^(-1) D
    // Vector Dx({m66[0][4], m66[1][4]});
    // Vector Dy({m66[2][4], m66[3][4]});
    // Matrix eye2({{1,0},{0,1}});
    // Matrix mx; m66.getM(mx, 2, 2, 0, 0); mx.linear_combination(1.0,eye2,-1.0,mx); mx.inverse();
    // Matrix my; m66.getM(my, 2, 2, 2, 2); my.linear_combination(1.0,eye2,-1.0,my); my.inverse();
    // twiss0.etax.multiplication(mx, Dx);
    // twiss0.etay.multiplication(my, Dy);

    // Dispersion Function based on tracking:
    Status::type status = track_findorbit4(accelerator, codp, fpp);
    if (status != Status::success) return status;
    twiss0.etax[0] = (codp[0].rx - closed_orbit[0].rx) / dpp;
    twiss0.etax[1] = (codp[0].px - closed_orbit[0].px) / dpp;
    twiss0.etay[0] = (codp[0].ry - closed_orbit[0].ry) / dpp;
    twiss0.etay[1] = (codp[0].py - closed_orbit[0].py) / dpp;
  } else {
    fpp.rx += twiss0.etax[0] * dpp;
    fpp.px += twiss0.etax[1] * dpp;
    fpp.ry += twiss0.etay[0] * dpp;
    fpp.py += twiss0.etay[1] * dpp;
    Status::type status = track_linepass(accelerator, fpp, codp, element_offset, lost_plane, true);
    if (status != Status::success) return status;
  }

  if (not closed_flag) {
    closed_orbit.pop_back();
    atm.pop_back();
    codp.pop_back();
  }

  twiss.push_back(twiss0);

#ifdef TIMEIT
  end = std::chrono::steady_clock::now(); diff = end - start;
  std::cout << "calc_twiss::twiss0: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
#endif


#ifdef TIMEIT
  start = std::chrono::steady_clock::now();
#endif

  for(unsigned int i=1; i<atm.size(); ++i) {

    const Matrix& tm = atm[i];
    Twiss tw;
    tw.spos = twiss.back().spos + accelerator.lattice[i-1].length;
    // --- beta functions
    tw.betax = (std::pow(tm[0][0] * twiss0.betax - tm[0][1] * twiss0.alphax, 2) + pow(tm[0][1], 2)) / twiss0.betax;
    tw.betay = (std::pow(tm[2][2] * twiss0.betay - tm[2][3] * twiss0.alphay, 2) + pow(tm[2][3], 2)) / twiss0.betay;
    tw.alphax = -((tm[0][0] * twiss0.betax - tm[0][1] * twiss0.alphax) * (tm[1][0] * twiss0.betax - tm[1][1] * twiss0.alphax) + tm[0][1]*tm[1][1])/twiss0.betax;
    tw.alphay = -((tm[2][2] * twiss0.betay - tm[2][3] * twiss0.alphay) * (tm[3][2] * twiss0.betay - tm[3][3] * twiss0.alphay) + tm[2][3]*tm[3][3])/twiss0.betay;
    tw.mux = std::atan(tm[0][1]/(tm[0][0] * twiss0.betax - tm[0][1] * twiss0.alphax));
    tw.muy = std::atan(tm[2][3]/(tm[2][2] * twiss0.betay - tm[2][3] * twiss0.alphay));

    // --- closed orbit
    tw.co = closed_orbit[i];

    // // --- dispersion function based on propagation
    // Matrix t1(atm[i-1]); t1.inverse();
    // Matrix T; T.multiplication(atm[i],t1);
    // Vector Dx({T[0][4], T[1][4]});
    // Vector Dy({T[2][4], T[3][4]});
    // Matrix mx; T.getMx(mx);
    // Matrix my; T.getMy(my);
    // tw.etax.multiplication(mx, twiss[i-1].etax);
    // tw.etay.multiplication(my, twiss[i-1].etay);
    // tw.etax = Dx + tw.etax.multiplication(mx, twiss[i-1].etax);
    // tw.etay = Dy + tw.etay.multiplication(my, twiss[i-1].etay);

    // Dispersion Function based on tracking
    tw.etax[0] = (codp[i].rx - closed_orbit[i].rx) / dpp;
    tw.etax[1] = (codp[i].px - closed_orbit[i].px) / dpp;
    tw.etay[0] = (codp[i].ry - closed_orbit[i].ry) / dpp;
    tw.etay[1] = (codp[i].py - closed_orbit[i].py) / dpp;

    twiss.push_back(tw);
  }

  // unwraps betatron phases
  std::vector<double> jumpx(twiss.size(),0);
  std::vector<double> jumpy(twiss.size(),0);
  for(unsigned int i=1; i<twiss.size(); ++i) {
    if (twiss[i].mux < twiss[i-1].mux) jumpx[i] = jumpx[i-1] + M_PI; else jumpx[i] = jumpx[i-1];
    if (twiss[i].muy < twiss[i-1].muy) jumpy[i] = jumpy[i-1] + M_PI; else jumpy[i] = jumpy[i-1];
    twiss[i-1].mux += jumpx[i-1];
    twiss[i-1].muy += jumpy[i-1];
  }
  twiss.back().mux += jumpx.back();
  twiss.back().muy += jumpy.back();

#ifdef TIMEIT
  end = std::chrono::steady_clock::now(); diff = end - start;
  std::cout << "calc_twiss::propagation: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
#endif

  return Status::success;
}
