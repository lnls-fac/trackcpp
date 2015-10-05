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

#include <trackcpp/optics.h>
#include <trackcpp/linalg.h>

double get_magnetic_rigidity(const double energy) {
    double gamma = (energy/1e6) / (electron_rest_energy_MeV);
    double beta  = sqrt(1 - 1/(gamma*gamma));
    double b_rho = beta * energy / light_speed; // [T.m]
    return b_rho;
}

//#include <chrono>
Status::type calc_twiss(const Accelerator& accelerator, const Pos<double>& fixed_point, Matrix& m66, std::vector<Twiss>& twiss) {

  std::vector<Matrix> tmv;
  Status::type status;

  // auto start = std::chrono::steady_clock::now();
  // auto end = std::chrono::steady_clock::now();
  // auto diff = end - start;

  // start = std::chrono::steady_clock::now();
  Pos<double> fp = fixed_point;
  std::vector<Pos<double>> closed_orbit;
  Plane::type lost_plane;
  unsigned int element_offset = 0;
  track_linepass(accelerator, fp, closed_orbit, element_offset, lost_plane, true);
  closed_orbit.pop_back();
  // end = std::chrono::steady_clock::now();
  // diff = end - start;
  // std::cout << "closed_orbit: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

  // start = std::chrono::steady_clock::now();
  status = track_findm66 (accelerator, closed_orbit, tmv, m66);
  if (status != Status::success) return status;
  m66 = tmv.back();
  // end = std::chrono::steady_clock::now();
  // diff = end - start;
  // std::cout << "findm66: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;


  // start = std::chrono::steady_clock::now();
  Twiss tw0;
  //tw0.etax.
  double sin_mux = sgn(m66[0][1]) * std::sqrt(-m66[0][1]*m66[1][0]-pow(m66[0][0]-m66[1][1],2)/4.0);
  double sin_muy = sgn(m66[2][3]) * std::sqrt(-m66[2][3]*m66[3][2]-pow(m66[2][2]-m66[3][3],2)/4.0);
  tw0.alphax = (m66[0][0]-m66[1][1])/2.0/sin_mux;
  tw0.alphay = (m66[2][2]-m66[3][3])/2.0/sin_muy;
  tw0.betax  =  m66[0][1]/sin_mux;
  tw0.betay  =  m66[2][3]/sin_muy;

  // dispersion function based on eta = (1 - M)^(-1) D
  Vector2 Dx = {m66[0][4], m66[1][4]};
  Vector2 Dy = {m66[2][4], m66[3][4]};
  Matrix2 eye; matrix2_eye(eye);
  Matrix2 mx, eye_mx; getmx(m66, mx); matrix2_lc(eye_mx, 1, eye, -1, mx);
  Matrix2 my, eye_my; getmy(m66, my); matrix2_lc(eye_my, 1, eye, -1, my);
  Vector2 etax; linalg_solve2(etax, eye_mx, Dx);
  Vector2 etay; linalg_solve2(etay, eye_my, Dy);
  tw0.etax = Pos<double>(etax[0],etax[1],0,0,1,0);
  tw0.etay = Pos<double>(0,0,etay[0],etay[1],1,0);

  std::cout << tw0.etax << std::endl;
  //std::vector<Pos<double>> eye_mx; matrix6_lc(eye_mx, 1, eye, -1, m66);
  //linalg_solve2(const std::vector<Pos<double> >& M, const Pos<double>& B) {

  //linalg_solve2(const std::vector<Pos<double> >& M, const Pos<double>& B) {


      // t.etax = _np.linalg.solve(_np.eye(2,2) - mx, Dx)
      // t.etay = _np.linalg.solve(_np.eye(2,2) - my, Dy)

  twiss.push_back(tw0);


  for(unsigned int i=1; i<tmv.size(); ++i) {
    const Matrix& tm = tmv[i];
    Twiss tw;
    tw.betax = (std::pow(tm[0][0] * tw0.betax - tm[0][1] * tw0.alphax, 2) + pow(tm[0][1], 2)) / tw0.betax;
    tw.betay = (std::pow(tm[2][2] * tw0.betay - tm[2][3] * tw0.alphay, 2) + pow(tm[2][3], 2)) / tw0.betay;
    tw.alphax = -((tm[0][0] * tw0.betax - tm[0][1] * tw0.alphax) * (tm[1][0] * tw0.betax - tm[1][1] * tw0.alphax) + tm[0][1]*tm[1][1])/tw0.betax;
    tw.alphay = -((tm[2][2] * tw0.betay - tm[2][3] * tw0.alphay) * (tm[3][2] * tw0.betay - tm[3][3] * tw0.alphay) + tm[2][3]*tm[3][3])/tw0.betay;
    tw.mux = std::atan(tm[0][1]/(tm[0][0] * tw0.betax - tm[0][1] * tw0.alphax));
    tw.muy = std::atan(tm[2][3]/(tm[2][2] * tw0.betay - tm[2][3] * tw0.alphay));
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

  // end = std::chrono::steady_clock::now();
  // diff = end - start;
  // std::cout << "calc_twiss: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

  return Status::success;
}
