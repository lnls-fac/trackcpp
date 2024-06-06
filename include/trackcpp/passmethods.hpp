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

// Trackcpp passmethods are based on SLAC Andrei Terebilo AT version 1.3
// <http://www.slac.stanford.edu/grp/ssrl/spear/at/>.

#ifndef _PASS_METHOD_AT_H
#define _PASS_METHOD_AT_H

// Obs: these are a c++ implementation of the AT passmethods.
//      apart from discrepancies between math library implementations and TWOPI,CGAMMA constants
//      these passmethods agree with AT passmethods up to machine 64-bit precision.

#include "passmethods.h"
#include "accelerator.h"
#include "elements.h"
#include "pos.h"
#include "auxiliary.h"
#include "tpsa.h"
#include "linalg.h"
#include <cmath>

template <typename T> inline T SQR(const T& X) { return X*X; }
template <typename T> inline T POW3(const T& X) { return X*X*X; }

template <typename T>
inline void drift(Pos<T>& pos, const double& length) {

  T pnorm = 1 / (1 + pos.de);
  T norml = length * pnorm;
  pos.rx += norml * pos.px;
  pos.ry += norml * pos.py;
  pos.dl += 0.5 * norml * pnorm * (pos.px*pos.px + pos.py*pos.py);
}

//template <typename T>
//inline void drift(Pos<T> &pos, const double& length) {
//  drift(pos, length);
//}

template <typename T>
inline void calcpolykick(const Pos<T> &pos, const std::vector<double>& polynom_a,
                         const std::vector<double>& polynom_b,
                         T& real_sum, T& imag_sum) {

  const int n = std::min(polynom_b.size(), polynom_a.size());
  if (n == 0) {
    real_sum = imag_sum = 0;
  } else {
    real_sum = polynom_b[n-1];
    imag_sum = polynom_a[n-1];
    for(int i=n-2;i>=0;--i) {
      T real_sum_tmp = real_sum * pos.rx - imag_sum * pos.ry + polynom_b[i];
      imag_sum = imag_sum * pos.rx + real_sum * pos.ry + polynom_a[i];
      real_sum = real_sum_tmp;
    }
  }
}

template <typename T>
void fastdrift(Pos<T> &pos, const T& norml) {

  T dx = norml * pos.px;
  T dy = norml * pos.py;
  pos.rx += dx;
  pos.ry += dy;
  pos.dl += 0.5 * norml * (pos.px*pos.px + pos.py*pos.py) / (1 + pos.de);
}

template <typename T>
T b2_perp(const T& bx, const T& by, const T& px, const T& py, const T& curv=1) {

  // Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity
  const T& curv2 = SQR(curv);
  const T& v_norm2_inv = (curv2 + SQR(px) + SQR(py));
  T&& b2p = SQR(by);
  b2p += SQR(bx);
  b2p *= curv2;
  b2p += SQR(bx*py - by*px);
  b2p /= v_norm2_inv;
  return b2p;
}

template <typename T>
Status::type kicktablethinkick(Pos<T>& pos, const int& kicktable_idx,
                               const double& brho, const int nr_steps, const double& rescale_kicks) {

  T hkick, vkick;
  Status::type status = kicktable_getkicks(kicktable_idx, pos.rx, pos.ry, hkick, vkick);
  pos.px += rescale_kicks * hkick / (brho * brho) / nr_steps;
  pos.py += rescale_kicks * vkick / (brho * brho) / nr_steps;
  if (status == Status::kicktable_out_of_range) {
    if (not isfinite(pos.px)) {
      pos.rx = nan("");
    }
    if (not isfinite(pos.py)) {
      pos.ry = nan("");
    }
  }
  return status;
}

template <typename T>
void matthinkick(Pos<T> &pos, const Matrix &m) {

  T rx = pos.rx, px = pos.px;
  T ry = pos.ry, py = pos.py;
  T de = pos.de, dl = pos.dl;

  pos.rx = m[0][0]*rx+m[0][1]*px+m[0][2]*ry+m[0][3]*py+m[0][4]*de+m[0][5]*dl;
  pos.px = m[1][0]*rx+m[1][1]*px+m[1][2]*ry+m[1][3]*py+m[1][4]*de+m[1][5]*dl;
  pos.ry = m[2][0]*rx+m[2][1]*px+m[2][2]*ry+m[2][3]*py+m[2][4]*de+m[2][5]*dl;
  pos.py = m[3][0]*rx+m[3][1]*px+m[3][2]*ry+m[3][3]*py+m[3][4]*de+m[3][5]*dl;
  pos.de = m[4][0]*rx+m[4][1]*px+m[4][2]*ry+m[4][3]*py+m[4][4]*de+m[4][5]*dl;
  pos.dl = m[5][0]*rx+m[5][1]*px+m[5][2]*ry+m[5][3]*py+m[5][4]*de+m[5][5]*dl;
}

template <typename T>
void strthinkick(Pos<T>& pos, const double& length,
                 const std::vector<double>& polynom_a,
                 const std::vector<double>& polynom_b,
                 const Accelerator& accelerator,
                 const double rad_const = 0,
                 const double qexcit_const = 0) {

  T real_sum, imag_sum;
  calcpolykick<T>(pos, polynom_a, polynom_b, real_sum, imag_sum);

  if (rad_const != 0) {
    T&& pnorm = 1 / (1 + pos.de);
    const T& rx = pos.rx;
    const T&  px = pos.px * pnorm;
    const T& ry = pos.ry;
    const T&  py = pos.py * pnorm;
    const T& b2p = b2_perp(imag_sum, real_sum, px, py);
    const T& delta_factor = SQR(1+pos.de);
    const T& dl_ds = (1+(px*px + py*py)/2);
    pos.de -= rad_const*delta_factor*b2p*dl_ds*length;

    if (qexcit_const != 0) {
      // quantum excitation kick
      const T& d = delta_factor * qexcit_const * sqrt(POW3(sqrt(b2p)) * dl_ds);
      pos.de += d * gen_random_number();
    }

    pnorm = (1 + pos.de);  // actually this is the inverse of pnorm
    pos.px = px * pnorm;
    pos.py = py * pnorm;
  }
  pos.px -= length * real_sum;
  pos.py += length * imag_sum;
}

template <typename T>
void bndthinkick(Pos<T>& pos, const double& length,
                 const std::vector<double>& polynom_a,
                 const std::vector<double>& polynom_b,
                 const double& irho,
                 const Accelerator& accelerator,
                 const double rad_const = 0,
                 const double qexcit_const = 0) {

  T real_sum, imag_sum;
  calcpolykick<T>(pos, polynom_a, polynom_b, real_sum, imag_sum);
  T de = pos.de;

  if (rad_const != 0) {
    T&& pnorm = 1 / (1 + pos.de);
    const T& rx = pos.rx;
    const T& px = pos.px * pnorm;
    const T& ry = pos.ry;
    const T& py = pos.py * pnorm;
    const T& curv = 1 + irho*rx;
    const T& b2p = b2_perp(imag_sum, real_sum+irho, px, py, curv);
    const T& delta_factor = SQR(1 + pos.de);
    const T& dl_ds = (curv + (px*px+py*py)/2);
    pos.de -= rad_const*delta_factor*b2p*dl_ds*length;

    if (qexcit_const != 0) {
      // quantum excitation kick
      const T& d = delta_factor * qexcit_const * sqrt(POW3(sqrt(b2p)) * dl_ds);
      pos.de += d * gen_random_number();
    }

    pnorm = (1 + pos.de);  // actually this is the inverse of pnorm
    pos.px = px * pnorm;
    pos.py = py * pnorm;
  }
  pos.px -= length * (real_sum - (de - pos.rx * irho) * irho);
  pos.py += length * imag_sum;
  pos.dl += length * irho * pos.rx;
}

template <typename T>
void edge_fringe(Pos<T>& pos, const double& inv_rho,
                 const double& edge_angle, const double& fint,
                const double& gap) {

  const T &rx = pos.rx, &ry = pos.ry, &de = pos.de;
  T       &px = pos.px, &py = pos.py;
  T fx      = inv_rho * std::tan(edge_angle)/(1 + de);
  T psi_bar = edge_angle - inv_rho * gap * fint *
    (1 + std::sin(edge_angle) * std::sin(edge_angle))
    / std::cos(edge_angle) / (1 + de);
  T fy      = inv_rho * tan(psi_bar) / (1 + de);
  px       += rx * fx;
  py       -= ry * fy;
}

template <typename T>
inline void translate_pos(Pos<T> &pos, const double* t) {

  pos.rx += t[0]; pos.px += t[1];
  pos.ry += t[2]; pos.py += t[3];
  pos.de += t[4]; pos.dl += t[5];
}

template <typename T>
inline void rotate_pos(Pos<T> &pos, const double* R) {

  const T rx0 = pos.rx, px0 = pos.px;
  const T ry0 = pos.ry, py0 = pos.py;
  const T de0 = pos.de, dl0 = pos.dl;
  pos.rx = R[0*6+0] * rx0 + R[0*6+1] * px0 + R[0*6+2] * ry0 + R[0*6+3] * py0 + R[0*6+4] * de0 + R[0*6+5] * dl0;
  pos.px = R[1*6+0] * rx0 + R[1*6+1] * px0 + R[1*6+2] * ry0 + R[1*6+3] * py0 + R[1*6+4] * de0 + R[1*6+5] * dl0;
  pos.ry = R[2*6+0] * rx0 + R[2*6+1] * px0 + R[2*6+2] * ry0 + R[2*6+3] * py0 + R[2*6+4] * de0 + R[2*6+5] * dl0;
  pos.py = R[3*6+0] * rx0 + R[3*6+1] * px0 + R[3*6+2] * ry0 + R[3*6+3] * py0 + R[3*6+4] * de0 + R[3*6+5] * dl0;
  pos.de = R[4*6+0] * rx0 + R[4*6+1] * px0 + R[4*6+2] * ry0 + R[4*6+3] * py0 + R[4*6+4] * de0 + R[4*6+5] * dl0;
  pos.dl = R[5*6+0] * rx0 + R[5*6+1] * px0 + R[5*6+2] * ry0 + R[5*6+3] * py0 + R[5*6+4] * de0 + R[5*6+5] * dl0;
}

template <typename T>
void global_2_local(Pos<T> &pos, const Element &elem) {
  if (elem.has_t_in) {
      translate_pos(pos, elem.t_in);
  }
  if (elem.has_r_in) {
      rotate_pos(pos, elem.r_in);
  }
}

template <typename T>
void local_2_global(Pos<T> &pos, const Element &elem) {
  if (elem.has_r_out) {
      rotate_pos(pos, elem.r_out);
  }
  if (elem.has_t_out) {
      translate_pos(pos, elem.t_out);
  }
}

template <typename T>
Status::type pm_identity_pass(Pos<T> &pos, const Element &elem,
                              const Accelerator& accelerator) {
  return Status::success;
}

template <typename T>
Status::type pm_drift_g2l_pass(Pos<T> &pos, const Element &elem,
                                const Accelerator& accelerator) {

  global_2_local(pos, elem);
  drift<T>(pos, elem.length);
  local_2_global(pos, elem);
  return Status::success;
}

template <typename T>
Status::type pm_drift_pass(Pos<T> &pos, const Element &elem,
                           const Accelerator& accelerator) {
  drift<T>(pos, elem.length);
  return Status::success;
}

template <typename T>
Status::type pm_str_mpole_symplectic4_pass(Pos<T> &pos, const Element &elem,
                                           const Accelerator& accelerator) {

  global_2_local(pos, elem);
  double sl = elem.length / float(elem.nr_steps);
  double l1 = sl * DRIFT1;
  double l2 = sl * DRIFT2;
  double k1 = sl * KICK1;
  double k2 = sl * KICK2;
  const std::vector<double> &polynom_a = elem.polynom_a;
  const std::vector<double> &polynom_b = elem.polynom_b;
  double rad_const = 0;
  double qexcit_const = 0; // quantum excitation scale factor

  if (accelerator.radiation_on){
    rad_const = CGAMMA*POW3(accelerator.energy/1e9)/(TWOPI); /*[m] M.Sands(4.1)*/
  }

  if (accelerator.radiation_on == RadiationState::full){
    qexcit_const = CQEXT*SQR(accelerator.energy)*sqrt(accelerator.energy*sl);
  }

  for(unsigned int i=0; i<elem.nr_steps; ++i) {
    drift<T>(pos, l1);
    strthinkick<T>(pos, k1, polynom_a, polynom_b, accelerator, rad_const, 0);
    drift<T>(pos, l2);
    strthinkick<T>(pos, k2, polynom_a, polynom_b, accelerator, rad_const, qexcit_const);
    drift<T>(pos, l2);
    strthinkick<T>(pos, k1, polynom_a, polynom_b, accelerator, rad_const, 0);
    drift<T>(pos, l1);
  }
  local_2_global(pos, elem);
  return Status::success;
}

template <typename T>
Status::type pm_bnd_mpole_symplectic4_pass(Pos<T> &pos, const Element &elem,
                                           const Accelerator& accelerator) {

  double sl = elem.length / float(elem.nr_steps);
  double l1 = sl * DRIFT1;
  double l2 = sl * DRIFT2;
  double k1 = sl * KICK1;
  double k2 = sl * KICK2;
  double irho = elem.angle / elem.length;
  const std::vector<double> &polynom_a = elem.polynom_a;
  const std::vector<double> &polynom_b = elem.polynom_b;
  double rad_const = 0;
  double qexcit_const = 0; // quantum excitation scale factor

  if (accelerator.radiation_on){
    rad_const = CGAMMA*POW3(accelerator.energy/1e9)/(TWOPI);
  }

  if (accelerator.radiation_on == RadiationState::full) {
    qexcit_const = CQEXT*SQR(accelerator.energy)*sqrt(accelerator.energy*sl);
  }

  global_2_local(pos, elem);
  edge_fringe(pos, irho, elem.angle_in, elem.fint_in, elem.gap);
  for(unsigned int i=0; i<elem.nr_steps; ++i) {
    drift<T>(pos, l1);
    bndthinkick<T>(pos, k1, polynom_a, polynom_b, irho, accelerator, rad_const, 0);
    drift<T>(pos, l2);
    bndthinkick<T>(pos, k2, polynom_a, polynom_b, irho, accelerator, rad_const, qexcit_const);
    drift<T>(pos, l2);
    bndthinkick<T>(pos, k1, polynom_a, polynom_b, irho, accelerator, rad_const, 0);
    drift<T>(pos, l1);
  }
  edge_fringe(pos, irho, elem.angle_out, elem.fint_out, elem.gap);
  local_2_global(pos, elem);

  return Status::success;
}

template <typename T>
Status::type pm_corrector_pass(Pos<T> &pos, const Element &elem,
                               const Accelerator& accelerator) {

  global_2_local(pos, elem);
  const double& xkick = elem.hkick;
  const double& ykick = elem.vkick;
  if (elem.length == 0) {
    T &px = pos.px, &py = pos.py;
    px += xkick;
    py += ykick;
  } else {
    T &rx = pos.rx, &px = pos.px;
    T &ry = pos.ry, &py = pos.py;
    T &de = pos.de, &dl = pos.dl;
    T pnorm   = 1 / (1 + de);
    T norml   = elem.length * pnorm;
    dl += norml * pnorm * 0.5 * (
        xkick * xkick/3.0 + ykick * ykick/3.0 +
        px*px + py*py +
        px * xkick + py * ykick
        );
    rx += norml * (px + 0.5 * xkick);
    px += xkick;
    ry += norml * (py + 0.5 * ykick);
    py += ykick;
  }
  local_2_global(pos, elem);
  return Status::success;
}

template <typename T>
Status::type pm_cavity_0_pass(Pos<T> &pos, const Element &elem,
                            const Accelerator& accelerator) {

  if (not accelerator.cavity_on) return pm_drift_pass(pos, elem, accelerator);

  global_2_local(pos, elem);
  double nv = elem.voltage / accelerator.energy;
  if (elem.length == 0) {
    T &de = pos.de, &dl = pos.dl;
    de +=  -nv * sin(TWOPI*elem.frequency * dl/ light_speed - elem.phase_lag);
    } else {
    T &rx = pos.rx, &px = pos.px;
    T &ry = pos.ry, &py = pos.py;
    T &de = pos.de, &dl = pos.dl;
    // drift half length
    T pnorm   = 1 / (1 + de);
    T norml   = (0.5 * elem.length) * pnorm;
    rx += norml * px;
    ry += norml * py;
    dl += 0.5 * norml * pnorm * (px*px + py*py);
    // longitudinal momentum kick
    de += -nv * sin(TWOPI*elem.frequency*dl/light_speed - elem.phase_lag);
    // drift half length
    pnorm   = 1.0 / (1.0 + de);
    norml   = (0.5 * elem.length) * pnorm;
    rx += norml * px;
    ry += norml * py;
    dl += 0.5 * norml * pnorm * (px*px + py*py);
  }
  local_2_global(pos, elem);
  return Status::success;
}

// template <typename T>
// Status::type pm_cavity_1comp_pass(Pos<T> &pos, const Element &elem,
//                             const Accelerator& accelerator) {

//   if (not accelerator.cavity_on) return pm_drift_pass(pos, elem, accelerator);

//   global_2_local(pos, elem);
//   double nv = elem.voltage / accelerator.energy;
//   double frf = elem.frequency;
//   double harmonic_number = accelerator.harmonic_number;
//   double L0 = accelerator.get_length();
//   double s = 0.0;
//   double accum_s = 0.0;
//   unsigned int last_cav_idx = 0;
//   for (unsigned int i = 0; i < accelerator.lattice.size(); i++){
//     auto elem2 = accelerator.lattice[i];
//     if (elem2 == elem){
//       accum_s += s;
//     }
//     if (elem2.frequency != 0.0){
//       s = 0.0;
//       last_cav_idx = i;
//     }
//     s += elem2.length;
//   }
//   double ddl = (light_speed*harmonic_number/frf - L0);
//   if (elem.length == 0) {
//     T &de = pos.de, &dl = pos.dl;
//     dl -= ddl*(accum_s/L0);
//     de +=  -nv * sin(TWOPI*frf *dl/ light_speed - elem.phase_lag);
//     if (elem == accelerator.lattice[last_cav_idx]){
//       dl -= ddl*(s/L0);
//     }
//     } else {
//     T &rx = pos.rx, &px = pos.px;
//     T &ry = pos.ry, &py = pos.py;
//     T &de = pos.de, &dl = pos.dl;
//     // drift half length
//     T pnorm   = 1 / (1 + de);
//     T norml   = (0.5 * elem.length) * pnorm;
//     rx += norml * px;
//     ry += norml * py;
//     dl += 0.5 * norml * pnorm * (px*px + py*py);
//     // longitudinal momentum kick
//     dl -= ddl*(accum_s/L0);
//     de += -nv * sin(TWOPI*frf *dl/ light_speed - elem.phase_lag);
//     if (elem == accelerator.lattice[last_cav_idx]){
//       dl -= ddl*(s/L0);
//     }
//     // drift half length
//     pnorm   = 1.0 / (1.0 + de);
//     norml   = (0.5 * elem.length) * pnorm;
//     rx += norml * px;
//     ry += norml * py;
//     dl += 0.5 * norml * pnorm * (px*px + py*py);
//   }
//   local_2_global(pos, elem);
//   return Status::success;
// }

template <typename T>
Status::type pm_cavity_1comp_pass(Pos<T> &pos, const Element &elem,
                            const Accelerator& accelerator) {

  if (not accelerator.cavity_on) return pm_drift_pass(pos, elem, accelerator);

  global_2_local(pos, elem);
  double nv = elem.voltage / accelerator.energy;
  double frf = elem.frequency;
  double harmonic_number = accelerator.harmonic_number;
  double L0 = accelerator.get_length();
  unsigned int last_cav_idx = 0;
  for (unsigned int i = 0; i < accelerator.lattice.size(); i++){
    auto elem2 = accelerator.lattice[i];
    if (elem2.frequency != 0.0){
      last_cav_idx = i;
    }
  }
  double ddl = (light_speed*harmonic_number/frf - L0);
  if (elem.length == 0) {
    T &de = pos.de, &dl = pos.dl;

    de +=  -nv * sin(TWOPI*frf *dl/ light_speed - elem.phase_lag);
    if (elem == accelerator.lattice[last_cav_idx]){
      dl -= ddl;
    }
    } else {
    T &rx = pos.rx, &px = pos.px;
    T &ry = pos.ry, &py = pos.py;
    T &de = pos.de, &dl = pos.dl;
    // drift half length
    T pnorm   = 1 / (1 + de);
    T norml   = (0.5 * elem.length) * pnorm;
    rx += norml * px;
    ry += norml * py;
    dl += 0.5 * norml * pnorm * (px*px + py*py);
    // longitudinal momentum kick

    de += -nv * sin(TWOPI*frf *dl/ light_speed - elem.phase_lag);
    if (elem == accelerator.lattice[last_cav_idx]){
      dl -= ddl;
    }
    // drift half length
    pnorm   = 1.0 / (1.0 + de);
    norml   = (0.5 * elem.length) * pnorm;
    rx += norml * px;
    ry += norml * py;
    dl += 0.5 * norml * pnorm * (px*px + py*py);
  }
  local_2_global(pos, elem);
  return Status::success;
}

template <typename T>
Status::type pm_cavity_1frac_pass(Pos<T> &pos, const Element &elem,
                            const Accelerator& accelerator) {

  if (not accelerator.cavity_on) return pm_drift_pass(pos, elem, accelerator);

  global_2_local(pos, elem);
  double nv = elem.voltage / accelerator.energy;
  double frf = elem.frequency;
  double harmonic_number = accelerator.harmonic_number;
  double L0 = accelerator.get_length();
  double time_aware_frac = accelerator.get_time_aware_frac();
  double ddl = (light_speed*harmonic_number/frf - L0)*time_aware_frac;
  if (elem.length == 0) {
    T &de = pos.de, &dl = pos.dl;
    dl -= ddl;
    de +=  -nv * sin(TWOPI*frf *dl/ light_speed - elem.phase_lag);
    } else {
    T &rx = pos.rx, &px = pos.px;
    T &ry = pos.ry, &py = pos.py;
    T &de = pos.de, &dl = pos.dl;
    // drift half length
    T pnorm   = 1 / (1 + de);
    T norml   = (0.5 * elem.length) * pnorm;
    rx += norml * px;
    ry += norml * py;
    dl += 0.5 * norml * pnorm * (px*px + py*py);
    // longitudinal momentum kick
    dl -= ddl;
    de += -nv * sin(TWOPI*frf *dl/ light_speed - elem.phase_lag);
    // drift half length
    pnorm   = 1.0 / (1.0 + de);
    norml   = (0.5 * elem.length) * pnorm;
    rx += norml * px;
    ry += norml * py;
    dl += 0.5 * norml * pnorm * (px*px + py*py);
  }
  local_2_global(pos, elem);
  return Status::success;
}

template <typename T>
Status::type pm_cavity_2_pass(Pos<T> &pos, const Element &elem,
                            const Accelerator& accelerator, const double nturn) {

  if (not accelerator.cavity_on) return pm_drift_pass(pos, elem, accelerator);

  global_2_local(pos, elem);
  double nv = elem.voltage / accelerator.energy;
  double frf = elem.frequency;
  double harmonic_number = accelerator.harmonic_number;
  double L0 = accelerator.get_length();
  // double time_aware_frac = accelerator.get_time_aware_frac();
  double ddl =  (light_speed*harmonic_number/frf - L0)* nturn;// * time_aware_frac ;
  if (elem.length == 0) {
    T &de = pos.de, &dl = pos.dl;
    de +=  -nv * sin(TWOPI*frf *(dl-ddl)/ light_speed - elem.phase_lag);
    } else {
    T &rx = pos.rx, &px = pos.px;
    T &ry = pos.ry, &py = pos.py;
    T &de = pos.de, &dl = pos.dl;
    // drift half length
    T pnorm   = 1 / (1 + de);
    T norml   = (0.5 * elem.length) * pnorm;
    rx += norml * px;
    ry += norml * py;
    dl += 0.5 * norml * pnorm * (px*px + py*py);
    // longitudinal momentum kick
    de += -nv * sin(TWOPI*frf *(dl-ddl)/ light_speed - elem.phase_lag);
    // drift half length
    pnorm   = 1.0 / (1.0 + de);
    norml   = (0.5 * elem.length) * pnorm;
    rx += norml * px;
    ry += norml * py;
    dl += 0.5 * norml * pnorm * (px*px + py*py);
  }
  local_2_global(pos, elem);
  return Status::success;
}

template <typename T>
Status::type pm_thinquad_pass(Pos<T> &pos, const Element &elem,
                              const Accelerator& accelerator) {

  return Status::passmethod_not_implemented;
}

template <typename T>
Status::type pm_thinsext_pass(Pos<T> &pos, const Element &elem,
                              const Accelerator& accelerator) {

  return Status::passmethod_not_implemented;
}


template <typename T>
Status::type pm_kickmap_pass(Pos<T> &pos, const Element &elem,
                             const Accelerator& accelerator) {

  if (elem.kicktable_idx < 0 or elem.kicktable_idx >= kicktable_list.size()) return Status::kicktable_not_defined;

  Status::type status = Status::success;

  double sl   = elem.length / float(elem.nr_steps);
  const double brho = get_magnetic_rigidity(accelerator.energy);

  global_2_local(pos, elem);
  if (elem.kicktable_idx < 0) {
    drift<T>(pos, sl);
  } else {
    for(unsigned int i=0; i<elem.nr_steps; ++i) {
      drift<T>(pos, sl / 2);
      Status::type status = kicktablethinkick(pos, elem.kicktable_idx, brho, elem.nr_steps, elem.rescale_kicks);
      if (status != Status::success) return status;
      drift<T>(pos, sl / 2);
    }
  }
  local_2_global(pos, elem);

  return status;
}

template <typename T>
Status::type pm_matrix_pass(Pos<T> &pos, const Element &elem,
                            const Accelerator& accelerator) {

  global_2_local(pos, elem);
  if (elem.length > 0) {
    double sl = elem.length / float(elem.nr_steps);
    double l1 = sl * DRIFT1;
    double l2 = sl * DRIFT2;
    double k1 = KICK1 / float(elem.nr_steps);
    double k2 = KICK2 / float(elem.nr_steps);
    Matrix mat1 = elem.matrix66;
    Matrix mat2 = elem.matrix66;
    multiply_transf_matrix66(mat1, k1);
    multiply_transf_matrix66(mat2, k2);
    for(unsigned int i=0; i<elem.nr_steps; ++i) {
      drift<T>(pos, l1);
      matthinkick<T>(pos, mat1);
      drift<T>(pos, l2);
      matthinkick<T>(pos, mat2);
      drift<T>(pos, l2);
      matthinkick<T>(pos, mat1);
      drift<T>(pos, l1);
    }
  } else {
    matthinkick<T>(pos, elem.matrix66);
  }
  local_2_global(pos, elem);
  return Status::success;
}

#endif
