#ifndef _PASS_METHOD_AT_H
#define _PASS_METHOD_AT_H

// TRACKC++
// ========
// Author:     Ximenes R. Resende
// email:      xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:  LNLS - Laboratorio Nacional de Luz Sincrotron
// Date:     Tue Dec 10 17:57:20 BRST 2013
//
// Obs: these are a c++ implementation of the AT passmethods.
//      apart from discrepancies between math library implementations and TWOPI,CGAMMA constants
//      these passmethods agree with AT passmethods up to machine 64-bit precision.

#define ATCOMPATIBLE 1

#include "passmethods.h"
#include "accelerator.h"
#include "elements.h"
#include "pos.h"
#include "auxiliary.h"

// constants for 4th-order symplectic integrator
#define DRIFT1 ( 0.6756035959798286638e00)
#define DRIFT2 (-0.1756035959798286639e00)
#define KICK1  ( 0.1351207191959657328e01)
#define KICK2  (-0.1702414383919314656e01)

#ifdef ATCOMPATIBLE
  #define TWOPI   6.28318530717959 // AT implementation of 2*PI...
  #define CGAMMA  8.846056192e-05  // AT implementation
#else
  //#define TWOPI   2*M_PI
  //static const double CGAMMA = 4*M_PI*electron_radius/pow(electron_rest_energy/electron_charge/1e9,3)/3;
#endif


template <typename T> inline T SQR(const T& X) { return X*X; }
template <typename T> inline T POW3(const T& X) { return X*X*X; }

template <typename T>
inline void drift(Pos<T>& pos, const double& length)
{
  T pnorm = 1 / (1 + pos.de);
  T norml = length * pnorm;
  pos.rx += norml * pos.px;
  pos.ry += norml * pos.py;
  pos.dl += 0.5 * norml * pnorm * (pos.px*pos.px + pos.py*pos.py);
}

//template <typename T>
//inline void drift(Pos<T> &pos, const double& length)
//{
//  drift(pos, length);
//}

template <typename T>
inline void calcpolykick(const Pos<T> &pos, const std::vector<double>& polynom_a,
                         const std::vector<double>& polynom_b,
                         T& real_sum, T& imag_sum)
{
  int n = polynom_b.size();
  real_sum = polynom_b[n-1];
  imag_sum = polynom_a[n-1];
  for(int i=n-2;i>=0;--i) {
    T real_sum_tmp = real_sum * pos.rx - imag_sum * pos.ry + polynom_b[i];
    imag_sum = imag_sum * pos.rx + real_sum * pos.ry + polynom_a[i];
    real_sum = real_sum_tmp;
  }
}

template <typename T>
void fastdrift(Pos<T> &pos, const T& norml)
{
  T dx = norml * pos.px;
  T dy = norml * pos.py;
  pos.rx += dx;
  pos.ry += dy;
  pos.dl += 0.5 * norml * (pos.px*pos.px + pos.py*pos.py) / (1 + pos.de);
}

template <typename T>
T b2_perp(const T& bx, const T& by, const T& rx, const T& px, const T& ry, const T& py, const double& irho = 0)
{
  // Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity
  T v_norm2 = 1 /(SQR(1+irho*rx) + SQR(px) + SQR(py));
  return((SQR(by*(1+irho*rx)) + SQR(bx*(1+irho*rx)) + SQR(bx*py - by*px))*v_norm2);
}

template <typename T>
Status::type kicktablethinkick(Pos<T>& pos, const Kicktable* kicktable,
                               const double& brho)
{
  T hkick, vkick;
  Status::type status = kicktable_getkicks(kicktable, pos.rx, pos.ry, hkick, vkick);
  pos.px += hkick / (brho * brho);
  pos.py += vkick / (brho * brho);
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
void strthinkick(Pos<T>& pos, const double& length,
                 const std::vector<double>& polynom_a,
                 const std::vector<double>& polynom_b,
                 const Accelerator& accelerator)
{
  T real_sum, imag_sum;
  calcpolykick<T>(pos, polynom_a, polynom_b, real_sum, imag_sum);
  if (accelerator.radiation_on) {
    T pnorm = 1 / (1 + pos.de);
    const T& rx = pos.rx;
    T  px = pos.px * pnorm;
    const T& ry = pos.ry;
    T  py = pos.py * pnorm;
    T b2p = b2_perp(imag_sum, real_sum, rx, px, ry, py, 0);
    double radiation_constant =
      CGAMMA*POW3(accelerator.energy/1e9)/(TWOPI); /*[m]/[GeV^3] M.Sands(4.1)*/
    pos.de -=
      radiation_constant*SQR(1+pos.de)*b2p*(1+(px*px + py*py)/2)*length;
    pnorm  = 1 / (1 + pos.de);
    pos.px = px / pnorm;
    pos.py = py / pnorm;
  }
  pos.px -= length * real_sum;
  pos.py += length * imag_sum;
}

template <typename T>
void bndthinkick(Pos<T>& pos, const double& length,
                 const std::vector<double>& polynom_a,
                 const std::vector<double>& polynom_b,
                 const double& irho,
                 const Accelerator& accelerator)
{
  T real_sum, imag_sum;
  calcpolykick<T>(pos, polynom_a, polynom_b, real_sum, imag_sum);
  T de = pos.de;
  if (accelerator.radiation_on) {
    T pnorm = 1 / (1 + pos.de);
    const T& rx = pos.rx;
    T  px = pos.px * pnorm;
    const T& ry = pos.ry;
    T  py = pos.py * pnorm;
    T b2p = b2_perp(imag_sum, real_sum + irho, rx, px, ry, py, irho);
    double radiation_constant =
      CGAMMA*POW3(accelerator.energy/1e9)/(TWOPI); /*[m]/[GeV^3] M.Sands(4.1)*/
    pos.de -=
      radiation_constant*SQR(1+pos.de)*b2p*(1+irho*rx + (px*px+py*py)/2)*length;
    pnorm = 1 / (1 + pos.de);
    pos.px = px / pnorm;
    pos.py = py / pnorm;
  }
  pos.px -= length * (real_sum - (de - pos.rx * irho) * irho);
  pos.py += length * imag_sum;
  pos.dl += length * irho * pos.rx;
}

template <typename T>
void edge_fringe(Pos<T>& pos, const double& inv_rho,
                 const double& edge_angle, const double& fint,
                const double& gap)
{
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
inline void translate_pos(Pos<T> &pos, const double* t)
{
  pos.rx += t[0]; pos.px += t[1];
  pos.ry += t[2]; pos.py += t[3];
  pos.de += t[4]; pos.dl += t[5];
}

template <typename T>
inline void rotate_pos(Pos<T> &pos, const double* R)
{
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
void global_2_local(Pos<T> &pos, const Element &elem)
{
  translate_pos(pos, elem.t_in);
  rotate_pos(pos, elem.r_in);
}

template <typename T>
void local_2_global(Pos<T> &pos, const Element &elem)
{
  rotate_pos(pos, elem.r_out);
  translate_pos(pos, elem.t_out);
}


template <typename T>
Status::type pm_identity_pass(Pos<T> &pos, const Element &elem,
                              const Accelerator& accelerator)
{
  return Status::success;
}

template <typename T>
Status::type pm_drift_pass(Pos<T> &pos, const Element &elem,
                           const Accelerator& accelerator)
{
  drift<T>(pos, elem.length);
  return Status::success;
}

template <typename T>
Status::type pm_str_mpole_symplectic4_pass(Pos<T> &pos, const Element &elem,
                                           const Accelerator& accelerator)
{
  global_2_local(pos, elem);
  double sl = elem.length / float(elem.nr_steps);
  double l1 = sl * DRIFT1;
  double l2 = sl * DRIFT2;
  double k1 = sl * KICK1;
  double k2 = sl * KICK2;
  const std::vector<double> &polynom_a = elem.polynom_a;
  const std::vector<double> &polynom_b = elem.polynom_b;
  for(unsigned int i=0; i<elem.nr_steps; ++i) {
    drift(pos, l1);
    strthinkick<T>(pos, k1, polynom_a, polynom_b, accelerator);
    drift(pos, l2);
    strthinkick<T>(pos, k2, polynom_a, polynom_b, accelerator);
    drift<T>(pos, l2);
    strthinkick<T>(pos, k1, polynom_a, polynom_b, accelerator);
    drift<T>(pos, l1);
  }
  local_2_global(pos, elem);
  return Status::success;
}

template <typename T>
Status::type pm_bnd_mpole_symplectic4_pass(Pos<T> &pos, const Element &elem,
                                           const Accelerator& accelerator)
{
  double sl = elem.length / float(elem.nr_steps);
  double l1 = sl * DRIFT1;
  double l2 = sl * DRIFT2;
  double k1 = sl * KICK1;
  double k2 = sl * KICK2;
  double irho = elem.angle / elem.length;
  const std::vector<double> &polynom_a = elem.polynom_a;
  const std::vector<double> &polynom_b = elem.polynom_b;

  global_2_local(pos, elem);
  edge_fringe(pos, irho, elem.angle_in, elem.fint_in, elem.gap);
  for(unsigned int i=0; i<elem.nr_steps; ++i) {
    drift<T>(pos, l1);
    bndthinkick<T>(pos, k1, polynom_a, polynom_b, irho, accelerator);
    drift<T>(pos, l2);
    bndthinkick<T>(pos, k2, polynom_a, polynom_b, irho, accelerator);
    drift<T>(pos, l2);
    bndthinkick<T>(pos, k1, polynom_a, polynom_b, irho, accelerator);
    drift<T>(pos, l1);
  }
  edge_fringe(pos, irho, elem.angle_out, elem.fint_out, elem.gap);
  local_2_global(pos, elem);

  return Status::success;
}


template <typename T>
Status::type pm_corrector_pass(Pos<T> &pos, const Element &elem,
                               const Accelerator& accelerator)
{
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
Status::type pm_cavity_pass(Pos<T> &pos, const Element &elem,
                            const Accelerator& accelerator)
{
  if (not accelerator.cavity_on) return pm_drift_pass(pos, elem, accelerator);

  global_2_local(pos, elem);
  double nv = elem.voltage / accelerator.energy;
  if (elem.length == 0) {
    T &de = pos.de, &dl = pos.dl;
    de +=  -nv * sin(TWOPI*elem.frequency * dl/ light_speed);
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
    de += -nv * sin(TWOPI*elem.frequency*dl/light_speed);
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
                              const Accelerator& accelerator)
{
  return Status::passmethod_not_implemented;
}

template <typename T>
Status::type pm_thinsext_pass(Pos<T> &pos, const Element &elem,
                              const Accelerator& accelerator)
{
  return Status::passmethod_not_implemented;
}


template <typename T>
Status::type pm_kicktable_pass(Pos<T> &pos, const Element &elem,
                               const Accelerator& accelerator)
{
  if (elem.kicktable == nullptr) return Status::kicktable_not_defined;

  Status::type status = Status::success;

  double sl   = elem.length / float(elem.nr_steps);
  double brho = get_magnetic_rigidity(accelerator.energy);

  global_2_local(pos, elem);
  if (elem.kicktable == nullptr) {
    drift<T>(pos, sl);
  } else {
    for(unsigned int i=0; i<elem.nr_steps; ++i) {
      drift<T>(pos, sl / 2);
      Status::type status = kicktablethinkick(pos, elem.kicktable, brho);
      if (status != Status::success) return status;
      drift<T>(pos, sl / 2);
    }
  }
  local_2_global(pos, elem);

  return status;
}

#undef DRIFT1
#undef DRIFT2
#undef KICK1
#undef KICK2
#undef SQR
#undef TWOPI
#undef CGAMMA

#endif
