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

#include <trackcpp/diffusion_matrix.h>

// NOTE: I don't understand this pass method.
// There is a discrepancy between edge_fringe from the dipole
// pass method and this one.
// I decided to keep this unchanged so far, to reproce AT results.
// Besides, I don't understand the linearized matriz used to propagate B.
void edge_fringe_b(Pos<double>& pos, Matrix& bdiff, const double& irho,
                   const double& edge_angle, const double& fint,
                   const double& gap) {

  const double &rx = pos.rx, &ry = pos.ry, &de = pos.de;

  const double psi = irho * gap * fint / std::cos(edge_angle) *
    (1 + std::sin(edge_angle) * std::sin(edge_angle));
  const double fx = irho * std::tan(edge_angle);  // 1/(1+de) ???
  const double fy = irho * std::tan(edge_angle - psi/(1+de));  // 1/(1+de) ???
  const double fac = ry*(irho*irho + fy*fy)*psi/(1+de)/(1+de)/irho;

  //  Propagate B
  for(unsigned int m=0; m<6; ++m){
    bdiff[1][m] += fx * bdiff[0][m];
    bdiff[3][m] -= fy * bdiff[2][m];
  }
  if(fint > 0 && gap > 0)
    for(unsigned int m=0; m<6; ++m)
        bdiff[3][m] -= bdiff[4][m] * fac;

  for(unsigned int m=0; m<6; ++m) {
    bdiff[m][1] += fx * bdiff[m][0];
    bdiff[m][3] -= fy * bdiff[m][2];
  }

  if(fint > 0 && gap > 0)
    for(unsigned int m=0; m<6; m++)
      bdiff[m][3] -= bdiff[m][4] * fac;

  pos.px += rx * fx;
  pos.py -= ry * fy;
}


void drift_propagate(Pos<double>& pos, Matrix& bdiff, const double L){

  Matrix drift (6);
  drift.eye(1);

  const double ptot = 1 + pos.de;
  const double lop = L/ptot;

  // This matrix is the linearization of the drift map.
  // If you derive the map in relation to the coordinates
  // you will get this matrix
  drift[0][1] = lop;
  drift[2][3] = lop;
  drift[0][4] = -lop*pos.px/ptot;
  drift[2][4] = -lop*pos.py/ptot;
  drift[5][1] = lop*pos.px/ptot;
  drift[5][3] = lop*pos.py/ptot;
  drift[5][4] = -lop * (pos.px*pos.px + pos.py*pos.py)/(ptot*ptot);

  pos.rx += lop * pos.px;
  pos.ry += lop * pos.py;
  pos.dl += lop / ptot * (pos.px*pos.px + pos.py*pos.py) / 2;

  bdiff.sandwichme_with_matrix(drift);
}


void thinkick_symplectic(const Pos<double>& pos, Matrix& bdiff,
                         const Vector& pol_a, const Vector& pol_b,
                         const double frac, const double irho,
                         const int max_order) {

  double im_sum = max_order * pol_a[max_order];
  double re_sum = max_order * pol_b[max_order];

  // Recursively calculate the derivatives
  //   ReSumN = (irho/B0)*Re(d(By + iBx)/dx)
  //   ImSumN = (irho/B0)*Im(d(By + iBx)/dy)
  for(unsigned int i=max_order-1; i>0; --i){
    const double& re_tmp = re_sum*pos.rx - im_sum*pos.ry + i*pol_b[i];
    im_sum = im_sum*pos.rx + re_sum*pos.ry + i*pol_a[i];
    re_sum = re_tmp;
  }

  Matrix m66 (6);
  m66.eye(1);
  m66[1][0] = -frac*re_sum;
  m66[1][2] = frac*im_sum;
  m66[3][0] = frac*im_sum;
  m66[3][2] = frac*re_sum;
  m66[1][4] = frac*irho;
  m66[1][0] += -frac*irho*irho;
  m66[5][0] = frac*irho;
  bdiff.sandwichme_with_matrix(m66);
}


void thinkick_rad(Pos<double>& pos, Matrix& bdiff, const Vector& pol_a,
                const Vector& pol_b, const double frac, const double irho,
                const int max_order, const double energy, const bool radon){

  const double CRAD = CGAMMA * energy * energy * energy / (TWOPI*1e27);
  const double p_norm = (1+pos.de);
  const double p_norm2 = p_norm*p_norm;
  double im_sum = pol_a[max_order];
  double re_sum = pol_b[max_order];
  // save a copy of the initial value of dp/p
  const double dp_0 = pos.de;
  const double& xl = pos.px / p_norm;
  const double& yl = pos.py / p_norm;

  // recursively calculate the local transverse magnetic field
  // re_sum = irho*By/B0
  // im_sum = irho*Bx/B0
  for(int i=max_order-1; i>=0; --i){
    const double& re_tmp = re_sum*pos.rx - im_sum*pos.ry + pol_b[i];
    im_sum = im_sum*pos.rx + re_sum*pos.ry + pol_a[i];
    re_sum = re_tmp;
  }

  // calculate |B x n|^3 - the third power of the B field component
  // orthogonal to the normalized velocity vector n
  const double& b2p = b2_perp(
    im_sum, re_sum+irho, pos.rx, xl, pos.ry, yl, irho);
  const double& b3p = b2p*std::sqrt(b2p);

  const double& gamma = energy/M0C2;
  const double& gamma2 = gamma*gamma;
  const double& gamma4 = gamma2*gamma2;
  const double& bb_ = CU*CER*LAMBDABAR*gamma*gamma4*frac*b3p*p_norm2*p_norm2 *
                      (1 + pos.rx*irho + (xl*xl + yl*yl)/2);

  // Add rad kick to bdiff
  bdiff[1][1] += bb_*xl*xl;
  bdiff[1][3] += bb_*xl*yl;
  bdiff[1][4] += bb_*xl;
  bdiff[3][3] += bb_* yl*yl;
  bdiff[3][4] += bb_*yl;
  bdiff[3][1] += bb_*xl*yl;
  bdiff[4][1] += bb_*xl;
  bdiff[4][3] += bb_*yl;
  bdiff[4][4] += bb_;

  if (radon){
    // Update trajectory
    pos.de -= CRAD*frac*p_norm2*b2p*(1 + pos.rx*irho + (xl*xl + yl*yl)/2);

    // recalculate momenta from angles after losing energy to radiation
    const double& p_after = 1 + pos.de;
    pos.px = xl * p_after;
    pos.py = yl * p_after;
  }

  pos.px -= frac*(re_sum - (dp_0-pos.rx*irho)*irho);
  pos.py += frac*im_sum;
  pos.dl += frac*irho*pos.rx;
}

// Find Ohmi's diffusion matrix b_diff of a thick multipole.
void propagate_b_diff(const Element& ele, const Pos<double>& pos,
                      const double energy, const bool radon, Matrix& bdiff) {

  const double irho = ele.angle / ele.length;
  const std::vector<double>& pola = ele.polynom_a;
  const std::vector<double>& polb = ele.polynom_b;
  const int mord = std::min(pola.size(), polb.size()) - 1;
  // const int mord = std::max(std::min(pola.size(), polb.size()) - 1, 1);

  // 4-th order symplectic integrator constants
  const double frac = ele.length / float(ele.nr_steps);
  const double len1 = frac * DRIFT1;
  const double len2 = frac * DRIFT2;
  const double kick1 = frac * KICK1;
  const double kick2 = frac * KICK2;

  Pos<double> rin = pos;
  // Transform rin to a local coordinate system of the element
  global_2_local(rin, ele);

  // Propagate rin and b_diff through the entrance edge
  edge_fringe_b(rin, bdiff, irho, ele.angle_in, ele.fint_in, ele.gap);

  // Propagate rin and b_diff through a 4-th order integrator
  for(unsigned int i=0; i<ele.nr_steps; ++i) {
    drift_propagate(rin, bdiff, len1);

    thinkick_symplectic(rin, bdiff, pola, polb, kick1, irho, mord);
    thinkick_rad(rin, bdiff, pola, polb, kick1, irho, mord, energy, radon);

    drift_propagate(rin, bdiff, len2);

    thinkick_symplectic(rin, bdiff, pola, polb, kick2, irho, mord);
    thinkick_rad(rin, bdiff, pola, polb, kick2, irho, mord, energy, radon);

    drift_propagate(rin, bdiff, len2);

    thinkick_symplectic(rin, bdiff, pola, polb, kick1, irho, mord);
    thinkick_rad(rin, bdiff, pola, polb, kick1, irho, mord, energy, radon);

    drift_propagate(rin, bdiff, len1);
  }

  edge_fringe_b(rin, bdiff, irho, ele.angle_out, ele.fint_out, ele.gap);

  // Transform orbit to the global coordinate system
  local_2_global(rin, ele);
}


Status::type track_diffusionmatrix (const Accelerator& accelerator,
                            const Pos<double>& fixed_point,
                            const std::vector<Matrix>& tm,
                            std::vector<Matrix>& bmat) {

  Status::type status  = Status::success;
  const std::vector<Element>& lattice = accelerator.lattice;
  const double energy = accelerator.energy;
  const bool radon = accelerator.radiation_on;

  Pos<double> fp = fixed_point;

  const int nr_elements  = lattice.size();
	std::vector<bool> indcs;

  bmat.clear(); bmat.reserve(nr_elements+1);
  Matrix bdiff (6);
  bmat.push_back(bdiff);
  for(unsigned int i=0; i<nr_elements; ++i) {
    const Element& ele = lattice[i];
    if (ele.pass_method == PassMethod::pm_str_mpole_symplectic4_pass)
      propagate_b_diff(ele, fp, energy, radon, bdiff);
    else if (ele.pass_method == PassMethod::pm_bnd_mpole_symplectic4_pass)
      propagate_b_diff(ele, fp, energy, radon, bdiff);
    else
      bdiff.sandwichme_with_matrix(tm[i]);
    bmat.push_back(bdiff);
    // track through element
    track_elementpass (ele, fp, accelerator);
  }
  return status;
}
