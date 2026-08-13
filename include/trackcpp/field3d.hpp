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
// This code applies the symplectic integrator for 3D field tracking as described in the paper:
// "Explicit symplectic integrator for s-dependent static magnetic field"
//  DOI: https://doi.org/10.1103/PhysRevE.68.046502
#ifndef _FIELD3D_H
#define _FIELD3D_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include "auxiliary.h"
#include "pos.h"
#include "coefmatrix.h"


inline void calc_sdependence(
    const CoefMatrix& hori_coefs_cos,
    const CoefMatrix& hori_coefs_sin,
    const double& hori_ks,
    const double& s,
    CoefMatrix& sdependence
)
{
    size_t M = hori_coefs_cos.rows();
    size_t N = hori_coefs_cos.cols();
    for (int n = 1; n <= N; ++n) {
        const double nks = n * hori_ks;
        const double coss = cos(nks * s);
        const double sins = sin(nks * s);
        for (int m = 1; m <= M; ++m) {
            const double fac1 = hori_coefs_cos[m - 1][n - 1];
            const double fac2 = -hori_coefs_sin[m - 1][n - 1];
            sdependence[m-1][n-1] = fac1 * sins + fac2 * coss;
        }
    }
}

template <typename T>
void ay_inty_day_dx(
    const double& brho, const double& hori_kx, const double& hori_ks,
    const CoefMatrix& sdependence,
    const T& x, const T& y, const double& s,
    T& ayn, T& idayn
)
{
    size_t M = sdependence.rows();
    size_t N = sdependence.cols();

    std::vector<double> nks(N), nks2(N);
    for (int n = 1; n <= N; ++n) {
        nks[n-1] = n * hori_ks;
        nks2[n-1] = nks[n-1] * nks[n-1];
    }

    for (int m = 1; m <= M; ++m) {
        const double mkx = m * hori_kx;
        const double mkx2 = mkx * mkx;
        const T sinx = sin(mkx * x);
        const T cosx = cos(mkx * x);
        for (int n = 1; n <= N; ++n) {
            const double hori_ky = std::sqrt(mkx2 + nks2[n-1]);
            double ratio = mkx / (nks[n-1] * hori_ky);
            ayn += sinx * sinh(hori_ky * y) * ratio * sdependence[m-1][n-1];

            ratio = mkx2 / (nks[n-1] * hori_ky * hori_ky);
            idayn += cosx * (cosh(hori_ky * y) - 1.0) * ratio * sdependence[m-1][n-1];
        }
    }
    ayn /= brho;
    idayn /= brho;
}

template <typename T>
void ax_intx_dax_dy(
    const double& brho, const double& hori_kx, const double& hori_ks,
    const CoefMatrix& sdependence,
    const T& x, const T& y, const double& s,
    T& axn, T& idaxn
)
{
    size_t M = sdependence.rows();
    size_t N = sdependence.cols();

    std::vector<double> nks(N), nks2(N);
    for (int n = 1; n <= N; ++n) {
        nks[n-1] = n * hori_ks;
        nks2[n-1] = nks[n-1] * nks[n-1];
    }

    for (int m = 1; m <= M; ++m) {
        const double mkx = m * hori_kx;
        const double mkx2 = mkx * mkx;
        const T sinx = sin(mkx * x);
        const T cosx = cos(mkx * x);
        for (int n = 1; n <= N; ++n) {
            const double hori_ky = std::sqrt(mkx2 + nks2[n-1]);
            double ratio = 1 / nks[n-1];
            axn += cosx * cosh(hori_ky * y) * ratio * sdependence[m-1][n-1];

            ratio = hori_ky / (nks[n-1] * mkx);
            idaxn += sinx * sinh(hori_ky * y) * ratio * sdependence[m-1][n-1];
        }
    }
    axn /= brho;
    idaxn /= brho;
}

template <typename T>
void exp_h1_s(T& s, T step) {
    s += step / 2.0;
}

template <typename T>
void exp_h2_y(Pos<T>& map, const T& pnorm, double step) {
    map.ry += 0.5*step*pnorm*map.py;
}


template <typename T>
void exp_h2_z(Pos<T>& map, const T& pnorm, double step) {
    map.dl += 0.25*step*pnorm*pnorm*(map.py * map.py);
}


template <typename T>
void exp_h3_x(Pos<T>& map, const T& pnorm, double step) {
    map.rx += step*pnorm*map.px;
}

template <typename T>
void exp_h3_z(Pos<T>& map, const T& pnorm, double step) {
    map.dl += 0.5*step*pnorm*pnorm*(map.px * map.px);
}

template <typename T>
void prop_h1(Pos<T>& map, double& s, double step) {
    exp_h1_s(s, step);
}

template <typename T>
void prop_h2(Pos<T>& map, const T& pnorm, double step) {
    exp_h2_y(map, pnorm, step);
    exp_h2_z(map, pnorm, step);
}

template <typename T>
void prop_h3(Pos<T>& map, const T& pnorm, double step) {
    exp_h3_x(map, pnorm, step);
    exp_h3_z(map, pnorm, step);
}

template <typename T>
void prop_ix(const double& brho, const double& hori_kx, const double& hori_ks,
             const CoefMatrix& sdependence,
             Pos<T>& map, double s, int sign, double step) {
    T axn = 0.0;
    T idaxn = 0.0;
    ax_intx_dax_dy(
        brho, hori_kx, hori_ks, sdependence,
        map.rx, map.ry, s, axn, idaxn
    );
    map.px += axn * sign * -1.0;
    map.py += idaxn * sign * -1.0;
}

template <typename T>
void prop_iy(const double& brho, const double& hori_kx, const double& hori_ks,
             const CoefMatrix& sdependence,
             Pos<T>& map, double s, int sign, double step) {
    T ayn = 0.0;
    T idayn = 0.0;
    ay_inty_day_dx(
        brho, hori_kx, hori_ks, sdependence,
        map.rx, map.ry, s, ayn, idayn
    );
    map.px += idayn * sign * -1.0;
    map.py += ayn * sign * -1.0;
}

template <typename T>
void prop_step(const double& brho, const double& hori_kx, const double& hori_ks,
               const CoefMatrix& hori_coefs_cos,
               const CoefMatrix& hori_coefs_sin,
               Pos<T>& map, const T& pnorm, double& s, double step) {
    prop_h1(map, s, step);

    CoefMatrix sdependence(hori_coefs_cos.rows(), hori_coefs_cos.cols());
    calc_sdependence(hori_coefs_cos, hori_coefs_sin, hori_ks, s, sdependence);

    prop_iy(brho, hori_kx, hori_ks, sdependence, map, s, +1, step);
    prop_h2(map, pnorm, step);
    prop_iy(brho, hori_kx, hori_ks, sdependence, map, s, -1, step);
    prop_ix(brho, hori_kx, hori_ks, sdependence, map, s, +1, step);
    prop_h3(map, pnorm, step);
    prop_ix(brho, hori_kx, hori_ks, sdependence, map, s, -1, step);
    prop_iy(brho, hori_kx, hori_ks, sdependence, map, s, +1, step);
    prop_h2(map, pnorm, step);
    prop_iy(brho, hori_kx, hori_ks, sdependence, map, s, -1, step);
    prop_h1(map, s, step);
}

#endif
