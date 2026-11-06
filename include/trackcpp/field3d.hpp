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


inline void calc_mn_factors(
    double hori_ks, double hori_kx,
    CoefMatrix& hori_ky,
    std::vector<double>& mkx,
    std::vector<double>& mkx2,
    std::vector<double>& nks
)
{
    // (M, N) is the shape of the "hori_coefs_cos/sin" matrices
    size_t M = mkx.size(), N = nks.size();
    double _mkx, _nks, _nks2;
    for (size_t m = 0; m < M; ++m) {
        _mkx = (m + 1) * hori_kx;
        mkx[m] = _mkx;
        mkx2[m] = _mkx * _mkx;
    }
    for (size_t n = 0; n < N; ++n) {
        _nks = (n + 1) * hori_ks;
        nks[n] = _nks;
        _nks2 = _nks * _nks; 
        for (size_t m = 0; m < M; ++m) {
            hori_ky[m][n] = std::sqrt(mkx2[m] + _nks2);
        }
    }
}

inline void calc_sdependence(
    const CoefMatrix& hori_coefs_cos,
    const CoefMatrix& hori_coefs_sin,
    const std::vector<double>& nks,
    double s,
    CoefMatrix& sdependence
)
{
    size_t M = hori_coefs_cos.rows();
    size_t N = hori_coefs_cos.cols();
    double coss, sins, fac1, fac2;

    for (size_t n = 0; n < N; ++n) {
        coss = cos(nks[n] * s);
        sins = sin(nks[n] * s);
        for (size_t m = 0; m < M; ++m) {
            fac1 = hori_coefs_cos[m][n];
            fac2 = -hori_coefs_sin[m][n];
            sdependence[m][n] = fac1 * sins + fac2 * coss;
        }
    }
}

template <typename T>
void ay_inty_day_dx(
    double brho,
    const T& x, const T& y,
    const CoefMatrix& sdependence, const CoefMatrix& hori_ky,
    const std::vector<double>& mkx, const std::vector<double>& mkx2,
    const std::vector<double>& nks,
    T& ayn, T& idayn
)
{
    size_t M = sdependence.rows();
    size_t N = sdependence.cols();

    double _mkx, _mkx2, _nks, _hori_ky, ratio, _sdep;
    for (size_t m = 0; m < M; ++m) {
        _mkx = mkx[m];
        _mkx2 = mkx2[m];
        const T sinx = sin(_mkx * x);
        const T cosx = cos(_mkx * x);
        for (size_t n = 0; n < N; ++n) {
            _nks = nks[n];
            _hori_ky = hori_ky[m][n];
            _sdep = sdependence[m][n];
            ratio = _mkx / (_nks * _hori_ky);
            ayn += sinx * sinh(_hori_ky * y) * ratio * _sdep;
            ratio = _mkx2 / (_nks * _hori_ky * _hori_ky);
            idayn += cosx * (cosh(_hori_ky * y) - 1.0) * ratio * _sdep;
        }
    }
    ayn /= brho;
    idayn /= brho;
}

template <typename T>
void ax_intx_dax_dy(
    double brho,
    const T& x, const T& y,
    const CoefMatrix& sdependence, const CoefMatrix& hori_ky,
    const std::vector<double>& mkx,
    const std::vector<double>& nks,
    T& axn, T& idaxn
)
{
    size_t M = sdependence.rows();
    size_t N = sdependence.cols();

    double _mkx, _nks, _hori_ky, ratio, _sdep;
    for (size_t m = 0; m < M; ++m) {
        _mkx = mkx[m];
        const T sinx = sin(_mkx * x);
        const T cosx = cos(_mkx * x);
        for (size_t n = 0; n < N; ++n) {
            _nks = nks[n];
            _hori_ky = hori_ky[m][n];
            _sdep = sdependence[m][n];
            ratio = 1 / _nks;
            axn += cosx * cosh(_hori_ky * y) * ratio * _sdep;
            ratio = _hori_ky / (_nks * _mkx);
            idaxn += sinx * sinh(_hori_ky * y) * ratio * _sdep;
        }
    }
    axn /= brho;
    idaxn /= brho;
}

inline void prop_h1(double& s, double step) {
    s += step / 2.0;
}

template <typename T>
void prop_h2(Pos<T>& map, const T& pnorm, double step) {
    map.ry += 0.5*step*pnorm*map.py;
    map.dl += 0.25*step*pnorm*pnorm*(map.py * map.py);
}

template <typename T>
void prop_h3(Pos<T>& map, const T& pnorm, double step) {
    map.rx += step*pnorm*map.px;
    map.dl += 0.5*step*pnorm*pnorm*(map.px * map.px);
}

template <typename T>
void prop_ix(Pos<T>& map, double brho, int sign,
            const CoefMatrix& sdependence, const CoefMatrix& hori_ky,
            const std::vector<double>& mkx, const std::vector<double>& nks)
{
    T axn = 0.0, idaxn = 0.0;
    ax_intx_dax_dy(brho, map.rx, map.ry, sdependence, hori_ky, mkx, nks, axn, idaxn);
    map.px += axn * sign * -1.0;
    map.py += idaxn * sign * -1.0;
}

template <typename T>
void prop_iy(Pos<T>& map, double brho, int sign,
            const CoefMatrix& sdependence, const CoefMatrix& hori_ky,
            const std::vector<double>& mkx, const std::vector<double>& mkx2,
            const std::vector<double>& nks)
{
    T ayn = 0.0, idayn = 0.0;
    ay_inty_day_dx(brho, map.rx, map.ry, sdependence, hori_ky, mkx, mkx2, nks, ayn, idayn);
    map.px += idayn * sign * -1.0;
    map.py += ayn * sign * -1.0;
}

template <typename T>
void prop_step(
    Pos<T>& map, double& s,
    double brho, double step, 
    const T& pnorm,
    const CoefMatrix& hori_coefs_cos,
    const CoefMatrix& hori_coefs_sin,
    const CoefMatrix& hori_ky,
    CoefMatrix& sdependence,
    const std::vector<double>& mkx, 
    const std::vector<double>& mkx2,
    const std::vector<double>& nks
) {
    prop_h1(s, step);

    calc_sdependence(hori_coefs_cos, hori_coefs_sin, nks, s, sdependence);

    prop_iy(map, brho, +1, sdependence, hori_ky, mkx, mkx2, nks);
    prop_h2(map, pnorm, step);
    prop_iy(map, brho, -1, sdependence, hori_ky, mkx, mkx2, nks);
    prop_ix(map, brho, +1, sdependence, hori_ky, mkx, nks);
    prop_h3(map, pnorm, step);
    prop_ix(map, brho, -1, sdependence, hori_ky, mkx, nks);
    prop_iy(map, brho, +1, sdependence, hori_ky, mkx, mkx2, nks);
    prop_h2(map, pnorm, step);
    prop_iy(map, brho, -1, sdependence, hori_ky, mkx, mkx2, nks);
    prop_h1(s, step);
}

#endif
