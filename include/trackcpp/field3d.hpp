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

#include <trackcpp/auxiliary.h>
#include <string>
#include <fstream>
#include <cmath>


template <typename T>
T ay(const double& brho, const double& hori_kx, const double& hori_ks, 
              const std::vector<std::vector<double>>& hori_coefs_cos,
              const std::vector<std::vector<double>>& hori_coefs_sin,
              const T& x, const T& y, const double& s) 
{
    T ay_ = 0.0;
    int M = hori_coefs_cos.size();
    int N = hori_coefs_cos[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double hori_ky = std::sqrt(std::pow(m * hori_kx, 2) + std::pow(n * hori_ks, 2));
            double fac1 = hori_coefs_cos[m - 1][n - 1] * (m * hori_kx) / (n * hori_ks * hori_ky);
            double fac2 = -hori_coefs_sin[m - 1][n - 1] * (m * hori_kx) / (n * hori_ks * hori_ky);
            ay_ += fac1 * sin(m * hori_kx * x) * sinh(hori_ky * y) * std::sin(n * hori_ks * s);
            ay_ += fac2 * sin(m * hori_kx * x) * sinh(hori_ky * y) * std::cos(n * hori_ks * s);
        }
    }

    return 1 * ay_/brho;
}

template <typename T>
T ax(const double& brho, const double& hori_kx, const double& hori_ks, 
              const std::vector<std::vector<double>>& hori_coefs_cos,
              const std::vector<std::vector<double>>& hori_coefs_sin,
              const T& x, const T& y, const double& s)
{
    T ax_ = 0.0;
    const int M = hori_coefs_cos.size();
    const int N = hori_coefs_cos[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double hori_ky = std::sqrt(std::pow(m * hori_kx, 2) + std::pow(n * hori_ks, 2));
            double fac1 = hori_coefs_cos[m - 1][n - 1] / (n * hori_ks);
            double fac2 = -hori_coefs_sin[m - 1][n - 1] / (n * hori_ks);
            ax_ += fac1 * cos(m * hori_kx * x) * cosh(hori_ky * y) * std::sin(n * hori_ks * s);
            ax_ += fac2 * cos(m * hori_kx * x) * cosh(hori_ky * y) * std::cos(n * hori_ks * s);
        }
    }

    return 1 * ax_/brho;
}

template <typename T>
T inty_day_dx(const double& brho, const double& hori_kx, const double& hori_ks, 
              const std::vector<std::vector<double>>& hori_coefs_cos,
              const std::vector<std::vector<double>>& hori_coefs_sin,
              const T& x, const T& y, const double& s)
{
    T day_dx = 0.0;
    int M = hori_coefs_cos.size();
    int N = hori_coefs_cos[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double hori_ky = std::sqrt(std::pow(m * hori_kx, 2) + std::pow(n * hori_ks, 2));
            double fac1 = hori_coefs_cos[m - 1][n - 1] * std::pow(m * hori_kx, 2) / (n * hori_ks * std::pow(hori_ky, 2));
            double fac2 = -hori_coefs_sin[m - 1][n - 1] * std::pow(m * hori_kx, 2) / (n * hori_ks * std::pow(hori_ky, 2));
            day_dx += fac1 * cos(m * hori_kx * x) * std::sin(n * hori_ks * s) * (cosh(hori_ky * y) - 1.0);
            day_dx += fac2 * cos(m * hori_kx * x) * std::cos(n * hori_ks * s) * (cosh(hori_ky * y) - 1.0);
        }
    }

    return 1 * day_dx/brho;
}

template <typename T>
T intx_dax_dy(const double& brho, const double& hori_kx, const double& hori_ks, 
              const std::vector<std::vector<double>>& hori_coefs_cos,
              const std::vector<std::vector<double>>& hori_coefs_sin,
              const T& x, const T& y, const double& s)
{
    T dax_dy = 0.0;
    int M = hori_coefs_cos.size();
    int N = hori_coefs_cos[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double hori_ky = std::sqrt(std::pow(m * hori_kx, 2) + std::pow(n * hori_ks, 2));
            double fac1 = hori_coefs_cos[m - 1][n - 1] * hori_ky / (n * hori_ks * m * hori_kx);
            double fac2 = -hori_coefs_sin[m - 1][n - 1] * hori_ky / (n * hori_ks * m * hori_kx);
            dax_dy += fac1 * sin(m * hori_kx * x) * std::sin(n * hori_ks * s) * sinh(hori_ky * y);
            dax_dy += fac2 * sin(m * hori_kx * x) * std::cos(n * hori_ks * s) * sinh(hori_ky * y);
        }
    }

    return 1* dax_dy/brho;
}

template <typename T>
void exp_h1_s(T& s, T step) {
    s += step / 2.0;
}

template <typename T>
void exp_iy_px(const double& brho, const double& hori_kx, const double& hori_ks, 
               const std::vector<std::vector<double>>& hori_coefs_cos,
               const std::vector<std::vector<double>>& hori_coefs_sin,
               Pos<T>& map, double s, int sign, double step) {
    T factor = inty_day_dx(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map.rx, map.ry, s) * sign * -1.0;
    map.px += factor;
}


template <typename T>
void exp_iy_py(const double& brho, const double& hori_kx, const double& hori_ks, 
               const std::vector<std::vector<double>>& hori_coefs_cos,
               const std::vector<std::vector<double>>& hori_coefs_sin,
               Pos<T>& map, double s, int sign, double step) {
    T factor = ay(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map.rx, map.ry, s) * sign * -1.0;
    map.py += factor;
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
void exp_ix_px(const double& brho, const double& hori_kx, const double& hori_ks, 
               const std::vector<std::vector<double>>& hori_coefs_cos,
               const std::vector<std::vector<double>>& hori_coefs_sin,
               Pos<T>& map, double s, int sign, double step) {
    T factor = ax(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map.rx, map.ry, s) * sign * -1.0;
    map.px += factor;
}


template <typename T>
void exp_ix_py(const double& brho, const double& hori_kx, const double& hori_ks, 
               const std::vector<std::vector<double>>& hori_coefs_cos,
               const std::vector<std::vector<double>>& hori_coefs_sin,
               Pos<T>& map, double s, int sign, double step) {
    T factor = intx_dax_dy(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map.rx, map.ry, s) * sign * -1.0;
    map.py += factor;
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
             const std::vector<std::vector<double>>& hori_coefs_cos,
             const std::vector<std::vector<double>>& hori_coefs_sin,
             Pos<T>& map, double s, int sign, double step) {
    exp_ix_px(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, sign, step);
    exp_ix_py(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, sign, step);

}

template <typename T>
void prop_iy(const double& brho, const double& hori_kx, const double& hori_ks, 
             const std::vector<std::vector<double>>& hori_coefs_cos,
             const std::vector<std::vector<double>>& hori_coefs_sin,
             Pos<T>& map, double s, int sign, double step) {
    exp_iy_px(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, sign, step);
    exp_iy_py(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, sign, step);

}

template <typename T>
void prop_step(const double& brho, const double& hori_kx, const double& hori_ks,
               const std::vector<std::vector<double>>& hori_coefs_cos,
               const std::vector<std::vector<double>>& hori_coefs_sin,
               Pos<T>& map, const T& pnorm, double& s, double step) {
    prop_h1(map, s, step);
    prop_iy(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, +1, step);
    prop_h2(map, pnorm, step);
    prop_iy(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, -1, step);
    prop_ix(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, +1, step);
    prop_h3(map, pnorm, step);
    prop_ix(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, -1, step);
    prop_iy(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, +1, step);
    prop_h2(map, pnorm, step);
    prop_iy(brho, hori_kx, hori_ks, hori_coefs_cos, hori_coefs_sin, map, s, -1, step);
    prop_h1(map, s, step);
}
